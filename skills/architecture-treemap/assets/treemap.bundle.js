/**
 * treemap.bundle.js — architecture-treemap interactive renderer
 *
 * Expects a global COMPONENTS_DATA object (injected by render_treemap.py).
 * Requires D3 v7 to be loaded before this script.
 *
 * Three top-level views:
 *   Logical  — treemap sized by size_estimate_loc, coloured by classification
 *   Physical — treemap sized by size_loc, coloured by primary logical-owner classification
 *   Graph    — d3.forceSimulation with three-cluster layout (core left, seam centre, removable right)
 *
 * Edge overlay (Logical and Graph views):
 *   Click a cell/node → edges light up (outgoing solid, incoming dashed by direction).
 *   Shift-click adds more edges to the overlay.
 *   Background click clears.
 *   Header toggle "Show all edges (faint)" — off by default.
 *
 * Evidence class styling (orthogonal to direction):
 *   evidence_class="static"         → solid stroke (no dash)
 *   evidence_class="audit-asserted" → medium-dash stroke (4 4)
 *   absent (direct-call default: solid, others: audit-asserted behaviour)
 *
 *   Direction encoding uses OPACITY + COLOR-LIGHTNESS:
 *   outgoing edges are drawn at full colour + opacity.
 *   Incoming edges are drawn lighter (mixed toward white) + 70% opacity.
 *   This keeps two orthogonal channels (evidence_class=dasharray, direction=opacity/tone)
 *   visually distinct without fighting each other.
 *
 * Pill-accumulation bug from prototype fixed:
 *   The prototype called g.each() after data join, accumulating SVG elements on
 *   re-render.  This renderer calls svg.selectAll('*').remove() at the top of
 *   each render function and rebuilds from scratch.  No accumulation possible.
 *
 * Metric badges (top-left of each tile):
 *   Up to 3 badges drawn for whichever of these metrics are present:
 *   high fan-out (fan_out >= 5), high cyclomatic (cyclomatic >= 10),
 *   low test_ratio (test_ratio < 0.3), high churn (churn_90d >= 10).
 *   Full metrics breakdown shown in the side panel.
 */

// ── Constants ─────────────────────────────────────────────────────────────────

const CLS_BASE = {
  core:      '#2a7a2a',
  seam:      '#c87800',
  removable: '#a02020',
};

// Edge type stroke colours.
// The graph view (dark #111 background) uses ETYPE_COLOR_GRAPH for better legibility;
// the treemap views (slightly lighter background) use ETYPE_COLOR for the overlay.
const ETYPE_COLOR = {
  'direct-call':          '#888',
  'shared-state':         '#3a86c0',
  'background-knowledge': '#8a5cb0',
};

// Lightened palette for the graph view's dark background.
const ETYPE_COLOR_GRAPH = {
  'direct-call':          '#aab',
  'shared-state':         '#60b0e8',
  'background-knowledge': '#b080d8',
};

// Metric badge thresholds and labels
const METRIC_BADGES = [
  { field: 'fan_out',    test: v => v >= 5,   label: 'hi fan-out', color: '#6060a0' },
  { field: 'cyclomatic', test: v => v >= 10,  label: 'hi CC',      color: '#a06020' },
  { field: 'test_ratio', test: v => v < 0.3,  label: 'low tests',  color: '#a04040' },
  { field: 'churn_90d',  test: v => v >= 10,  label: 'hi churn',   color: '#405080' },
];

// Force-graph layout zones (logical x-centre targets, normalised 0-1)
const CLUSTER_X = { core: 0.2, seam: 0.5, removable: 0.8 };
const CLUSTER_Y = 0.5;

// ── Data maps (built once at load time) ───────────────────────────────────────

const logicalById  = new Map(
  (COMPONENTS_DATA.logical_components || []).map(c => [c.id, c])
);
const physicalById = new Map(
  (COMPONENTS_DATA.physical_components || []).map(c => [c.id, c])
);
const pruneByComp  = new Map(
  (COMPONENTS_DATA.prune_candidates || []).map(p => [p.logical_component, p])
);
const findingIndex = COMPONENTS_DATA.finding_index || {};

// ── State ─────────────────────────────────────────────────────────────────────

let currentView       = 'logical';
let selectedNodeIds   = new Set();   // ids of nodes whose edges are shown
let showAllEdges      = false;       // "Show all edges (faint)" toggle
let graphSimulation   = null;        // d3 force simulation (Graph view)

// ── DOM refs ─────────────────────────────────────────────────────────────────

const svg       = d3.select('#treemap');
const tip       = document.getElementById('tooltip');
const panelTitle = document.getElementById('panel-title');
const panelBody  = document.getElementById('panel-body');

// ── Helpers ───────────────────────────────────────────────────────────────────

function getViewportSize() {
  const area = document.getElementById('treemap-area');
  return { w: area.clientWidth, h: area.clientHeight };
}

function clsColor(cls, loc, maxLoc) {
  // Interpolate from a near-neutral dark toward the classification base colour.
  // t=1 at max LOC, t=0.25 at zero LOC, so even tiny components stay legible.
  const t = 0.25 + 0.75 * (loc / Math.max(maxLoc, 1));
  const base = d3.color(CLS_BASE[cls] || '#444');
  return d3.interpolateRgb(d3.color('#2a2a2a'), base)(t);
}

function primaryClassForPhysical(pc) {
  const owners = pc.logical_owners || [];
  if (owners.length === 0) return 'core';
  // Priority: removable > seam > core (most architecturally interesting wins)
  const rank = { removable: 3, seam: 2, core: 1 };
  let best = 'core', bestRank = 0;
  for (const id of owners) {
    const lc = logicalById.get(id);
    if (!lc) continue;
    const r = rank[lc.classification] || 0;
    if (r > bestRank) { bestRank = r; best = lc.classification; }
  }
  return best;
}

function wrapText(str, maxChars) {
  if (!str) return [''];
  const words = str.replace(/_/g, ' ').split(' ');
  const lines = [];
  let cur = '';
  for (const w of words) {
    if ((cur + ' ' + w).trim().length > maxChars) {
      if (cur) lines.push(cur);
      cur = w;
    } else {
      cur = (cur + ' ' + w).trim();
    }
  }
  if (cur) lines.push(cur);
  return lines.length ? lines : [str.substring(0, maxChars)];
}

function badgesForMetrics(metrics) {
  if (!metrics) return [];
  const badges = [];
  for (const def of METRIC_BADGES) {
    if (metrics[def.field] !== undefined && def.test(metrics[def.field])) {
      badges.push(def);
      if (badges.length >= 3) break;
    }
  }
  return badges;
}

// Draw up to 3 metric badges in the top-left of a cell (SVG g element).
function drawMetricBadges(gEl, metrics, cellW, cellH) {
  if (!metrics || cellW < 40 || cellH < 16) return;
  const badges = badgesForMetrics(metrics);
  if (badges.length === 0) return;

  const badgeH = 9;
  const pad = 2;
  let bx = 3;
  const by = 3;

  for (const def of badges) {
    const label = def.label;
    const bw = label.length * 5 + 6;
    if (bx + bw > cellW - 3) break;
    gEl.append('rect')
      .attr('class', 'metric-badge')
      .attr('x', bx).attr('y', by)
      .attr('width', bw).attr('height', badgeH)
      .attr('rx', 2)
      .attr('fill', def.color)
      .attr('opacity', 0.85);
    gEl.append('text')
      .attr('class', 'metric-badge')
      .attr('x', bx + 3).attr('y', by + 7)
      .attr('font-size', 6)
      .attr('fill', '#fff')
      .text(label);
    bx += bw + pad;
  }
}

// ── Edge overlay helpers ──────────────────────────────────────────────────────

// Compute the dash-array for an edge based on evidence_class (orthogonal channel).
// static → no dash (solid). audit-asserted → 4 4 dash.
// Absent: direct-call defaults to solid, everything else to audit-asserted.
function edgeDashArray(edge) {
  const ec = edge.evidence_class;
  if (ec === 'static') return null;
  if (ec === 'audit-asserted') return '4 4';
  // fallback default
  return edge.type === 'direct-call' ? null : '4 4';
}

// Return edges relevant to the current view (both endpoints must exist in view).
function edgesForView(viewIds) {
  const viewIdSet = new Set(viewIds);
  return (COMPONENTS_DATA.edges || []).filter(
    e => viewIdSet.has(e.from) && viewIdSet.has(e.to)
  );
}

// ── Treemap — shared rendering helpers ────────────────────────────────────────

function buildTreemapHierarchy(groups, nodeKey) {
  const { w, h } = getViewportSize();
  const root = d3.hierarchy({ id: 'root', children: groups })
    .sum(d => d.value || 0)
    .sort((a, b) => (b.value || 0) - (a.value || 0));
  d3.treemap().size([w, h]).paddingOuter(4).paddingInner(2).paddingTop(18)(root);
  return root;
}

// Draw classification group labels (CORE, SEAM, REMOVABLE).
function drawGroupLabels(svg, groupData) {
  const groupNodes = svg.selectAll('g.group')
    .data(groupData)
    .join('g').attr('class', 'group');

  groupNodes.append('rect')
    .attr('x', d => d.x0).attr('y', d => d.y0)
    .attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0)
    .attr('fill', 'none')
    .attr('stroke', d => CLS_BASE[d.data.id] || '#444')
    .attr('stroke-width', 1);

  groupNodes.append('text')
    .attr('x', d => d.x0 + 4).attr('y', d => d.y0 + 13)
    .attr('fill', d => CLS_BASE[d.data.id] || '#888')
    .attr('font-size', 10)
    .attr('font-weight', 700)
    .text(d => d.data.name ? d.data.name.toUpperCase() : '');

  return groupNodes;
}

// ── SVG edge overlay (shared by Logical and Physical views) ──────────────────

// Per-view node-centre maps, populated during each render call.
let logicalNodeCenters  = new Map();
let physicalNodeCenters = new Map();

// Return the node-centre map for the active view.
function nodeCentersForView(view) {
  return view === 'physical' ? physicalNodeCenters : logicalNodeCenters;
}

// Draw (or refresh) the Bezier edge overlay for a treemap view.
//   viewIds    — array of component ids that exist in this view
//   nodeCenters — Map<id, {cx, cy}> built during the render pass
function buildEdgeOverlaySvg(viewIds, nodeCenters) {
  // Remove previous overlay layer
  svg.select('g.edge-overlay').remove();

  const activeEdges = edgesForView(viewIds);
  const shownEdges  = showAllEdges
    ? activeEdges
    : activeEdges.filter(e => selectedNodeIds.has(e.from) || selectedNodeIds.has(e.to));

  if (shownEdges.length === 0) return;

  const overlay = svg.append('g').attr('class', 'edge-overlay').attr('pointer-events', 'none');

  // Arrowhead markers in <defs>: one per edge-type colour, outgoing and incoming.
  // Outgoing: solid filled arrowhead on marker-end.
  // Incoming: open chevron on marker-start (distinguishable hue-shift via separate marker).
  const defs = overlay.append('defs');
  const ARROW_TYPES = [
    { id: 'direct-call',          color: '#888' },
    { id: 'shared-state',         color: '#3a86c0' },
    { id: 'background-knowledge', color: '#8a5cb0' },
  ];
  for (const t of ARROW_TYPES) {
    // marker-end arrowhead (outgoing): solid triangle, points along path direction
    defs.append('marker')
      .attr('id', `arrow-out-${t.id}`)
      .attr('viewBox', '0 -4 10 8')
      .attr('refX', 9)
      .attr('refY', 0)
      .attr('markerWidth', 6)
      .attr('markerHeight', 6)
      .attr('orient', 'auto')
      .append('path')
        .attr('d', 'M0,-4L10,0L0,4Z')
        .attr('fill', t.color)
        .attr('opacity', 0.9);

    // marker-start arrowhead (incoming): open chevron, lighter fill
    const inColor = d3.interpolateRgb(d3.color(t.color), d3.color('#ddd'))(0.45);
    defs.append('marker')
      .attr('id', `arrow-in-${t.id}`)
      .attr('viewBox', '0 -4 10 8')
      .attr('refX', 1)
      .attr('refY', 0)
      .attr('markerWidth', 6)
      .attr('markerHeight', 6)
      .attr('orient', 'auto-start-reverse')
      .append('path')
        .attr('d', 'M10,-4L0,0L10,4')
        .attr('fill', 'none')
        .attr('stroke', inColor)
        .attr('stroke-width', 1.5);
  }

  for (const edge of shownEdges) {
    const src = nodeCenters.get(edge.from);
    const snk = nodeCenters.get(edge.to);
    if (!src || !snk) continue;

    // Direction cue: is this outgoing from a selected node, or incoming?
    const isOutgoing = selectedNodeIds.has(edge.from);
    const { w, h } = getViewportSize();
    const cx = w / 2;
    const cy = h / 2;

    // Cubic-bezier: pull control points toward viewport centre for a gentle arc.
    const cp1x = src.cx + (cx - src.cx) * 0.4;
    const cp1y = src.cy + (cy - src.cy) * 0.4;
    const cp2x = snk.cx + (cx - snk.cx) * 0.4;
    const cp2y = snk.cy + (cy - snk.cy) * 0.4;

    const pathStr = `M ${src.cx},${src.cy} C ${cp1x},${cp1y} ${cp2x},${cp2y} ${snk.cx},${snk.cy}`;

    const safeType = edge.type || 'direct-call';
    const baseColor = ETYPE_COLOR[safeType] || '#888';
    // Outgoing: full colour, opacity 0.9. Incoming: lighter + 0.6 opacity.
    const strokeColor = isOutgoing
      ? baseColor
      : d3.interpolateRgb(d3.color(baseColor), d3.color('#ddd'))(0.45);
    const strokeOpacity = showAllEdges && !selectedNodeIds.has(edge.from) && !selectedNodeIds.has(edge.to)
      ? 0.15
      : (isOutgoing ? 0.9 : 0.6);

    const dashArr = edgeDashArray(edge);
    const path = overlay.append('path')
      .attr('d', pathStr)
      .attr('fill', 'none')
      .attr('stroke', strokeColor)
      .attr('stroke-width', 1.5)
      .attr('stroke-opacity', strokeOpacity);

    if (dashArr) path.attr('stroke-dasharray', dashArr);

    // Directional arrowheads: outgoing gets marker-end, incoming gets marker-start.
    // Applied only on non-faint (selected) edges for clarity.
    const isFaint = showAllEdges && !selectedNodeIds.has(edge.from) && !selectedNodeIds.has(edge.to);
    if (!isFaint) {
      const typeKey = ETYPE_COLOR[safeType] ? safeType : 'direct-call';
      if (isOutgoing) {
        path.attr('marker-end', `url(#arrow-out-${typeKey})`);
      } else {
        path.attr('marker-start', `url(#arrow-in-${typeKey})`);
      }
    }
  }
}

// ── Logical view ──────────────────────────────────────────────────────────────

function renderLogical() {
  const { w, h } = getViewportSize();
  svg.attr('viewBox', `0 0 ${w} ${h}`);
  svg.selectAll('*').remove();
  logicalNodeCenters.clear();
  physicalNodeCenters.clear();

  const logicals = COMPONENTS_DATA.logical_components || [];

  if (logicals.length === 0) {
    svg.append('text')
      .attr('x', w / 2).attr('y', h / 2)
      .attr('text-anchor', 'middle')
      .attr('fill', '#666').attr('font-size', 14)
      .text('No logical model yet — run /synthesize-audit to build it.');
    return;
  }

  const maxLoc = d3.max(logicals, d => d.size_estimate_loc) || 1;
  const groups = ['core', 'seam', 'removable'].map(cls => ({
    id: cls,
    name: cls,
    children: logicals
      .filter(c => c.classification === cls)
      .map(c => ({ ...c, value: Math.max(c.size_estimate_loc, 30) })),
  })).filter(g => g.children.length > 0);

  const root = buildTreemapHierarchy(groups);
  d3.treemap().size([w, h]).paddingOuter(4).paddingInner(2).paddingTop(18)(root);

  drawGroupLabels(svg, root.children);

  const leaves = root.leaves();

  // Record cell centres for edge overlay
  for (const d of leaves) {
    logicalNodeCenters.set(d.data.id, {
      cx: (d.x0 + d.x1) / 2,
      cy: (d.y0 + d.y1) / 2,
    });
  }

  const g = svg.selectAll('g.node')
    .data(leaves)
    .join('g')
    .attr('class', 'node')
    .attr('transform', d => `translate(${d.x0},${d.y0})`);

  g.append('rect')
    .attr('width',  d => Math.max(0, d.x1 - d.x0))
    .attr('height', d => Math.max(0, d.y1 - d.y0))
    .attr('fill',   d => clsColor(d.data.classification, d.data.size_estimate_loc, maxLoc))
    .attr('rx', 2)
    .on('mousemove', (event, d) => showTip(event, d.data, 'logical'))
    .on('mouseleave', hideTip)
    .on('click', (event, d) => {
      event.stopPropagation();
      if (event.shiftKey) {
        selectedNodeIds.add(d.data.id);
      } else {
        selectedNodeIds = new Set([d.data.id]);
        showPanel(d.data.id, 'logical');
      }
      buildEdgeOverlaySvg(logicals.map(c => c.id), logicalNodeCenters);
    });

  // Cell labels and metric badges
  g.each(function(d) {
    const cw = d.x1 - d.x0;
    const ch = d.y1 - d.y0;
    if (cw < 30 || ch < 14) return;

    const gEl = d3.select(this);
    const metrics = d.data.metrics;

    // Metric badges at top-left
    drawMetricBadges(gEl, metrics, cw, ch);

    // Name label — offset down if badges drawn
    const hasBadges = metrics && badgesForMetrics(metrics).length > 0;
    const textY0 = hasBadges ? 16 : 13;

    const lines = wrapText(d.data.name, Math.floor(cw / 6.5));
    lines.slice(0, 3).forEach((line, i) => {
      gEl.append('text')
        .attr('x', 4).attr('y', textY0 + i * 12)
        .attr('class', i === 0 ? '' : 'sub')
        .text(line);
    });

    if (ch > 32 && cw > 50) {
      gEl.append('text')
        .attr('class', 'sub')
        .attr('x', 4).attr('y', ch - 4)
        .text(d.data.size_estimate_loc + ' loc');
    }
  });

  // Draw any pre-existing selection
  buildEdgeOverlaySvg(logicals.map(c => c.id), logicalNodeCenters);
}

// ── Physical view ─────────────────────────────────────────────────────────────

function renderPhysical() {
  const { w, h } = getViewportSize();
  svg.attr('viewBox', `0 0 ${w} ${h}`);
  svg.selectAll('*').remove();
  physicalNodeCenters.clear();
  logicalNodeCenters.clear();

  const physicals = COMPONENTS_DATA.physical_components || [];
  if (physicals.length === 0) {
    svg.append('text')
      .attr('x', w / 2).attr('y', h / 2)
      .attr('text-anchor', 'middle')
      .attr('fill', '#666').attr('font-size', 14)
      .text('No physical components in this manifest.');
    return;
  }

  const maxLoc = d3.max(physicals, d => d.size_loc) || 1;
  const groups = ['core', 'seam', 'removable'].map(cls => ({
    id: cls,
    name: cls,
    children: physicals
      .filter(pc => primaryClassForPhysical(pc) === cls)
      .map(pc => ({ ...pc, value: Math.max(pc.size_loc, 30), _primaryCls: cls })),
  })).filter(g => g.children.length > 0);

  const root = buildTreemapHierarchy(groups);
  d3.treemap().size([w, h]).paddingOuter(4).paddingInner(2).paddingTop(18)(root);

  drawGroupLabels(svg, root.children);

  const leaves = root.leaves();

  // Record cell centres for edge overlay (keyed on physical component ids)
  for (const d of leaves) {
    physicalNodeCenters.set(d.data.id, {
      cx: (d.x0 + d.x1) / 2,
      cy: (d.y0 + d.y1) / 2,
    });
  }

  const physicalIds = physicals.map(pc => pc.id);

  const g = svg.selectAll('g.node')
    .data(leaves)
    .join('g')
    .attr('class', 'node')
    .attr('transform', d => `translate(${d.x0},${d.y0})`);

  g.append('rect')
    .attr('width',  d => Math.max(0, d.x1 - d.x0))
    .attr('height', d => Math.max(0, d.y1 - d.y0))
    .attr('fill',   d => clsColor(d.data._primaryCls, d.data.size_loc, maxLoc))
    .attr('rx', 2)
    .on('mousemove', (event, d) => showTip(event, d.data, 'physical'))
    .on('mouseleave', hideTip)
    .on('click', (event, d) => {
      event.stopPropagation();
      if (event.shiftKey) {
        selectedNodeIds.add(d.data.id);
      } else {
        selectedNodeIds = new Set([d.data.id]);
        showPanel(d.data.id, 'physical');
      }
      buildEdgeOverlaySvg(physicalIds, physicalNodeCenters);
    });

  // Per-cell labels, owner pills, metric badges
  g.each(function(d) {
    const cw = d.x1 - d.x0;
    const ch = d.y1 - d.y0;
    if (cw < 30 || ch < 14) return;

    const gEl = d3.select(this);
    const metrics = d.data.metrics;

    // Metric badges at top-left
    drawMetricBadges(gEl, metrics, cw, ch);

    const hasBadges = metrics && badgesForMetrics(metrics).length > 0;
    const textY0 = hasBadges ? 16 : 13;

    const base = d.data.path.split('/').pop();
    const lines = wrapText(base, Math.floor(cw / 6.2));
    lines.slice(0, 2).forEach((line, i) => {
      gEl.append('text')
        .attr('x', 4).attr('y', textY0 + i * 11)
        .attr('class', i === 0 ? '' : 'sub')
        .text(line);
    });

    // Owner pills at bottom — fix for pill-accumulation bug:
    // pills are drawn exactly once per g.each iteration (no re-join/re-append
    // accumulation possible because svg.selectAll('*').remove() wipes the SVG
    // before each render call).
    if (ch > 36 && d.data.logical_owners && d.data.logical_owners.length > 0) {
      const owners = d.data.logical_owners;
      const pillH = 10;
      const pillPad = 2;
      let px = 3;
      const py = ch - pillH - 3;

      for (const oid of owners) {
        const lc = logicalById.get(oid);
        if (!lc) continue;
        const label = lc.name.substring(0, 10);
        const pw = label.length * 5.5 + 6;
        if (px + pw > cw - 3) break;  // overflow guard

        const pillColor = d3.color(CLS_BASE[lc.classification] || '#444');
        pillColor.opacity = 0.7;

        gEl.append('rect').attr('class', 'badge-stripe')
          .attr('x', px).attr('y', py)
          .attr('width', pw).attr('height', pillH)
          .attr('rx', 2)
          .attr('fill', pillColor.toString());

        gEl.append('text').attr('class', 'badge-stripe')
          .attr('x', px + 3).attr('y', py + 8)
          .attr('font-size', 7).attr('fill', '#fff')
          .text(label);

        px += pw + pillPad;
      }
    }

    if (ch > 50 && cw > 50) {
      const offsetFromBottom = (d.data.logical_owners && d.data.logical_owners.length > 0) ? 16 : 4;
      gEl.append('text').attr('class', 'sub')
        .attr('x', 4).attr('y', ch - offsetFromBottom)
        .text(d.data.size_loc + ' loc');
    }
  });

  // Draw any pre-existing selection (e.g. after a resize or view switch back)
  buildEdgeOverlaySvg(physicalIds, physicalNodeCenters);
}

// ── Graph view ────────────────────────────────────────────────────────────────

function renderGraph() {
  const { w, h } = getViewportSize();
  svg.attr('viewBox', `0 0 ${w} ${h}`);
  svg.selectAll('*').remove();

  if (graphSimulation) {
    graphSimulation.stop();
    graphSimulation = null;
  }

  const logicals = COMPONENTS_DATA.logical_components || [];
  if (logicals.length === 0) {
    svg.append('text')
      .attr('x', w / 2).attr('y', h / 2)
      .attr('text-anchor', 'middle')
      .attr('fill', '#666').attr('font-size', 14)
      .text('No logical model yet — run /synthesize-audit to build it.');
    return;
  }

  const viewIds = new Set(logicals.map(c => c.id));
  const activeEdges = edgesForView(Array.from(viewIds));

  const maxLoc = d3.max(logicals, d => d.size_estimate_loc) || 1;

  // Node radius = sqrt of normalised LOC, scaled to range [8, 32]
  const rScale = d3.scaleSqrt()
    .domain([0, maxLoc])
    .range([8, 32]);

  // Seed positions by classification zone. The seeded PRNG (from the URL hash)
  // drives BOTH the base jitter here and the perturbation below, so the layout
  // is reproducible per seed — same hash gives the same starting positions.
  const rng = mulberry32(getSeedFromHash());
  const nodes = logicals.map(c => ({
    ...c,
    x: CLUSTER_X[c.classification] * w + (rng() * 60 - 30),
    y: CLUSTER_Y * h + (rng() * 60 - 30),
    r: rScale(c.size_estimate_loc),
  }));

  const nodeById = new Map(nodes.map(n => [n.id, n]));

  const edgeLinks = activeEdges.map(e => ({
    ...e,
    source: e.from,
    target: e.to,
  }));

  // Draw classification zone labels in background
  const zoneLabels = [
    { cls: 'core',      label: 'CORE',      tx: CLUSTER_X.core * w },
    { cls: 'seam',      label: 'SEAM',      tx: CLUSTER_X.seam * w },
    { cls: 'removable', label: 'REMOVABLE', tx: CLUSTER_X.removable * w },
  ];

  const zonelayer = svg.append('g').attr('class', 'zones');
  for (const z of zoneLabels) {
    zonelayer.append('text')
      .attr('x', z.tx).attr('y', 20)
      .attr('text-anchor', 'middle')
      .attr('fill', CLS_BASE[z.cls])
      .attr('font-size', 11)
      .attr('font-weight', 700)
      .attr('opacity', 0.5)
      .text(z.label);
  }

  // Edge layer (below nodes)
  const edgeLayer = svg.append('g').attr('class', 'graph-edges');

  const edgePaths = edgeLayer.selectAll('line.graph-edge')
    .data(edgeLinks)
    .join('line')
    .attr('class', 'graph-edge')
    .attr('stroke', e => ETYPE_COLOR_GRAPH[e.type] || '#aab')
    .attr('stroke-width', 2)
    .attr('stroke-opacity', e => {
      if (!showAllEdges && selectedNodeIds.size > 0) {
        return (selectedNodeIds.has(e.from) || selectedNodeIds.has(e.to)) ? 0.9 : 0.15;
      }
      return showAllEdges ? 0.3 : 0.7;
    })
    .each(function(e) {
      // Widened dash gap so static (solid) vs audit-asserted (dashed) is unmistakable.
      const dashArr = edgeDashArray(e);
      if (dashArr) d3.select(this).attr('stroke-dasharray', '5 4');
    });

  // Node layer
  const nodeLayer = svg.append('g').attr('class', 'graph-nodes');

  const nodeGroups = nodeLayer.selectAll('g.gnode')
    .data(nodes)
    .join('g')
    .attr('class', 'gnode')
    .style('cursor', 'pointer')
    .call(d3.drag()
      .on('start', (event, d) => {
        if (!event.active) graphSimulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
      })
      .on('drag', (event, d) => {
        d.fx = event.x;
        d.fy = event.y;
      })
      .on('end', (event, d) => {
        if (!event.active) graphSimulation.alphaTarget(0);
        // Keep pinned; double-click to unpin
      })
    )
    .on('dblclick', (event, d) => {
      d.fx = null;
      d.fy = null;
      graphSimulation.alpha(0.3).restart();
    })
    .on('click', (event, d) => {
      event.stopPropagation();
      if (event.shiftKey) {
        selectedNodeIds.add(d.id);
      } else {
        selectedNodeIds = new Set([d.id]);
        showPanel(d.id, 'logical');
      }
      refreshGraphEdges(edgePaths);
    })
    .on('mousemove', (event, d) => showTip(event, d, 'logical'))
    .on('mouseleave', hideTip);

  // Halo (classification colour ring)
  nodeGroups.append('circle')
    .attr('r', d => d.r + 3)
    .attr('fill', 'none')
    .attr('stroke', d => CLS_BASE[d.classification] || '#444')
    .attr('stroke-width', 1.5)
    .attr('opacity', 0.6);

  // Node fill
  nodeGroups.append('circle')
    .attr('r', d => d.r)
    .attr('fill', d => clsColor(d.classification, d.size_estimate_loc, maxLoc));

  // Node label (only for nodes large enough)
  nodeGroups.each(function(d) {
    if (d.r < 14) return;
    d3.select(this).append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '0.35em')
      .attr('font-size', Math.min(d.r * 0.55, 11))
      .attr('fill', '#fff')
      .attr('pointer-events', 'none')
      .text(d.name.substring(0, 12));
  });

  // Force simulation — seeded for reproducibility via URL hash. Reuse the same
  // seeded rng created for the base positions above so the whole layout is a
  // pure function of the seed.
  for (const n of nodes) {
    n.x += (rng() - 0.5) * 40;
    n.y += (rng() - 0.5) * 40;
  }

  graphSimulation = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(edgeLinks).id(d => d.id).distance(80).strength(0.3))
    .force('charge', d3.forceManyBody().strength(-200))
    .force('collide', d3.forceCollide(d => d.r + 6))
    // Classification gravity: pull each node toward its zone centre
    .force('cluster-x', d3.forceX(d => CLUSTER_X[d.classification] * w).strength(0.15))
    .force('cluster-y', d3.forceY(h * CLUSTER_Y).strength(0.08))
    .on('tick', () => {
      nodeGroups.attr('transform', d => `translate(${d.x},${d.y})`);
      edgePaths
        .attr('x1', e => nodeById.get(e.from)?.x ?? 0)
        .attr('y1', e => nodeById.get(e.from)?.y ?? 0)
        .attr('x2', e => nodeById.get(e.to)?.x ?? 0)
        .attr('y2', e => nodeById.get(e.to)?.y ?? 0);
    });
}

function refreshGraphEdges(edgePaths) {
  edgePaths.attr('stroke-opacity', e => {
    if (selectedNodeIds.size > 0) {
      const active = selectedNodeIds.has(e.from) || selectedNodeIds.has(e.to);
      return active ? 0.9 : 0.15;
    }
    return showAllEdges ? 0.3 : 0.7;
  });
}

// ── Seeded PRNG (Mulberry32) ──────────────────────────────────────────────────
// Determinism contract: same seed → same layout. Seed lives in URL hash.

function getSeedFromHash() {
  const hash = window.location.hash.replace('#', '');
  const parsed = parseInt(hash, 10);
  return isNaN(parsed) ? 42 : parsed;
}

function mulberry32(seed) {
  let s = seed >>> 0;
  return function() {
    s += 0x6D2B79F5;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ── Tooltip ───────────────────────────────────────────────────────────────────

function showTip(event, data, mode) {
  const name = mode === 'logical' ? data.name : data.path.split('/').pop();
  const cls  = mode === 'logical' ? data.classification : primaryClassForPhysical(data);
  const loc  = mode === 'logical' ? data.size_estimate_loc : data.size_loc;
  const clr  = CLS_BASE[cls] || '#888';
  let extra = '';
  if (mode === 'physical' && data.logical_owners && data.logical_owners.length) {
    const shown = data.logical_owners.slice(0, 3).join(', ');
    const more  = data.logical_owners.length > 3 ? '…' : '';
    extra = `<div style="margin-top:3px;color:#aaa">owners: ${shown}${more}</div>`;
  }
  tip.innerHTML = `<strong>${name}</strong><span class="tt-cls" style="color:${clr}">${cls}</span><div class="tt-size">${loc} loc</div>${extra}`;
  tip.style.display = 'block';
  moveTip(event);
}

function moveTip(event) {
  const x = event.clientX + 12;
  const y = event.clientY + 12;
  tip.style.left = Math.min(x, window.innerWidth - 280) + 'px';
  tip.style.top  = Math.min(y, window.innerHeight - 100) + 'px';
}

function hideTip() {
  tip.style.display = 'none';
}

document.addEventListener('mousemove', e => {
  if (tip.style.display === 'block') moveTip(e);
});

// ── Side panel ────────────────────────────────────────────────────────────────

function showPanel(id, mode) {
  if (mode === 'logical') showLogicalPanel(id);
  else showPhysicalPanel(id);
}

function showLogicalPanel(id) {
  const lc = logicalById.get(id);
  if (!lc) return;

  const edgesOut = (COMPONENTS_DATA.edges || []).filter(e => e.from === id);
  const edgesIn  = (COMPONENTS_DATA.edges || []).filter(e => e.to   === id);
  const prune    = pruneByComp.get(id);

  // Determine "judgement only" — no supporting evidence in manifest
  const hasEvidence = (
    (lc.smoke_findings && lc.smoke_findings.length > 0) ||
    (lc.audit_slices && lc.audit_slices.length > 0) ||
    prune
  );

  panelTitle.textContent = lc.name;

  const classHtml = `<span class="badge ${lc.classification}">${lc.classification}</span>`;
  const judgementPill = hasEvidence ? '' : '<span class="badge" style="background:#222;color:#888;border-color:#444;margin-left:4px">judgement only</span>';

  document.querySelectorAll('#panel-header .badge').forEach(b => b.remove());
  panelTitle.insertAdjacentHTML('afterend', classHtml + judgementPill);

  // View-in-other-mode link
  const otherView = currentView === 'logical' ? 'physical' : 'logical';
  const otherLabel = otherView === 'physical' ? 'Physical' : 'Logical';
  const crossLink = `<div style="font-size:10px;color:#5a9;margin-top:4px;cursor:pointer" onclick="jumpToOtherView('${id}','${otherView}')">→ View in ${otherLabel} mode</div>`;

  let html = '';

  html += `<div class="panel-section"><h3>Description</h3><p>${escHtml(lc.description)}</p>${crossLink}</div>`;

  if (prune) {
    html += `<div class="panel-section"><h3>Prune candidate</h3>`;
    html += `<div class="prune-box">`;
    html += `<span class="loc">${prune.loc_to_remove}</span> <span class="loc-label">loc to remove</span>`;
    html += `<div style="font-size:10px;color:#888;margin-top:2px">${prune.bounded ? 'Bounded removal' : 'Unbounded'}</div>`;
    if (prune.unblocks_simplification_of && prune.unblocks_simplification_of.length) {
      html += `<div class="unblocks"><strong style="color:#888">Unblocks:</strong><ul>`;
      prune.unblocks_simplification_of.forEach(u => { html += `<li>${escHtml(u)}</li>`; });
      html += `</ul></div>`;
    }
    if (prune.latent_bug_fixed) {
      html += `<div class="bug"><strong>Latent bug fixed:</strong> ${escHtml(prune.latent_bug_fixed)}</div>`;
    }
    html += `</div></div>`;
  }

  // Physical files
  html += `<div class="panel-section"><h3>Physical files</h3>`;
  (lc.physical_files || []).forEach(f => {
    html += `<div class="file-row">
      <span class="fp" title="${escHtml(f.path)}">${escHtml(f.path)}</span>
      <span class="lr">${f.line_range ? f.line_range : ''}${f.ownership ? ' (' + f.ownership + ')' : ''}</span>
    </div>`;
  });
  html += `</div>`;

  // Metrics breakdown
  if (lc.metrics && Object.keys(lc.metrics).length) {
    html += metricsBreakdownHtml(lc.metrics);
  }

  // Smoke findings
  if (lc.smoke_findings && lc.smoke_findings.length) {
    html += `<div class="panel-section"><h3>Smoke findings</h3><div>`;
    lc.smoke_findings.forEach(f => {
      const fi = findingIndex[f];
      const title = fi ? fi.title : null;
      const sev   = fi ? fi.severity : null;
      const url   = fi && fi.url ? fi.url : null;
      const sevColor = { info: '#5a9', warn: '#fc0', error: '#f66' }[sev] || '#888';
      if (url) {
        html += `<a href="${escHtml(url)}" target="_blank" style="text-decoration:none">`;
      }
      html += `<span class="finding-pill" title="${title ? escHtml(title) : f}">${f}`;
      if (sev) html += `<span style="color:${sevColor};margin-left:3px;font-size:8px">${sev.toUpperCase()}</span>`;
      html += `</span>`;
      if (url) html += `</a>`;
    });
    html += `</div></div>`;
  }

  // Audit slices
  if (lc.audit_slices && lc.audit_slices.length) {
    html += `<div class="panel-section"><h3>Audit slices</h3><div>`;
    lc.audit_slices.forEach(s => {
      html += `<span class="slice-pill">${escHtml(s)}</span>`;
    });
    html += `</div></div>`;
  }

  // Edges out
  if (edgesOut.length) {
    html += `<div class="panel-section"><h3>Edges out (${edgesOut.length})</h3>`;
    edgesOut.forEach(e => { html += edgeRowHtml(e, 'out', 'logical'); });
    html += `</div>`;
  }

  // Edges in
  if (edgesIn.length) {
    html += `<div class="panel-section"><h3>Edges in (${edgesIn.length})</h3>`;
    edgesIn.forEach(e => { html += edgeRowHtml(e, 'in', 'logical'); });
    html += `</div>`;
  }

  panelBody.innerHTML = html;
}

function showPhysicalPanel(id) {
  const pc = physicalById.get(id);
  if (!pc) return;
  const cls = primaryClassForPhysical(pc);

  panelTitle.textContent = pc.path.split('/').pop();
  const classHtml = `<span class="badge ${cls}">${cls}</span>`;
  document.querySelectorAll('#panel-header .badge').forEach(b => b.remove());
  panelTitle.insertAdjacentHTML('afterend', classHtml);

  let html = '';
  html += `<div class="panel-section"><h3>Path</h3><p style="font-family:monospace;color:#7ab8e0;word-break:break-all">${escHtml(pc.path)}</p></div>`;
  html += `<div class="panel-section"><h3>Size</h3><p>${pc.size_loc} loc &nbsp;·&nbsp; kind: ${pc.kind}</p></div>`;

  if (pc.metrics && Object.keys(pc.metrics).length) {
    html += metricsBreakdownHtml(pc.metrics);
  }

  if (pc.logical_owners && pc.logical_owners.length) {
    html += `<div class="panel-section"><h3>Logical owners</h3>`;
    html += `<div style="font-size:10px;color:#5a9;margin-bottom:4px;cursor:pointer" onclick="this.nextElementSibling.style.display=this.nextElementSibling.style.display==='none'?'flex':'none'">(click owner to jump to Logical view)</div>`;
    html += `<div class="owners-list">`;
    (pc.logical_owners || []).forEach(oid => {
      const lc = logicalById.get(oid);
      if (!lc) {
        html += `<span class="owner-pill core" style="background:#1a1a3a;color:#888">${escHtml(oid)} (not in model)</span>`;
        return;
      }
      // "→ View in Logical mode" is the cross-link affordance from Q-1.3
      html += `<span class="owner-pill ${lc.classification}" style="cursor:pointer"
                 title="→ View in Logical mode"
                 onclick="jumpToOtherView('${oid}','logical')">${escHtml(lc.name)}</span>`;
    });
    html += `</div></div>`;
  }

  // Edges out / in — resolved against physical component ids
  const edgesOut = (COMPONENTS_DATA.edges || []).filter(e => e.from === id);
  const edgesIn  = (COMPONENTS_DATA.edges || []).filter(e => e.to   === id);

  if (edgesOut.length) {
    html += `<div class="panel-section"><h3>Edges out (${edgesOut.length})</h3>`;
    edgesOut.forEach(e => { html += edgeRowHtml(e, 'out', 'physical'); });
    html += `</div>`;
  }

  if (edgesIn.length) {
    html += `<div class="panel-section"><h3>Edges in (${edgesIn.length})</h3>`;
    edgesIn.forEach(e => { html += edgeRowHtml(e, 'in', 'physical'); });
    html += `</div>`;
  }

  panelBody.innerHTML = html;
}

function metricsBreakdownHtml(metrics) {
  const rows = [];
  const def = [
    ['loc',        'LOC'],
    ['fan_in',     'Fan-in'],
    ['fan_out',    'Fan-out'],
    ['cyclomatic', 'Cyclomatic'],
    ['churn_90d',  'Churn 90d'],
    ['test_ratio', 'Test ratio'],
  ];
  for (const [field, label] of def) {
    if (metrics[field] !== undefined) {
      rows.push(`<div class="file-row"><span style="color:#aaa">${label}</span><span style="color:#ddd">${metrics[field]}</span></div>`);
    }
  }
  if (!rows.length) return '';
  return `<div class="panel-section"><h3>Metrics</h3>${rows.join('')}</div>`;
}

// Resolve a display name for an edge peer id, checking physical then logical maps.
// In physical mode the ids are physical file paths; in logical mode they are logical names.
function peerDisplayName(peerId, mode) {
  if (mode === 'physical') {
    const pc = physicalById.get(peerId);
    if (pc) return pc.path.split('/').pop();
  }
  const lc = logicalById.get(peerId);
  if (lc) return lc.name;
  return peerId;
}

// Render a single edge row for the side panel.
//   e    — edge object
//   dir  — 'out' or 'in'
//   mode — 'logical' or 'physical' (controls id-space for peer name resolution)
function edgeRowHtml(e, dir, mode) {
  const peer   = dir === 'out' ? e.to : e.from;
  const peerNm = peerDisplayName(peer, mode || 'logical');
  const arrow  = dir === 'out' ? '→' : '←';
  const panelMode = mode || 'logical';

  // Evidence class annotation
  const ecLabel = e.evidence_class
    ? `<span style="font-size:9px;color:#555;margin-left:4px">[${e.evidence_class}]</span>`
    : '';

  let html = `<div class="edge-row ${e.type}">`;
  html += `<span class="peer" style="cursor:pointer" onclick="showPanel('${peer}','${panelMode}')">${arrow} ${escHtml(peerNm)}</span>`;
  html += `<span class="etype" title="${escHtml(EDGE_TYPE_HELP[e.type] || '')}" style="cursor:help">${e.type}</span>`;
  html += ecLabel;
  if (e.smell) {
    html += `<div class="smell">⚠ ${escHtml(e.smell)}</div>`;
  }
  if (e.smoke_finding) {
    const fi = findingIndex[e.smoke_finding];
    const fiTitle = fi ? fi.title : e.smoke_finding;
    html += `<div class="smell" style="color:#5cd65c">↳ ${escHtml(e.smoke_finding)}${fi ? ' — ' + escHtml(fiTitle) : ''}</div>`;
  }
  if (e.evidence) {
    html += `<div class="evidence">${escHtml(e.evidence)}</div>`;
  }
  html += `</div>`;
  return html;
}

// Cross-view jump (Q-1.3 cross-link affordance)
// setView() resets selectedNodeIds; re-select AFTER setView so the target's
// edges light up immediately on arrival.
function jumpToOtherView(id, targetView) {
  setView(targetView);
  // setView cleared selectedNodeIds — re-select the target node now.
  selectedNodeIds = new Set([id]);
  // Rebuild the edge overlay for the destination view.
  if (targetView === 'logical') {
    const viewIds = (COMPONENTS_DATA.logical_components || []).map(c => c.id);
    buildEdgeOverlaySvg(viewIds, logicalNodeCenters);
  } else if (targetView === 'physical') {
    const viewIds = (COMPONENTS_DATA.physical_components || []).map(c => c.id);
    buildEdgeOverlaySvg(viewIds, physicalNodeCenters);
  }
  showPanel(id, targetView === 'logical' ? 'logical' : 'physical');
}

// ── Edge-type pedagogy (Tier-1 popover content) ───────────────────────────────

// This content is sourced from references/edge-types.md and condensed for popovers.
const EDGE_TYPE_HELP = {
  'direct-call': [
    'One component calls another by name. The caller\'s source code contains the callee\'s identifier.',
    'Connascence of Name — the weakest form of coupling.',
    'Normal and necessary. Only a problem when the callee is itself removable.',
    'When removable: only when the callee is classified removable.',
  ].join(' | '),
  'shared-state': [
    'Two components both read and write the same data structure (global, singleton, DOM node).',
    'Connascence of Identity or Value.',
    'Tolerable at small scale; corrosive as the structure grows. Out-of-order writes cause silent bugs.',
    'When removable: when one of the two writers can be eliminated.',
  ].join(' | '),
  'background-knowledge': [
    'Component C depends on an invariant in component D that is not enforced in code, types, or tests.',
    'Connascence of Convention — the most expensive form.',
    'Always a smell: the compiler, type checker, and test suite cannot catch a violation.',
    'When removable: always, by pruning the component that carries the implicit invariant.',
  ].join(' | '),
};

function buildLegendPopovers() {
  // Attach popover titles to legend swatches by data-type attribute.
  // The legend HTML uses data-etype on popover-trigger elements.
  document.querySelectorAll('[data-etype]').forEach(el => {
    const etype = el.getAttribute('data-etype');
    if (EDGE_TYPE_HELP[etype]) {
      el.setAttribute('title', EDGE_TYPE_HELP[etype]);
    }
  });
}

// ── "How to read this" panel (Tier-2, localStorage state) ────────────────────

// ── Legend toggle (collapsed "?" by default) ─────────────────────────────────

const LEGEND_STORAGE_KEY = 'architecture-treemap-legend-open';

function toggleLegend() {
  const panel    = document.getElementById('legend');
  const launcher = document.getElementById('legend-launcher');
  if (!panel) return;
  const nowOpen = panel.style.display !== 'block';
  panel.style.display = nowOpen ? 'block' : 'none';
  if (launcher) launcher.style.display = nowOpen ? 'none' : 'flex';
  localStorage.setItem(LEGEND_STORAGE_KEY, String(nowOpen));
}

function initLegend() {
  const panel    = document.getElementById('legend');
  const launcher = document.getElementById('legend-launcher');
  if (!panel) return;
  const stored = localStorage.getItem(LEGEND_STORAGE_KEY);
  // Collapsed by default on first visit (stored === null → false)
  const isOpen = stored === 'true';
  panel.style.display   = isOpen ? 'block' : 'none';
  if (launcher) launcher.style.display = isOpen ? 'none' : 'flex';
}

function initHowToReadPanel() {
  const panel  = document.getElementById('how-to-read-panel');
  const toggle = document.getElementById('how-to-read-toggle');
  if (!panel || !toggle) return;

  const storageKey = 'architecture-treemap-how-to-read-open';
  const stored = localStorage.getItem(storageKey);
  // Closed by default on first visit
  const isOpen = stored === 'true';
  panel.style.display = isOpen ? 'block' : 'none';
  toggle.textContent  = isOpen ? '▾' : '▸';

  toggle.addEventListener('click', () => {
    const nowOpen = panel.style.display !== 'block';
    panel.style.display = nowOpen ? 'block' : 'none';
    toggle.textContent  = nowOpen ? '▾' : '▸';
    localStorage.setItem(storageKey, String(nowOpen));
  });

  // Populate example from data (background-knowledge edges only)
  const bgEdges = (COMPONENTS_DATA.edges || []).filter(e => e.type === 'background-knowledge');
  const exampleSlot = document.getElementById('how-to-read-example');
  if (exampleSlot) {
    if (bgEdges.length > 0) {
      const ex = bgEdges[0];
      const fromLc = logicalById.get(ex.from);
      const toLc   = logicalById.get(ex.to);
      const fromNm = fromLc ? fromLc.name : ex.from;
      const toNm   = toLc   ? toLc.name   : ex.to;
      const evText = ex.evidence ? ex.evidence.substring(0, 120) + (ex.evidence.length > 120 ? '…' : '') : '';
      exampleSlot.innerHTML = `Example: <strong>${escHtml(fromNm)} → ${escHtml(toNm)}</strong>${evText ? ' — ' + escHtml(evText) : ''}`;
    } else {
      exampleSlot.style.display = 'none';
    }
  }
}

// ── View toggle ───────────────────────────────────────────────────────────────

function setView(view) {
  currentView = view;
  selectedNodeIds = new Set();

  document.getElementById('btn-logical').classList.toggle('active', view === 'logical');
  document.getElementById('btn-physical').classList.toggle('active', view === 'physical');
  document.getElementById('btn-graph').classList.toggle('active', view === 'graph');

  // Reset panel
  panelTitle.textContent = 'Click a cell to inspect';
  document.querySelectorAll('#panel-header .badge').forEach(b => b.remove());
  panelBody.innerHTML = '<p class="empty-state">Select a component from the treemap to see its details, files, and edges.</p>';

  if (view === 'logical')  renderLogical();
  else if (view === 'physical') renderPhysical();
  else if (view === 'graph')    renderGraph();
}

// ── "Show all edges" toggle ───────────────────────────────────────────────────

function toggleAllEdges(checked) {
  showAllEdges = checked;
  if (currentView === 'logical') {
    const viewIds = (COMPONENTS_DATA.logical_components || []).map(c => c.id);
    buildEdgeOverlaySvg(viewIds, logicalNodeCenters);
  } else if (currentView === 'physical') {
    const viewIds = (COMPONENTS_DATA.physical_components || []).map(c => c.id);
    buildEdgeOverlaySvg(viewIds, physicalNodeCenters);
  }
}

// ── Background click clears edge selection ────────────────────────────────────

document.getElementById('treemap-area').addEventListener('click', () => {
  selectedNodeIds = new Set();
  if (currentView === 'logical') {
    buildEdgeOverlaySvg((COMPONENTS_DATA.logical_components || []).map(c => c.id), logicalNodeCenters);
  } else if (currentView === 'physical') {
    buildEdgeOverlaySvg((COMPONENTS_DATA.physical_components || []).map(c => c.id), physicalNodeCenters);
  }
});

// ── Escape-HTML helper ────────────────────────────────────────────────────────

function escHtml(str) {
  if (!str) return '';
  return String(str)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}

// ── Resize handler ────────────────────────────────────────────────────────────

window.addEventListener('resize', () => {
  if (currentView === 'logical')       renderLogical();
  else if (currentView === 'physical') renderPhysical();
  else if (currentView === 'graph')    renderGraph();
});

// ── Initialise ────────────────────────────────────────────────────────────────

buildLegendPopovers();
initLegend();
initHowToReadPanel();
renderLogical();
