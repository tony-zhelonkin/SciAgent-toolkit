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
  core:        '#2a7a2a',
  seam:        '#c87800',
  removable:   '#a02020',
  unclassified:'#555555',  // neutral grey: an unowned physical file inherits NO
                           // classification. It must NOT masquerade as core.
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

// ── Metric model (single source of truth, descriptor-driven) ──────────────────
// The metric registry (scripts/metric_registry.py) is serialised into the
// manifest's `metric_descriptors` block. The badge set, side-panel rows, and
// good/bad semantics are all DERIVED from it here, so the render is a pure
// function of the file and adding a metric never touches this renderer.
//
// Backward compatibility: pre-registry manifests omit the block; we fall back
// to a built-in default table matching the original six named metrics.
const DEFAULT_METRIC_DESCRIPTORS = {
  loc:         { label: 'LOC',        short_badge_label: null,         unit: 'loc',     higher_is_better: false, render: 'hidden', badge_threshold: null, badge_color: '#666' },
  fan_in:      { label: 'Fan-in',     short_badge_label: null,         unit: '',        higher_is_better: false, render: 'hidden', badge_threshold: null, badge_color: '#666' },
  fan_out:     { label: 'Fan-out',    short_badge_label: 'hi fan-out', unit: '',        higher_is_better: false, render: 'badge',  badge_threshold: 5,    badge_color: '#6060a0' },
  cyclomatic:  { label: 'Cyclomatic', short_badge_label: 'hi CC',      unit: '',        higher_is_better: false, render: 'badge',  badge_threshold: 10,   badge_color: '#a06020' },
  churn_90d:   { label: 'Churn 90d',  short_badge_label: 'hi churn',   unit: 'commits', higher_is_better: false, render: 'badge',  badge_threshold: 10,   badge_color: '#405080' },
  test_ratio:  { label: 'Test ratio', short_badge_label: 'low tests',  unit: '',        higher_is_better: true,  render: 'badge',  badge_threshold: 0.3,  badge_color: '#a04040' },
};

const METRIC_DESCRIPTORS = (COMPONENTS_DATA.metric_descriptors && Object.keys(COMPONENTS_DATA.metric_descriptors).length)
  ? COMPONENTS_DATA.metric_descriptors
  : DEFAULT_METRIC_DESCRIPTORS;

// A badge lights when a value is on the WORSE side of its threshold:
// higher_is_better metrics warn BELOW threshold; others warn AT/ABOVE.
function metricBadgeLights(desc, value) {
  if (!desc || desc.render !== 'badge' || desc.badge_threshold === null || desc.badge_threshold === undefined) return false;
  if (value === undefined || value === null) return false;
  return desc.higher_is_better ? (value < desc.badge_threshold) : (value >= desc.badge_threshold);
}

// ── Pedagogy model (single source of truth, manifest-driven) ──────────────────
// scripts/pedagogy_registry.py serialises a `pedagogy` block into the manifest:
// classification / edge-type / evidence-class explainers, the direction
// encoding, per-metric explainers, the epistemic-source legend, and the
// glossary. The renderer reads it so every (?) affordance and the glossary are
// a PURE FUNCTION of the file — the same invariant the metric descriptors hold.
// This REPLACES the old hardcoded EDGE_TYPE_HELP (JS) and the classification /
// evidence help that was inlined in the template HTML (the Connascence-of-Value
// smear the metrics ADR killed for metrics, now killed for pedagogy too).
//
// Backward compatibility: pre-pedagogy manifests omit the block; we fall back
// to a compact built-in default so older manifests still render their (?) help.
const DEFAULT_PEDAGOGY = {
  epistemic_sources: {
    'measured':             { label: 'measured', what: 'Deterministic: an AST/git/tool proved this directly.', color: '#3a7a3a' },
    'metric-anchored':      { label: 'metric-anchored judgment', what: 'Judgment biased by a metric — the metric points, it does not decide.', color: '#a07020' },
    'requires-your-intent': { label: 'requires your intent', what: 'Extrinsic: not in the code. Only a human naming a concern can ground it.', color: '#8a4a8a' },
  },
  classifications: {
    core:      { what: 'Load-bearing: the app cannot function without it.', how: 'Proposed from import graph + metrics (high fan-in, central).', teaches: 'Last to refactor. "Essential" is partly a product claim, not pure code.', epistemic_source: 'metric-anchored', color: '#2a7a2a', provisional_label: 'core (structural hypothesis)' },
    seam:      { what: 'A deliberately thin boundary between two larger parts.', how: 'Weakly detectable: small LOC + fan-in hint at it; deliberateness is not measurable.', teaches: 'Always a candidate — confirm intent. Never an asserted detection.', epistemic_source: 'requires-your-intent', color: '#c87800', provisional_label: 'seam? (confirm intent)' },
    removable: { what: 'A strip candidate: dead, superseded, or dead-on-arrival.', how: 'NOT in the code. Requires a user-stated concern (/audit-slice).', teaches: 'Meaningless without intent. A structural pass can never assert it.', epistemic_source: 'requires-your-intent', color: '#a02020', provisional_label: 'requires your judgment — run /audit-slice' },
  },
  edge_types: {
    'direct-call':          { what: "One component names another's identifier.", how: 'Detected statically by the AST import graph.', teaches: 'Connascence of Name — weakest, normal. Only a problem when the callee is removable.', epistemic_source: 'measured', color: '#888' },
    'shared-state':         { what: 'Two components both read/write the same structure.', how: 'Audit-asserted: not statically detectable.', teaches: 'Connascence of Identity/Value. Corrosive as the structure grows.', epistemic_source: 'metric-anchored', color: '#3a86c0' },
    'background-knowledge': { what: 'An unenforced invariant one assumes about another.', how: 'Requires intent: only a human who knows the contract can name it.', teaches: 'Connascence of Convention — most expensive; no tool can catch a violation.', epistemic_source: 'requires-your-intent', color: '#8a5cb0' },
  },
  evidence_classes: {
    'static':         { what: 'The AST/import graph proved this edge.', how: 'Emitted by the extractor. Renders SOLID.', teaches: 'Solid = provable. Cannot prove runtime/dynamic coupling.', epistemic_source: 'measured' },
    'audit-asserted': { what: 'A human asserted this edge from a slice.', how: 'Authored during synthesis. Renders DASHED.', teaches: 'Dashed = judgment, as trustworthy as the slice behind it.', epistemic_source: 'metric-anchored' },
  },
  direction_encoding: { what: 'Outgoing edges full-colour; incoming lighter.', how: 'A renderer convention, not data.', teaches: 'Direction shows who depends on whom — drives stability.', epistemic_source: 'measured' },
  metrics: {},
  glossary: [],
};

const PEDAGOGY = (COMPONENTS_DATA.pedagogy && Object.keys(COMPONENTS_DATA.pedagogy).length)
  ? COMPONENTS_DATA.pedagogy
  : DEFAULT_PEDAGOGY;

// Merge fallbacks so a partial manifest pedagogy still has the sub-blocks.
for (const key of Object.keys(DEFAULT_PEDAGOGY)) {
  if (PEDAGOGY[key] === undefined) PEDAGOGY[key] = DEFAULT_PEDAGOGY[key];
}

// Build a (?) explainer title string from a {what, how, teaches,
// epistemic_source} entry. The epistemic-source chip prefix is the load-bearing
// part; the prose is the trimmable a/b/c contract.
function explainerTitle(entry) {
  if (!entry) return '';
  const parts = [];
  const src = entry.epistemic_source;
  if (src && PEDAGOGY.epistemic_sources[src]) {
    parts.push('[' + PEDAGOGY.epistemic_sources[src].label + ']');
  }
  if (entry.what) parts.push(entry.what);
  if (entry.how) parts.push('How: ' + entry.how);
  if (entry.teaches) parts.push('Teaches: ' + entry.teaches);
  return parts.join(' — ');
}

// Small inline (?) affordance as an HTML string. Uniform across metric badges,
// classifications, edge-types, evidence-classes, and the direction legend.
// The entry is JSON-serialised into a data attribute so the delegated click
// handler can reconstruct the structured content for the popover (no id lookup).
function infoIcon(entry) {
  if (!entry) return '';
  // Fallback title for accessibility / environments where JS events are absent.
  const fallbackTitle = escHtml(explainerTitle(entry));
  const encoded = escHtml(JSON.stringify(entry));
  return `<span class="ped-info" data-ped-entry="${encoded}" title="${fallbackTitle}" tabindex="0" role="button" aria-label="More info">(?)</span>`;
}

// A small epistemic-source chip (coloured) for a classification/edge panel row.
function epistemicChip(source) {
  const meta = PEDAGOGY.epistemic_sources[source];
  if (!meta) return '';
  return `<span class="epi-chip" style="border-color:${meta.color};color:${meta.color}"
            title="${escHtml(meta.what)}">${escHtml(meta.label)}</span>`;
}

// ── Explainer popover (click-toggle floating panel for all (?) affordances) ───
//
// A single reusable <div id="ped-popover"> is shared by all (?) icons.
// Each .ped-info span carries its pedagogy entry JSON-encoded in data-ped-entry;
// openExplainerPopover() decodes it and builds the labelled structured panel.
// Only one popover is open at a time; click elsewhere or press Esc to dismiss.

function openExplainerPopover(anchorEl, entry) {
  const pop = document.getElementById('ped-popover');
  const body = document.getElementById('ped-popover-body');
  if (!pop || !body) return;

  // Build the structured content — labelled sections, not a mashed string.
  let html = '';

  const src = entry.epistemic_source;
  const srcMeta = src && PEDAGOGY.epistemic_sources && PEDAGOGY.epistemic_sources[src];
  if (srcMeta) {
    html += `<div class="pop-section">
      <span class="pop-epi" style="border-color:${srcMeta.color};color:${srcMeta.color}">${escHtml(srcMeta.label)}</span>
    </div>`;
  }
  if (entry.what) {
    html += `<div class="pop-section"><div class="pop-label">What</div><div class="pop-body">${escHtml(entry.what)}</div></div>`;
  }
  if (entry.how) {
    html += `<div class="pop-section"><div class="pop-label">How it is detected</div><div class="pop-body">${escHtml(entry.how)}</div></div>`;
  }
  if (entry.teaches) {
    html += `<div class="pop-section"><div class="pop-label">Architectural lesson</div><div class="pop-body">${escHtml(entry.teaches)}</div></div>`;
  }

  body.innerHTML = html || '<span style="color:#555">No detail available.</span>';

  // Position the popover near the anchor element, keeping it within the viewport.
  const rect = anchorEl.getBoundingClientRect();
  const popW = 330;  // slightly wider than max-width so the clamp below is the
                     // effective cap; avoids a race before the element has a real
                     // offsetWidth from the DOM.
  const spaceRight = window.innerWidth - rect.right;
  const left = spaceRight >= popW
    ? rect.right + 6
    : Math.max(4, rect.left - popW - 6);
  const top = Math.min(rect.top, window.innerHeight - 200);

  pop.style.left    = left + 'px';
  pop.style.top     = top  + 'px';
  pop.style.display = 'block';

  // Store the triggering element so a second click on the same icon closes it.
  pop._anchorEl = anchorEl;
}

function closeExplainerPopover() {
  const pop = document.getElementById('ped-popover');
  if (pop) {
    pop.style.display = 'none';
    pop._anchorEl = null;
  }
}

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
let graphAmbientIds   = new Set();   // cross-cutting sinks in the current graph
// Session-persistent node positions: populated on drag-end so that re-entering
// the Graph view restores the user's manual layout instead of re-seeding.
// Cleared only on page reload; first render (empty map) is byte-identical to
// the deterministic seeded layout.
const graphNodePositions = new Map(); // node id → {x, y}

// "Show tests" toggle (Physical view). Tests stay in the manifest (kind:test)
// but are excluded from the default backbone layout; this reveals them in a
// de-emphasised style. Persisted in localStorage like the other UI toggles.
const SHOW_TESTS_STORAGE_KEY = 'architecture-treemap-show-tests';
let showTests = (localStorage.getItem(SHOW_TESTS_STORAGE_KEY) === 'true');

// Backbone kinds shown in the default Physical view (everything but tests).
const BACKBONE_KINDS = new Set(['module', 'config', 'asset', 'other', 'doc']);

// Physical-view colour mode: 'role' (classification hue, default) or 'pressure'
// (a single-hue heat ramp keyed to refactor pressure). The toggle TIME-SHARES
// the colour channel — it never permanently overrides the classification hue —
// so the user can flip between "what role does this play?" and "where is the
// pressure?" without the two signals fighting. Persisted like other toggles.
const COLOR_MODE_STORAGE_KEY = 'architecture-treemap-color-mode';
let colorMode = (localStorage.getItem(COLOR_MODE_STORAGE_KEY) === 'pressure') ? 'pressure' : 'role';

// The pressure metric key: the genuine refactor_pressure when radon was present,
// else the degraded LOC proxy (clearly labelled in the UI). Resolved per render.
function pressureKeyFor(metrics) {
  if (!metrics) return null;
  if (metrics.refactor_pressure !== undefined) return 'refactor_pressure';
  if (metrics.refactor_pressure_loc_proxy !== undefined) return 'refactor_pressure_loc_proxy';
  return null;
}

// True when NO component carries a real refactor_pressure AND radon was flagged
// unavailable — the lens must render GREYED with an install hint, not blank.
function pressureUnavailable() {
  const unavail = (COMPONENTS_DATA.extractor && COMPONENTS_DATA.extractor.unavailable_metrics) || [];
  return unavail.indexOf('refactor_pressure') !== -1;
}

// Rank physical components by pressure (real or proxy), descending. Returns
// [{pc, value, key}] for components that HAVE a pressure value.
function topByPressure(physicals, limit) {
  const scored = [];
  for (const pc of physicals) {
    const key = pressureKeyFor(pc.metrics);
    if (!key) continue;
    scored.push({ pc, value: pc.metrics[key], key });
  }
  scored.sort((a, b) => b.value - a.value);
  return limit ? scored.slice(0, limit) : scored;
}

// Heat colour for a pressure value, interpolated cool→hot over [0, max].
function pressureColor(value, maxPressure) {
  const t = maxPressure > 0 ? Math.min(1, value / maxPressure) : 0;
  return d3.interpolateRgb(d3.color('#26323f'), d3.color('#ff5a36'))(0.15 + 0.85 * t);
}

// Cross-cutting / ambient sink detection. A node that many things depend on but
// that depends on nothing is a pure SINK — drawing every edge into it solid
// tangles the view into a hairball. When that sink is also cross-cutting
// (a shared config/kernel/util the whole codebase leans on), its incoming edges
// are drawn faintly so it reads as "everyone depends on this" without drowning
// the core/seam/removable structure. Edges still light to full opacity on
// selection, so the relationship is visible-on-demand, never hidden.
//
// Driven entirely by NODE PROPERTIES, never a hardcoded id:
//   high in-degree  (>= AMBIENT_MIN_FAN_IN, computed from the active view's edges)
//   + zero out-degree (a pure sink)
//   + a cross-cutting marker (config-kind file, all-config ownership, or a
//     config/util id hint — the judgment layer's signal for "shared kernel").
// In/out degree work in both id-spaces (logical-id and physical-id edges).
const AMBIENT_MIN_FAN_IN = 5;

function ambientNodeIds(viewIds, edges) {
  const inDeg = new Map();
  const outDeg = new Map();
  for (const id of viewIds) { inDeg.set(id, 0); outDeg.set(id, 0); }
  for (const e of edges) {
    if (outDeg.has(e.from)) outDeg.set(e.from, outDeg.get(e.from) + 1);
    if (inDeg.has(e.to))    inDeg.set(e.to,    inDeg.get(e.to)    + 1);
  }
  const ambient = new Set();
  for (const id of viewIds) {
    if (inDeg.get(id) >= AMBIENT_MIN_FAN_IN && outDeg.get(id) === 0 && isCrossCuttingNode(id)) {
      ambient.add(id);
    }
  }
  return ambient;
}

// True when a node carries a cross-cutting marker. Layered signals, strongest
// first: a config-kind physical file; a logical node owning only config files;
// otherwise a config/util id hint (the judgment layer names shared kernels this
// way). NOTE config-as-code lives in .py files (kind 'module'), so kind alone is
// insufficient — the id/ownership hint is what catches a config.py sink.
function isCrossCuttingNode(id) {
  const pc = physicalById.get(id);
  if (pc) {
    if (pc.kind === 'config') return true;
    return /(^|[-_/])config(\.py)?$/.test(pc.path) || /(^|[-_/])config$/.test(id);
  }
  const lc = logicalById.get(id);
  if (lc) {
    const files = lc.physical_files || [];
    const allConfig = files.length > 0 && files.every(f => {
      const owner = (COMPONENTS_DATA.physical_components || []).find(p => p.path === f.path);
      return owner && owner.kind === 'config';
    });
    if (allConfig) return true;
    return /(^|-)config$/.test(id) || /(^|-)util(s)?$/.test(id);
  }
  return false;
}

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
  // A physical file with no logical owner is UNCLASSIFIED — it has not been
  // grounded by any synthesis judgment. It renders neutral grey, never 'core'
  // (the old bug: unowned files silently inherited core and inflated the
  // load-bearing class). The non-source filter removes most junk; legit unowned
  // files (e.g. a freshly-extracted substrate before synthesis) stay honest.
  if (owners.length === 0) return 'unclassified';
  // Priority: removable > seam > core (most architecturally interesting wins)
  const rank = { removable: 3, seam: 2, core: 1 };
  let best = 'unclassified', bestRank = 0;
  for (const id of owners) {
    const lc = logicalById.get(id);
    if (!lc) continue;
    const r = rank[lc.classification] || 0;
    if (r > bestRank) { bestRank = r; best = lc.classification; }
  }
  return best;
}

// Approximate glyph width for the tile font (px per char). Used to convert a
// pixel budget to a character budget consistently across wrap + truncate.
const PX_PER_CHAR = 6.2;

// Truncate a single token to a character budget, appending an ellipsis when it
// overflows. The load-bearing fix for the bug where a long single token fell
// through wrapText untruncated and painted past the tile.
function truncateToChars(str, maxChars) {
  if (!str) return '';
  if (maxChars < 1) return '';
  if (str.length <= maxChars) return str;
  if (maxChars <= 1) return '…';
  return str.substring(0, maxChars - 1) + '…';
}

function wrapText(str, maxChars) {
  if (!str) return [''];
  if (maxChars < 1) return [];
  const words = str.replace(/_/g, ' ').split(' ');
  const lines = [];
  let cur = '';
  for (const w of words) {
    if ((cur + ' ' + w).trim().length > maxChars) {
      if (cur) lines.push(cur);
      // A single word longer than the budget is TRUNCATED to the budget with an
      // ellipsis (was: pushed whole, then painted past the tile — the bug).
      cur = (w.length > maxChars) ? truncateToChars(w, maxChars) : w;
    } else {
      cur = (cur + ' ' + w).trim();
    }
  }
  if (cur) lines.push(cur);
  // Fallthrough also truncates rather than emitting an over-budget substring.
  return lines.length ? lines : [truncateToChars(str, maxChars)];
}

// Define (idempotently) a per-cell clipPath sized to the tile, and return its
// url(#id) reference. The LOAD-BEARING overlap fix: clipping a cell's label +
// badge children to the tile box GUARANTEES no element ever paints onto a
// neighbour, independent of every wrap/truncate heuristic. Ids are namespaced
// per view so logical and physical clips never collide.
function ensureCellClip(defs, idPrefix, id, w, h) {
  const clipId = `clip-${idPrefix}-${id}`;
  const cp = defs.append('clipPath').attr('id', clipId);
  cp.append('rect')
    .attr('x', 0).attr('y', 0)
    .attr('width', Math.max(0, w)).attr('height', Math.max(0, h));
  return `url(#${clipId})`;
}

// Derive the lit badges for a metrics block from the descriptor registry.
// Returns up to 3 {label, color} chips, iterating descriptors in declared order.
function badgesForMetrics(metrics) {
  if (!metrics) return [];
  const badges = [];
  for (const field of Object.keys(METRIC_DESCRIPTORS)) {
    const desc = METRIC_DESCRIPTORS[field];
    if (!desc.short_badge_label) continue;
    if (metricBadgeLights(desc, metrics[field])) {
      badges.push({ label: desc.short_badge_label, color: desc.badge_color });
      if (badges.length >= 3) break;
    }
  }
  return badges;
}

// Draw up to 3 metric badges in the top-left of a cell (SVG g element).
// Returns the Y coordinate BELOW the badge row (the reserved bottom), so the
// caller can place the label beneath the badges instead of colliding with them.
// A badge that would overflow the tile width is dropped (not clipped mid-word),
// and the per-cell clipPath (ensureCellClip) is the backstop either way.
function drawMetricBadges(gEl, metrics, cellW, cellH) {
  const BADGE_TOP = 3;
  const badgeH = 9;
  if (!metrics || cellW < 40 || cellH < 16) return BADGE_TOP;
  const badges = badgesForMetrics(metrics);
  if (badges.length === 0) return BADGE_TOP;

  const pad = 2;
  let bx = 3;
  const by = BADGE_TOP;
  let drewAny = false;

  for (const def of badges) {
    const label = def.label;
    const bw = label.length * 5 + 6;
    if (bx + bw > cellW - 3) break;  // no room for this badge; stop the row.
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
    drewAny = true;
  }
  return drewAny ? (by + badgeH + 2) : BADGE_TOP;
}

// Define (once per render) the diagonal-hatch pattern used to mark test tiles
// as visually subordinate. Idempotent: removes any prior copy first.
function ensureTestHatchPattern() {
  svg.select('defs.test-hatch-defs').remove();
  const defs = svg.append('defs').attr('class', 'test-hatch-defs');
  const pattern = defs.append('pattern')
    .attr('id', 'test-hatch')
    .attr('patternUnits', 'userSpaceOnUse')
    .attr('width', 6).attr('height', 6)
    .attr('patternTransform', 'rotate(45)');
  pattern.append('rect')
    .attr('width', 6).attr('height', 6)
    .attr('fill', 'none');
  pattern.append('line')
    .attr('x1', 0).attr('y1', 0).attr('x2', 0).attr('y2', 6)
    .attr('stroke', '#000').attr('stroke-width', 2).attr('stroke-opacity', 0.35);
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

  // Shared defs for per-cell clip paths (overlap fix).
  const cellDefs = svg.append('defs').attr('class', 'cell-clip-defs');

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
    // Raised gate: need width for at least one truncated label + height for a
    // text row. A 31px tile cannot hold a name; skip it rather than overflow.
    if (cw < 46 || ch < 14) return;

    const gEl = d3.select(this);
    const metrics = d.data.metrics;

    // Clip everything in this cell to its own box (load-bearing overlap fix).
    const clipUrl = ensureCellClip(cellDefs, 'log', d.data.id, cw, ch);
    const inner = gEl.append('g').attr('clip-path', clipUrl);

    // Metric badges at top-left; returns the reserved bottom of the badge row.
    const badgeBottom = drawMetricBadges(inner, metrics, cw, ch);
    const hasBadges = badgeBottom > 3;

    // If badges were drawn but the tile is too short to clear them, show badges
    // only (no name) rather than colliding the name into the badge row.
    if (hasBadges && ch < 28) return;

    const textY0 = hasBadges ? badgeBottom + 9 : 13;
    const maxChars = Math.floor((cw - 8) / PX_PER_CHAR);
    const lines = wrapText(d.data.name, maxChars);
    // Bound line count to the vertical space remaining below the badge row.
    const maxLines = Math.max(1, Math.floor((ch - textY0) / 12) + 1);
    lines.slice(0, Math.min(3, maxLines)).forEach((line, i) => {
      inner.append('text')
        .attr('x', 4).attr('y', textY0 + i * 12)
        .attr('class', i === 0 ? '' : 'sub')
        .text(line);
    });

    if (ch > 32 && cw > 50) {
      inner.append('text')
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

  // Default view shows only backbone components; tests are excluded unless the
  // "Show tests" toggle is on (they remain in the manifest either way).
  const visible = physicals.filter(pc => pc.kind !== 'test' || showTests);

  const maxLoc = d3.max(visible, d => d.size_loc) || 1;
  // Group order includes 'unclassified' (neutral grey) so unowned files have a
  // home that is NOT core.
  const groups = ['core', 'seam', 'removable', 'unclassified'].map(cls => ({
    id: cls,
    name: cls,
    children: visible
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

  const physicalIds = visible.map(pc => pc.id);

  // Pressure-mode colour ramp domain (max pressure across visible tiles).
  const maxPressure = d3.max(visible, pc => {
    const k = pressureKeyFor(pc.metrics);
    return k ? pc.metrics[k] : 0;
  }) || 0;

  // Fill resolver: role mode = classification hue; pressure mode = heat ramp
  // (tiles with no pressure value go neutral so they recede). Time-shares the
  // channel — role hue is never permanently overwritten.
  function physicalFill(pc, primaryCls) {
    if (colorMode === 'pressure') {
      const k = pressureKeyFor(pc.metrics);
      if (!k) return '#262626';
      return pressureColor(pc.metrics[k], maxPressure);
    }
    return clsColor(primaryCls, pc.size_loc, maxLoc);
  }

  // Hatch pattern for de-emphasised test tiles (muted, striped, subordinate).
  ensureTestHatchPattern();

  // Shared defs for per-cell clip paths (overlap fix).
  const cellDefs = svg.append('defs').attr('class', 'cell-clip-defs');

  const g = svg.selectAll('g.node')
    .data(leaves)
    .join('g')
    .attr('class', 'node')
    .attr('data-pcid', d => d.data.id)
    .attr('transform', d => `translate(${d.x0},${d.y0})`);

  g.append('rect')
    .attr('class', 'pcell-rect')
    .attr('width',  d => Math.max(0, d.x1 - d.x0))
    .attr('height', d => Math.max(0, d.y1 - d.y0))
    .attr('fill',   d => physicalFill(d.data, d.data._primaryCls))
    .attr('rx', 2)
    .attr('opacity', d => d.data.kind === 'test' ? 0.4 : 1)
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

  // Hatch overlay on test tiles — drawn on top of the muted fill so tests read
  // as visually subordinate (striped) without being mistaken for backbone.
  g.filter(d => d.data.kind === 'test')
    .append('rect')
    .attr('class', 'test-hatch-overlay')
    .attr('width',  d => Math.max(0, d.x1 - d.x0))
    .attr('height', d => Math.max(0, d.y1 - d.y0))
    .attr('rx', 2)
    .attr('fill', 'url(#test-hatch)')
    .attr('pointer-events', 'none');

  // Cycle-membership marker: a dashed inner border on any tile in an import
  // cycle (Acyclic Dependencies Principle). Deterministic — read straight from
  // the per-component cycle_id/cycle_size the extractor stamped.
  g.filter(d => d.data.cycle_id)
    .append('rect')
    .attr('class', 'cycle-marker')
    .attr('x', 1.5).attr('y', 1.5)
    .attr('width',  d => Math.max(0, d.x1 - d.x0 - 3))
    .attr('height', d => Math.max(0, d.y1 - d.y0 - 3))
    .attr('rx', 2)
    .attr('fill', 'none')
    .attr('stroke', '#e08a3c')
    .attr('stroke-width', 1.5)
    .attr('stroke-dasharray', '3 2')
    .attr('pointer-events', 'none');

  // Per-cell labels, owner pills, metric badges (all clipped to the tile box).
  g.each(function(d) {
    const cw = d.x1 - d.x0;
    const ch = d.y1 - d.y0;
    // Raised gate: a tile narrower than ~46px cannot hold a filename; skip the
    // label rather than overflow it past the tile edge.
    if (cw < 46 || ch < 14) return;

    const gEl = d3.select(this);
    const metrics = d.data.metrics;

    // Clip everything in this cell to its own box (load-bearing overlap fix).
    const clipUrl = ensureCellClip(cellDefs, 'phys', d.data.id, cw, ch);
    const inner = gEl.append('g').attr('clip-path', clipUrl);

    // Metric badges at top-left; returns the reserved bottom of the badge row.
    const badgeBottom = drawMetricBadges(inner, metrics, cw, ch);
    const hasBadges = badgeBottom > 3;

    // Badges drawn but tile too short to clear them: badges only, no name.
    if (hasBadges && ch < 28) return;

    const textY0 = hasBadges ? badgeBottom + 8 : 13;
    const base = d.data.path.split('/').pop();
    const maxChars = Math.floor((cw - 8) / PX_PER_CHAR);
    const lines = wrapText(base, maxChars);
    const maxLines = Math.max(1, Math.floor((ch - textY0) / 11));
    lines.slice(0, Math.min(2, maxLines)).forEach((line, i) => {
      inner.append('text')
        .attr('x', 4).attr('y', textY0 + i * 11)
        .attr('class', i === 0 ? '' : 'sub')
        .text(line);
    });

    // Owner pills at bottom. Each owner label is truncated to the pill width and
    // the row is clipped by the cell clipPath, so a long owner name can never
    // bleed into a neighbour tile.
    if (ch > 36 && d.data.logical_owners && d.data.logical_owners.length > 0) {
      const owners = d.data.logical_owners;
      const pillH = 10;
      const pillPad = 2;
      let px = 3;
      const py = ch - pillH - 3;

      for (const oid of owners) {
        const lc = logicalById.get(oid);
        if (!lc) continue;
        const label = truncateToChars(lc.name, 10);
        const pw = label.length * 5.5 + 6;
        if (px + pw > cw - 3) break;  // overflow guard (clip backstops it too)

        const pillColor = d3.color(CLS_BASE[lc.classification] || '#444');
        pillColor.opacity = 0.7;

        inner.append('rect').attr('class', 'badge-stripe')
          .attr('x', px).attr('y', py)
          .attr('width', pw).attr('height', pillH)
          .attr('rx', 2)
          .attr('fill', pillColor.toString());

        inner.append('text').attr('class', 'badge-stripe')
          .attr('x', px + 3).attr('y', py + 8)
          .attr('font-size', 7).attr('fill', '#fff')
          .text(label);

        px += pw + pillPad;
      }
    }

    if (ch > 50 && cw > 50) {
      const offsetFromBottom = (d.data.logical_owners && d.data.logical_owners.length > 0) ? 16 : 4;
      inner.append('text').attr('class', 'sub')
        .attr('x', 4).attr('y', ch - offsetFromBottom)
        .text(d.data.size_loc + ' loc');
    }
  });

  // The "look here first" refactor-pressure side list (always present in the
  // Physical view; greyed with an install hint when radon was absent).
  renderPressureLens(visible, maxPressure);

  // Draw any pre-existing selection (e.g. after a resize or view switch back)
  buildEdgeOverlaySvg(physicalIds, physicalNodeCenters);
}

// ── Refactor-pressure "look here first" lens ──────────────────────────────────
// A collapsible side list ranking the top-N components by refactor pressure (the
// real metric, or the clearly-labelled LOC proxy when radon was absent). The
// honesty rule lives in the UI copy: "look here first", never "refactor this".
// Clicking an item highlights the tile and decomposes the factors.
//
// Collapse behaviour reuses the header.addEventListener('click') idiom from
// initGlossaryPanel(). Default: collapsed, so the list does not block tiles on
// first open. State is persisted in localStorage.
const LENS_STORAGE_KEY = 'architecture-treemap-lens-open';

function renderPressureLens(visible, maxPressure) {
  const host = document.getElementById('pressure-lens');
  if (!host) return;

  const usingProxy = !visible.some(pc => pc.metrics && pc.metrics.refactor_pressure !== undefined)
    && visible.some(pc => pc.metrics && pc.metrics.refactor_pressure_loc_proxy !== undefined);

  const lensHeadEntry = {
    what: 'Ranks attention, not verdicts.',
    how: 'Composite score: churn × complexity × fan-in, penalised by test coverage. High score = dense, volatile, load-bearing file.',
    teaches: 'A hot tile is an entry-point into investigation, not an instruction to refactor. Click a row to see WHY it is hot (the factor decomposition) and judge whether the heat is real.',
  };

  // Restore persisted collapse state (default: collapsed so the lens does not
  // block tiles on first visit).
  const isOpen = localStorage.getItem(LENS_STORAGE_KEY) === 'true';
  const toggleGlyph = isOpen ? '▾' : '▸';

  // The header row is always rendered; the body panel is toggled beneath it.
  let html = `
    <div id="lens-collapse-header" style="display:flex;align-items:center;gap:4px;cursor:pointer;user-select:none">
      <span class="lens-head" style="margin-bottom:0;flex:1">Look here first ${infoIcon(lensHeadEntry)}</span>
      <span id="lens-toggle-glyph" style="color:#888;font-size:10px;flex-shrink:0">${toggleGlyph}</span>
    </div>
    <div id="lens-body" style="display:${isOpen ? 'block' : 'none'}">`;

  if (pressureUnavailable() && !usingProxy) {
    // No pressure at all and radon flagged unavailable: grey the lens with a
    // reason, never an empty list that reads as "nothing to see".
    host.classList.add('lens-disabled');
    html += `<div class="lens-disabled-msg">Install <code>radon</code> to enable
      complexity &amp; refactor-pressure.<br>
      (<code>pip install radon</code>, then re-run /components-extract.)</div>`;
    html += '</div>';
    host.innerHTML = html;
    attachLensCollapseHandler(host);
    return;
  }

  host.classList.remove('lens-disabled');
  if (usingProxy) {
    html += `<div class="lens-proxy-note">Degraded LOC proxy (radon absent):
      LOC stands in for complexity. Install radon for the true score.</div>`;
  }

  const ranked = topByPressure(visible, 12);
  if (ranked.length === 0) {
    html += `<div class="lens-disabled-msg">No component has a computable
      refactor pressure (needs churn + complexity/LOC + fan-in).</div>`;
    html += '</div>';
    host.innerHTML = html;
    attachLensCollapseHandler(host);
    return;
  }

  for (const item of ranked) {
    const name = item.pc.path.split('/').pop();
    const t = maxPressure > 0 ? Math.min(1, item.value / maxPressure) : 0;
    const barColor = pressureColor(item.value, maxPressure);
    html += `<div class="lens-row" onclick="highlightPressureTile('${item.pc.id}')" title="Click to highlight + decompose">
      <div class="lens-row-top">
        <span class="lens-name">${escHtml(name)}</span>
        <span class="lens-val">${Math.round(item.value)}</span>
      </div>
      <div class="lens-bar"><div class="lens-bar-fill" style="width:${Math.round(t * 100)}%;background:${barColor}"></div></div>
      <div class="lens-factors" id="lens-factors-${escHtml(item.pc.id)}" style="display:none"></div>
    </div>`;
  }

  html += '</div>';  // close #lens-body
  host.innerHTML = html;
  attachLensCollapseHandler(host);
}

// Attach the collapse click handler after the lens HTML has been injected.
// Called each time renderPressureLens() rebuilds the DOM.
function attachLensCollapseHandler(host) {
  const header = host.querySelector('#lens-collapse-header');
  const body   = host.querySelector('#lens-body');
  const glyph  = host.querySelector('#lens-toggle-glyph');
  if (!header || !body) return;
  header.addEventListener('click', (e) => {
    // Do not collapse if the click was on a (?) icon — let the popover fire.
    if (e.target.closest('.ped-info')) return;
    const nowOpen = body.style.display !== 'block';
    body.style.display = nowOpen ? 'block' : 'none';
    if (glyph) glyph.textContent = nowOpen ? '▾' : '▸';
    localStorage.setItem(LENS_STORAGE_KEY, String(nowOpen));
  });
}

// Highlight a tile from the pressure list and decompose its pressure factors
// (churn x complexity x fan_in / test_ratio) so the user sees WHY it is hot.
function highlightPressureTile(pcid) {
  // Outline the tile.
  svg.selectAll('rect.pcell-rect').attr('stroke', null).attr('stroke-width', null);
  const cell = svg.select(`g.node[data-pcid="${pcid}"] rect.pcell-rect`);
  if (!cell.empty()) {
    cell.attr('stroke', '#fff').attr('stroke-width', 2.5);
  }
  selectedNodeIds = new Set([pcid]);
  buildEdgeOverlaySvg((COMPONENTS_DATA.physical_components || []).map(c => c.id), physicalNodeCenters);
  showPanel(pcid, 'physical');

  // Toggle the inline factor decomposition under the clicked row.
  const slot = document.getElementById(`lens-factors-${pcid}`);
  if (!slot) return;
  const pc = physicalById.get(pcid);
  const m = (pc && pc.metrics) || {};
  const usingProxy = m.refactor_pressure === undefined && m.refactor_pressure_loc_proxy !== undefined;
  const complexity = usingProxy ? `LOC/100 = ${((m.loc || 0) / 100).toFixed(1)} (proxy)` : (m.cyclomatic !== undefined ? m.cyclomatic : 'n/a');
  const factor = (label, val) => `<div class="lens-factor"><span>${label}</span><span>${val}</span></div>`;
  slot.innerHTML =
    factor('churn (90d)', m.churn_90d !== undefined ? m.churn_90d : 'n/a') +
    factor('complexity', complexity) +
    factor('fan-in', m.fan_in !== undefined ? m.fan_in : 'n/a') +
    factor('test ratio', m.test_ratio !== undefined ? m.test_ratio : '0 (untested → amplifies)');
  slot.style.display = slot.style.display === 'none' ? 'block' : 'none';
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
  // Identify ambient/cross-cutting sinks so their many incoming edges can be
  // drawn faintly instead of forming a hairball. Stored at module scope so
  // refreshGraphEdges (on click) reuses the same set.
  graphAmbientIds = ambientNodeIds(Array.from(viewIds), activeEdges);

  const maxLoc = d3.max(logicals, d => d.size_estimate_loc) || 1;

  // Node radius = sqrt of normalised LOC, scaled to range [8, 32]
  const rScale = d3.scaleSqrt()
    .domain([0, maxLoc])
    .range([8, 32]);

  // Seed positions by classification zone. The seeded PRNG (from the URL hash)
  // drives BOTH the base jitter here and the perturbation below, so the layout
  // is reproducible per seed — same hash gives the same starting positions.
  // Nodes with a stored position (from a prior drag this session) skip the RNG
  // so the RNG call sequence for unseen nodes is unchanged — determinism is
  // preserved for any first render where graphNodePositions is empty.
  const rng = mulberry32(getSeedFromHash());
  const nodes = logicals.map(c => {
    const stored = graphNodePositions.get(c.id);
    // Consume two RNG values regardless of whether the stored position is used,
    // so the sequence for later nodes is identical to the no-stored-positions case.
    const rx = rng() * 60 - 30;
    const ry = rng() * 60 - 30;
    const x = stored ? stored.x : CLUSTER_X[c.classification] * w + rx;
    const y = stored ? stored.y : CLUSTER_Y * h + ry;
    return {
      ...c,
      x,
      y,
      // Pin nodes whose positions were restored so the simulation doesn't move them.
      fx: stored ? stored.x : undefined,
      fy: stored ? stored.y : undefined,
      r: rScale(c.size_estimate_loc),
    };
  });

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
    .attr('stroke-opacity', e => graphEdgeOpacity(e, graphAmbientIds))
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
        // Keep pinned; double-click to unpin. Persist final position so that
        // re-entering the Graph view restores this arrangement.
        graphNodePositions.set(d.id, { x: d.fx, y: d.fy });
      })
    )
    .on('dblclick', (event, d) => {
      d.fx = null;
      d.fy = null;
      // Remove stored position so this node gets seed-layout on next render.
      graphNodePositions.delete(d.id);
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
  // Nodes with stored (restored) positions are already pinned; consume the RNG
  // values for them anyway to keep the sequence identical to a fresh render.
  for (const n of nodes) {
    const px = (rng() - 0.5) * 40;
    const py = (rng() - 0.5) * 40;
    if (n.fx === undefined) {
      n.x += px;
      n.y += py;
    }
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

// Opacity for a graph edge, honouring selection, the show-all toggle, and the
// ambient-sink damping. An edge into an ambient node is drawn faint by default
// (it still lights to full opacity when its endpoint is selected), so the
// cross-cutting relationship is visible-on-demand, not a permanent hairball.
function graphEdgeOpacity(e, ambientIds) {
  if (selectedNodeIds.size > 0) {
    const active = selectedNodeIds.has(e.from) || selectedNodeIds.has(e.to);
    return active ? 0.9 : 0.15;
  }
  if (ambientIds && ambientIds.has(e.to)) return 0.12;  // ambient sink: faint
  return showAllEdges ? 0.3 : 0.7;
}

function refreshGraphEdges(edgePaths) {
  edgePaths.attr('stroke-opacity', e => graphEdgeOpacity(e, graphAmbientIds));
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

  // Epistemic spine: a structural/LLM pass MUST NOT auto-assert removable. A
  // removable verdict with no grounding slice/finding/prune renders as a
  // PROVISIONAL "requires your judgment" marker rather than a confident red
  // badge — removable is the one verdict meaningless without intent.
  const clsEntry = (PEDAGOGY.classifications || {})[lc.classification];
  const provisionalRemovable = (lc.classification === 'removable' && !hasEvidence);

  let classHtml;
  if (provisionalRemovable) {
    classHtml = `<span class="badge provisional-removable" title="${escHtml((clsEntry && clsEntry.teaches) || '')}">requires your judgment — run /audit-slice</span>`;
  } else {
    classHtml = `<span class="badge ${lc.classification}" title="${escHtml(explainerTitle(clsEntry))}">${lc.classification}</span>`;
  }
  // Epistemic-source chip on the classification (manifest value or pedagogy default).
  const epiSource = lc.epistemic_source || (clsEntry ? clsEntry.epistemic_source : null);
  const epiHtml = epiSource ? epistemicChip(epiSource) : '';
  const judgementPill = hasEvidence ? '' : '<span class="badge" style="background:#222;color:#888;border-color:#444;margin-left:4px">judgement only</span>';

  document.querySelectorAll('#panel-header .badge, #panel-header .epi-chip').forEach(b => b.remove());
  panelTitle.insertAdjacentHTML('afterend', classHtml + epiHtml + judgementPill);

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
  const clsEntry = (PEDAGOGY.classifications || {})[cls];
  const classHtml = `<span class="badge ${cls}" title="${escHtml(explainerTitle(clsEntry))}">${cls}</span>`;
  const epiHtml = clsEntry ? epistemicChip(clsEntry.epistemic_source) : '';
  document.querySelectorAll('#panel-header .badge, #panel-header .epi-chip').forEach(b => b.remove());
  panelTitle.insertAdjacentHTML('afterend', classHtml + epiHtml);

  let html = '';
  html += `<div class="panel-section"><h3>Path</h3><p style="font-family:monospace;color:#7ab8e0;word-break:break-all">${escHtml(pc.path)}</p></div>`;
  html += `<div class="panel-section"><h3>Size</h3><p>${pc.size_loc} loc &nbsp;·&nbsp; kind: ${pc.kind}</p></div>`;

  // Import-cycle membership (deterministic, Acyclic Dependencies Principle).
  if (pc.cycle_id) {
    const cycleEntry = {
      what: 'This file is in an import cycle.',
      how: 'Detected by Tarjan SCC over the static import graph.',
      teaches: 'Every member can reach every other, so none can be changed in isolation. Acyclic Dependencies Principle: break the cycle by having one member depend on an abstraction instead.',
      epistemic_source: 'measured',
    };
    html += `<div class="panel-section"><h3>Import cycle</h3>
      <p style="color:#e08a3c">In a cycle of ${pc.cycle_size} modules
        ${infoIcon(cycleEntry)}
      </p></div>`;
  }

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

// Side-panel metric breakdown — descriptor-driven. Shows every metric present,
// in registry order, labelled and unit-suffixed from the descriptor. Any metric
// present but not described (forward-compat) is shown with its raw key.
function metricsBreakdownHtml(metrics) {
  const rows = [];
  const seen = new Set();
  const metricExplainers = PEDAGOGY.metrics || {};
  const emit = (field) => {
    if (metrics[field] === undefined || seen.has(field)) return;
    seen.add(field);
    const desc = METRIC_DESCRIPTORS[field];
    const label = desc ? desc.label : field;
    const unit  = desc && desc.unit ? ' ' + desc.unit : '';
    const help  = infoIcon(metricExplainers[field]);  // (?) on every metric row
    rows.push(`<div class="file-row"><span style="color:#aaa">${escHtml(label)} ${help}</span><span style="color:#ddd">${metrics[field]}${escHtml(unit)}</span></div>`);
  };
  for (const field of Object.keys(METRIC_DESCRIPTORS)) emit(field);
  for (const field of Object.keys(metrics)) emit(field);  // undescribed metrics
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

  // Edge-type help comes from the manifest pedagogy block (single source of
  // truth), NOT a hardcoded JS table. Same for the epistemic-source chip.
  const etypeEntry = (PEDAGOGY.edge_types || {})[e.type];
  const epiSource = e.epistemic_source
    || (etypeEntry ? etypeEntry.epistemic_source : null);

  let html = `<div class="edge-row ${e.type}">`;
  html += `<span class="peer" style="cursor:pointer" onclick="showPanel('${peer}','${panelMode}')">${arrow} ${escHtml(peerNm)}</span>`;
  html += `<span class="etype" title="${escHtml(explainerTitle(etypeEntry))}" style="cursor:help">${e.type}</span>`;
  html += infoIcon(etypeEntry);
  if (epiSource) html += epistemicChip(epiSource);
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

// ── Pedagogy-driven legend + glossary (no hardcoded help) ─────────────────────
//
// All explainer content is read from the manifest's `pedagogy` block (built by
// scripts/pedagogy_registry.py). The old hardcoded EDGE_TYPE_HELP table and the
// classification/evidence help that lived in the template HTML are GONE — that
// duplication was the Connascence-of-Value smear the metrics ADR killed for
// metrics; this kills it for pedagogy too.

// Populate the legend's classification / edge / evidence / direction rows and
// the glossary entirely from the pedagogy block. The template ships only empty
// host containers; this fills them so the render stays a pure function of the
// manifest.
function buildLegendPopovers() {
  buildClassificationLegend();
  buildEdgeTypeLegend();
  buildEvidenceLegend();
  buildDirectionLegend();
  buildEpistemicLegend();
  buildGlossary();
}

function legendRow(swatchColor, label, entry, extraMark) {
  const sw = swatchColor
    ? `<div class="leg-swatch" style="background:${swatchColor}"></div>`
    : '';
  const mark = extraMark ? ` ${extraMark}` : '';
  return `<div class="leg-row">${sw}<span>${escHtml(label)}${mark}</span>${infoIcon(entry)}</div>`;
}

function buildClassificationLegend() {
  const host = document.getElementById('leg-classifications');
  if (!host) return;
  const cls = PEDAGOGY.classifications || {};
  let html = '';
  for (const key of ['core', 'seam', 'removable']) {
    const entry = cls[key];
    const color = (entry && entry.color) || CLS_BASE[key] || '#444';
    html += legendRow(color, key, entry);
  }
  html += legendRow(CLS_BASE.unclassified, 'unclassified',
    { what: 'A physical file with no logical owner — not grounded by any synthesis judgment.',
      how: 'Rendered neutral grey; it does NOT inherit core.',
      teaches: 'Unowned != core. Run /synthesize-audit to classify it.',
      epistemic_source: 'measured' });
  host.innerHTML = html;
}

function buildEdgeTypeLegend() {
  const host = document.getElementById('leg-edge-types');
  if (!host) return;
  const et = PEDAGOGY.edge_types || {};
  let html = '';
  const marks = { 'shared-state': '⚠', 'background-knowledge': '★' };
  for (const key of ['direct-call', 'shared-state', 'background-knowledge']) {
    const entry = et[key];
    const color = (entry && entry.color) || ETYPE_COLOR[key] || '#888';
    html += legendRow(color, key, entry, marks[key] || '');
  }
  host.innerHTML = html;
}

function buildEvidenceLegend() {
  const host = document.getElementById('leg-evidence');
  if (!host) return;
  const ec = PEDAGOGY.evidence_classes || {};
  const dashRow = (label, dash, entry) =>
    `<div class="leg-row"><svg width="20" height="10"><line x1="0" y1="5" x2="20" y2="5"
       stroke="#888" stroke-width="1.5"${dash ? ` stroke-dasharray="${dash}"` : ''}/></svg>
       <span>${escHtml(label)}</span>${infoIcon(entry)}</div>`;
  let html = '';
  html += dashRow('static (AST-detected)', null, ec['static']);
  html += dashRow('audit-asserted', '4 4', ec['audit-asserted']);
  host.innerHTML = html;
}

function buildDirectionLegend() {
  const host = document.getElementById('leg-direction');
  if (!host) return;
  const entry = PEDAGOGY.direction_encoding;
  host.innerHTML =
    `<div style="font-size:9px;color:#777">Direction: outgoing = full colour;
       incoming = lighter + lower opacity ${infoIcon(entry)}</div>`;
}

function buildEpistemicLegend() {
  const host = document.getElementById('leg-epistemic');
  if (!host) return;
  const sources = PEDAGOGY.epistemic_sources || {};
  let html = '';
  for (const key of ['measured', 'metric-anchored', 'requires-your-intent']) {
    const meta = sources[key];
    if (!meta) continue;
    // Construct a popover-compatible entry from the epistemic-source descriptor.
    const popEntry = { what: meta.what, how: meta.how || null, teaches: meta.teaches || null };
    html += `<div class="leg-row">
      <span class="epi-chip" style="border-color:${meta.color};color:${meta.color}">${escHtml(meta.label)}</span>
      ${infoIcon(popEntry)}
    </div>`;
  }
  host.innerHTML = html;
}

// The glossary: each concept tethered to the proxy metric that surfaces it in
// THIS tool, with the honesty note where the proxy is weak. Content from the
// pedagogy block, never hardcoded prose.
function buildGlossary() {
  const host = document.getElementById('glossary-body');
  if (!host) return;
  const glossary = PEDAGOGY.glossary || [];
  if (glossary.length === 0) { host.innerHTML = '<p style="color:#555">No glossary in this manifest.</p>'; return; }
  let html = '';
  for (const g of glossary) {
    const anchor = g.anchor
      ? `<span class="gloss-anchor" title="Proxy metric / element in this tool">↳ ${escHtml(g.anchor)}</span>`
      : '';
    html += `<div class="gloss-entry">
      <div class="gloss-term">${escHtml(g.term)} ${anchor}</div>
      <div class="gloss-def">${escHtml(g.definition || '')}</div>`;
    if (g.proxy_honesty) {
      html += `<div class="gloss-honesty">Proxy honesty: ${escHtml(g.proxy_honesty)}</div>`;
    }
    html += `</div>`;
  }
  host.innerHTML = html;
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

// Glossary collapsible (replaces the old hardcoded "How to read this" panel;
// content now comes from the pedagogy block via buildGlossary()).
function initGlossaryPanel() {
  const panel  = document.getElementById('glossary-panel');
  const toggle = document.getElementById('glossary-toggle');
  const header = document.getElementById('glossary-header');
  if (!panel || !toggle || !header) return;

  const storageKey = 'architecture-treemap-glossary-open';
  const stored = localStorage.getItem(storageKey);
  const isOpen = stored === 'true';  // closed by default on first visit
  panel.style.display = isOpen ? 'block' : 'none';
  toggle.textContent  = isOpen ? '▾' : '▸';

  header.addEventListener('click', () => {
    const nowOpen = panel.style.display !== 'block';
    panel.style.display = nowOpen ? 'block' : 'none';
    toggle.textContent  = nowOpen ? '▾' : '▸';
    localStorage.setItem(storageKey, String(nowOpen));
  });
}

// ── View toggle ───────────────────────────────────────────────────────────────

function setView(view) {
  currentView = view;
  selectedNodeIds = new Set();

  document.getElementById('btn-logical').classList.toggle('active', view === 'logical');
  document.getElementById('btn-physical').classList.toggle('active', view === 'physical');
  document.getElementById('btn-graph').classList.toggle('active', view === 'graph');

  // Pressure lens + colour-mode toggle are Physical-view affordances only.
  const lens = document.getElementById('pressure-lens');
  const colorToggle = document.getElementById('color-mode-toggle');
  if (lens) lens.style.display = (view === 'physical') ? 'block' : 'none';
  if (colorToggle) colorToggle.style.display = (view === 'physical') ? 'flex' : 'none';

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

// ── "Show tests" toggle (Physical view) ───────────────────────────────────────

function toggleTests(checked) {
  showTests = checked;
  localStorage.setItem(SHOW_TESTS_STORAGE_KEY, String(checked));
  if (currentView === 'physical') renderPhysical();
}

// ── Physical colour-mode toggle: role | pressure ──────────────────────────────
// Time-shares the colour channel; never permanently overrides the role hue.
function setColorMode(mode) {
  colorMode = (mode === 'pressure') ? 'pressure' : 'role';
  localStorage.setItem(COLOR_MODE_STORAGE_KEY, colorMode);
  const roleBtn = document.getElementById('btn-color-role');
  const pressBtn = document.getElementById('btn-color-pressure');
  if (roleBtn) roleBtn.classList.toggle('active', colorMode === 'role');
  if (pressBtn) pressBtn.classList.toggle('active', colorMode === 'pressure');
  if (currentView === 'physical') renderPhysical();
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

// ── Explainer popover event wiring ────────────────────────────────────────────
// Single delegated handler on the document: intercepts clicks on any .ped-info
// span that carries a data-ped-entry attribute (set by infoIcon()). Only one
// popover is open at a time; a second click on the same icon closes it.

document.addEventListener('click', function(e) {
  const icon = e.target.closest('.ped-info[data-ped-entry]');
  if (icon) {
    e.stopPropagation();
    const pop = document.getElementById('ped-popover');
    // Toggle: if this icon is already the anchor, close; else open with new content.
    if (pop && pop._anchorEl === icon && pop.style.display === 'block') {
      closeExplainerPopover();
    } else {
      let entry = {};
      try { entry = JSON.parse(icon.dataset.pedEntry); } catch (_) { /* malformed — show empty */ }
      openExplainerPopover(icon, entry);
    }
    return;
  }

  // Close button inside the popover.
  if (e.target.closest('#ped-popover-close')) {
    closeExplainerPopover();
    return;
  }

  // Click anywhere outside the popover closes it.
  const pop = document.getElementById('ped-popover');
  if (pop && pop.style.display === 'block' && !e.target.closest('#ped-popover')) {
    closeExplainerPopover();
  }
});

// Keyboard dismiss.
document.addEventListener('keydown', function(e) {
  if (e.key === 'Escape') closeExplainerPopover();
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
initGlossaryPanel();

// Reflect persisted "Show tests" state into the checkbox on load.
const _showTestsCheckbox = document.getElementById('show-tests');
if (_showTestsCheckbox) _showTestsCheckbox.checked = showTests;

// Reflect persisted colour-mode into the toggle buttons (state only; render
// happens on first setView/renderLogical below — Physical view applies it).
const _roleBtn = document.getElementById('btn-color-role');
const _pressBtn = document.getElementById('btn-color-pressure');
if (_roleBtn) _roleBtn.classList.toggle('active', colorMode === 'role');
if (_pressBtn) _pressBtn.classList.toggle('active', colorMode === 'pressure');

// Pressure lens + colour toggle start hidden (default view is Logical).
const _lens = document.getElementById('pressure-lens');
if (_lens) _lens.style.display = 'none';
const _colorToggle = document.getElementById('color-mode-toggle');
if (_colorToggle) _colorToggle.style.display = 'none';

renderLogical();
