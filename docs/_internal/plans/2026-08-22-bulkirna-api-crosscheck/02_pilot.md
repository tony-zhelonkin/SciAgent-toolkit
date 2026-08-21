# Phase 02 — rewrite three skills as the pilot

**Repo:** scio · **Blocked by:** phase 01's floor decision · Read `00_INDEX.md`
§2 for the principle.

## The three

The subset the RNAseq-toolkit agent named, and the highest-confidence rows in the
inventory:

1. **`coresh-signature-search`** — ships four scripts under `scripts/`, all
   superseded: `coresh_batch.R` → `coresh_search`, `extract_gene_loadings.R` →
   `coresh_loadings`/`gsdb_coresh`, `symbols_to_entrez.R` →
   `gene_to_entrez`/`entrez_to_gene`, `validate_coresh_install.R` →
   `coresh_validate`. The SKILL.md names them at nine places.
2. **`bulk-rnaseq-gsea`** — four reference files describing a wrapper around
   `clusterProfiler::GSEA()` plus `01_scripts/RNAseq-toolkit/...` paths. The
   `gs_*` and `gsdb_*` layers cover execution, master tables, plots and database
   registration.
3. **`annotate-bulk-rnaseq-data`** — `SKILL.md:71` pins **RNAseq-toolkit v0.2.0**
   by name and lists four script files, all now exports.

## What to keep, explicitly

Deleting the wrong half destroys the asset. **Keep every trap.** Named examples
that must survive verbatim in substance:

- the **CORESH Entrez-integer trap** and the assertion that catches it
  (`coresh-signature-search` ~line 193) — the package cannot warn you that your
  symbols silently matched nothing;
- **species mismatch** and the ortholog remedy (~199);
- **`eps = 0`** or p-values truncate at ~1e-4 (`bulk-rnaseq-gsea:173`);
- the **`p.adjust` vs `qvalue`** divergence between GSEA runs
  (`references/master-tables.md:396`) — unless phase 01 shows `gs_stat_types`
  resolves it, in which case route to that and say so;
- timing figures that set expectations (~10–20 s variance-only, 2–5 min with
  p-values).

## Shape after the rewrite

*When to reach for this, why, the traps, and which verb does it.* Not how to
implement it. A skill that used to teach a procedure becomes a skill that teaches
a judgement and names a function.

Frontmatter `description` stays ≤ 350 chars and must stop promising what the file
no longer contains — `bulk-rnaseq-gsea`'s currently advertises "MSigDB execution
via clusterProfiler/fgsea".

## The deletion that needs a check first

Removing `coresh-signature-search/scripts/*.R` changes **mount shape**, and those
paths are live in every bound project as
`.claude/skills/coresh-signature-search/scripts/…`. Before deleting, grep the
fleet for anything that sources them:

```bash
grep -rn "coresh_batch\|extract_gene_loadings\|symbols_to_entrez" \
  /data1/users/antonz /data2/users/JCRLab --include="*.R" --include="*.py" \
  --include="*.qmd" --include="*.md" 2>/dev/null | grep -v SciAgent-toolkit
```

A hit means a real analysis depends on the path. Then the scripts stay as thin
shims that call the package, or they stay put and only the prose changes — decide
on the evidence, and say which in the report.

## If the floor is the installed API

Phase 01 may find the absorbing functions absent from the live image. Then the
skill cannot simply say `coresh_search(...)`. Options, in preference order:

1. Say the version the verb requires, and give the one-line check the reader can
   run: `"coresh_search" %in% getNamespaceExports("bulkiRNA")`. **Probe, do not
   assert** — the same rule `delegate-cli` now follows for codex flags.
2. Keep both paths, package-first, with the script path named as the fallback for
   older images. Costs length, and length is what this plan is trying to reduce.
3. Wait for the image pin. Then this phase blocks on a scbio-docker change and
   should say so rather than shipping text nobody can run.

Do not invent a version-detection helper for scio to ship. That is a fourth
delivery mechanism for one package's problem.

## Gates

```bash
bash tests/run-all.sh                     # 65 passing, 0 failing
bin/scio lint --check toolkit --strict    # clean; descriptions ≤ 350 chars
bash tests/test_mount_layout_identity.sh  # if any scripts/ file moved
bash tests/test_skill_link_integrity.sh
```

## Do not

- Do not write or change bulkiRNA code, and do not touch `BULKIRNA_SHA`.
- Do not touch the five TE skills or the three ATAC skills.
- Do not delete a trap because the API "probably handles it". Verify or keep.
- Do not run the fleet sweep. But do finish before it — see `00_INDEX.md` §5.
