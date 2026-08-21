Implemented Worker 2 only in [skills/delegate-cli/SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md).

The skill is now asset-led, reduced from 416 to 150 lines, and retains the required prompt, `agy`, implementer, and verification guidance. The description is 285 characters. No assets, tests, or `lib/scio/` files were changed. No commit was made.

The skill-creator guidance informed the concise 150-line structure.

Disagreement with the plan: the repository currently has 62 tests rather than the expected 63, and `skills/delegate-cli/assets/` is absent. I left that Worker 1 scope untouched. Strict lint also exits successfully despite those assets being absent.

Gate outputs:

`bash tests/run-all.sh` — exit 0

```text
PASS [test_block_drift.sh]
block_read: missing required <id> argument
block_read: missing required <id> argument
block_hash_check: missing required <id> argument
block_hash_check: missing required <id> argument
block_write: missing required <id> argument
block_remove: missing required <id> argument
block_line_range: missing required <id> argument
block_write: only SCIO:CRAFT is writable
PASS [test_block_id_required.sh]
PASS [test_block_marker_boundary.sh]
PASS [test_block_missing_marker.sh]
PASS [test_block_no_markers_append.sh]
PASS [test_block_preserves_mode.sh]
PASS [test_block_roundtrip.sh]
block_write: unsupported SCIO:CRAFT marker version v2
PASS [test_block_version_tolerance.sh]
PASS [test_block_write_io_failure.sh]
PASS [test_build_release_contents.sh]
PASS [test_build_release_determinism.sh]
PASS [test_build_release_refusals.sh]
PASS [test_collision_allowlist_blocks_unknown.sh]
PASS [test_collision_allowlist_passes_known.sh]
PASS [test_craft_block_budget.sh]
PASS [test_craft_block_roundtrip.sh]
PASS [test_craft_verb.sh]
  parity: all 13 contract names present in BOTH R and Python
  python: import + path/caption/master/round/cue execution OK (no heavy deps)
SKIP: matplotlib absent — plotting paths not executed (validated structurally only)
PASS [test_figure_helpers_contract.sh]
PASS [test_hook_body_propagation.sh]
PASS [test_install_atomic.sh]
PASS [test_install_coexist_uninstall.sh]
PASS [test_install_invariants.sh]
PASS [test_install_locality_precedence.sh]
PASS [test_install_verify_and_dryrun.sh]
  names: all 10 contract functions present
  python: imports with jscatter ABSENT + config/accessor/lifecycle execution OK
PASS [test_interactive_helpers_contract.sh]
PASS [test_link_analysis_helpers.sh]
PASS [test_link_gitignore.sh]
PASS [test_link_refuses_populated_directory.sh]
PASS [test_link_sweeps_legacy_mounts.sh]
PASS [test_link_verb.sh]
PASS [test_lint_comment_intent.sh]
PASS [test_lint_stage_layout.sh]
PASS [test_lint_stage_thinness.sh]
PASS [test_lint_verb.sh]
PASS [test_mount_layout_identity.sh]
PASS [test_no_duplicate_basenames.sh]
PASS [test_no_exit_in_libs.sh]
SKIP: Rscript absent — framework resolution not executed
PASS [test_peak_atlas_framework_resolution.sh]
PASS [test_plan_templates.sh]
PASS [test_skill_frontmatter_valid.sh]
PASS [test_skill_link_integrity.sh]
PASS [test_template_provenance.sh]
PASS [test_validate_collision_quiet.sh]
PASS [test_validate_collision_three_kinds.sh]
PASS [test_validate_collision_warn.sh]
PASS [test_validate_description_scalar_styles.sh]
PASS [test_validate_docs_layout.sh]
PASS [test_validate_figure_style.sh]
PASS [test_validate_frontmatter_shape.sh]
PASS [test_validate_hooks.sh]
PASS [test_validate_no_false_positives.sh]
PASS [test_validate_provenance.sh]
PASS [test_validate_results_layout.sh]
PASS [test_verb_help.sh]
PASS [test_verb_smoke_load_graph.sh]
[mllmct] $HOME cache not writable — using a temp UV_CACHE_DIR
[mllmct] smoke_check_versions + offline pytest in locked sandbox
smoke_check_versions: exact pins
  [PASS] mllmcelltype==2.0.7 — got 2.0.7
  [PASS] google-genai==2.6.0 — got 2.6.0
  [PASS] pydantic==2.13.3 — got 2.13.3
  [PASS] requests==2.33.1 — got 2.33.1
smoke_check_versions: monkeypatch seams
  [PASS] mllmcelltype.prompts.create_prompt callable
  [PASS] create_prompt accepts 'prompt_template'
  [PASS] interactive_consensus_annotation accepts 'prompt_template'
  [PASS] mllmcelltype.interactive_consensus_annotation
  [PASS] mllmcelltype.logger.setup_logging callable
  [PASS] mllmcelltype.providers.openrouter.requests IS requests
  [PASS] google.genai.models.Models.generate_content
  [PASS] GenerateContentConfig has 'temperature'
  [PASS] GenerateContentConfig has 'seed'
  [PASS] usage_metadata has 'prompt_token_count'
  [PASS] usage_metadata has 'candidates_token_count'
  [PASS] usage_metadata has 'total_token_count'

smoke_check_versions: ALL PASS — pins + seams intact.
.........................................................                [100%]
57 passed in 0.89s
[nfcore-rnaseq-execution] testing make_samplesheet.sh
  PASS: header
  PASS: total lines (header + 4 rows)
  PASS: distinct samples
  PASS: SAMPLEB row count (per-lane)
  PASS: fastq_1 absolute & R1
  PASS: fastq_2 absolute & R2
  PASS: strandedness arg
  PASS: default strandedness
  PASS: seq_platform header
  PASS: seq_platform value
  PASS: no seq_platform by default
  PASS: errors on missing R2
  PASS: errors on missing dir

[nfcore-rnaseq-execution] ALL TESTS PASS
SKIP: Rscript not found — peak-atlas tests need R + Bioconductor.
SKIP: Rscript not found — peak-atlas tests need R + Bioconductor.
SKIP: Rscript not found — peak-atlas tests need R + Bioconductor.
SKIP [te-gene-featurecounts] image 'te-fc:2.0.2' not built — run env/build.sh to enable smoke tests
[te-reference-saf-build] testing build_te_saf.sh on synthetic fixtures
  PASS: grouped SAF header columns
  PASS: grouped SAF row count (4 loci)
  PASS: L1Md_A:L1:LINE label + cols (1/100/200/+)
  PASS: B1_Mus1:Alu:SINE label + cols
  PASS: MER20:hAT:DNA label + cols
  PASS: unique groups (L1Md_A pooled -> 3 groups)
  PASS: retro-only whitelist drops DNA (2 groups)
  PASS: no-exon SAF header columns
  PASS: no residual TE interval overlaps exon 1:120-180
  PASS: groups preserved after subtract (3)
  PASS: no-exon Start is 1-based (min=100)

[te-reference-saf-build] ALL TESTS PASS

==== test summary ====
passed: 62
failed: 0
```

`bin/scio lint --check toolkit --strict` — exit 0

```text
```

`bash -n skills/delegate-cli/assets/*.sh` — exit 127

```text
bash: skills/delegate-cli/assets/*.sh: No such file or directory
```

`grep -n '\-\-search\|run_in_background' skills/delegate-cli/SKILL.md` — exit 1, expected for no matches

```text
```