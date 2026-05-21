# architect/ commands

The 14 slash commands that drive the architect role's design-before-code
pipeline. They cover the standard cadence (`/map` → `/review` → `/synthesize`
→ `/design` → `/plan` → `/implement` → `/verify`) plus meta-level variants
for changes to the pipeline itself, the `/diagram` rendering helper, the
`/status` phase-doc renderer, and the `/architect` consistency gate.

For the canonical sequence and worked examples see
`../../docs/workflows/architect/00-quickstart.md`. For the rationale behind
the pipeline shape see `../../docs/workflows/architect/01-architecture.md`.
