# Universal Project Scaffold

This is the single universal project template used by `scripts/init-project.sh`
to bootstrap a new analysis project. It is not meant to be used directly —
run `scripts/init-project.sh /path/to/new-project` from the scbio-docker repo
root and the script will copy this tree into place, substitute placeholders,
and wire up the devcontainer.

Layout mirrors what real analysis projects converge on:

- `00_data/{raw,processed,references}/` — data inputs and intermediates
- `01_modules/.ref/` — gitignored slot for reference codebases / submodules
- `02_analysis/{config,helpers,scripts,notebooks}/` — analysis code
- `03_results/{checkpoints,plots,tables}/` — outputs
- `docs/{plan,ai-generated,raw}/` — research docs, plan, AI notes
- `.devcontainer/`, `.vscode/` — container + editor config
