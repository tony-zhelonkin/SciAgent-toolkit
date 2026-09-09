# Base analysis agents

The seven agents that the `base` bioinformatics role binds. They cover the
day-to-day analysis loop: looking up tool documentation, interpreting
biological mechanisms, poking at data files, drafting figure legends, curating
repo docs, and reviewing refactors. None of these agents are part of the
architect design pipeline — they are on-demand helpers invoked by the analyst
at the keyboard.

Recording the session is the `/handoff` command rather than an agent: the
reasoning worth keeping lives in the context that did the work, and a fresh
reader can only reconstruct it from artifacts.

| Agent | Purpose |
|-------|---------|
| `docs-librarian` | Find tool/package documentation via web search |
| `bio-interpreter` | Research biological mechanisms |
| `insight-explorer` | Skeptical exploration of data files |
| `captions` | Generate figure legends |
| `figure-audit` | Inspect rendered figures against the visual contract |
| `doc-curator` | Clean up repo documentation |
| `code-reviewer` | Review refactored code |
