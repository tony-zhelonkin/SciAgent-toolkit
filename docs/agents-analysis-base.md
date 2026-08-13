# Base analysis agents

The eight agents that the `base` bioinformatics role binds. They cover the
day-to-day analysis loop: looking up tool documentation, interpreting
biological mechanisms, poking at data files, drafting figure legends, curating
repo docs, reviewing refactors, and writing session handoffs. None of these
agents are part of the architect design pipeline — they are on-demand helpers
invoked by the analyst at the keyboard.

| Agent | Purpose |
|-------|---------|
| `docs-librarian` | Find tool/package documentation via web search |
| `bio-interpreter` | Research biological mechanisms |
| `insight-explorer` | Skeptical exploration of data files |
| `captions` | Generate figure legends |
| `figure-audit` | Inspect rendered figures against the visual contract |
| `doc-curator` | Clean up repo documentation |
| `code-reviewer` | Review refactored code |
| `handoff` | Create session handoff docs |
