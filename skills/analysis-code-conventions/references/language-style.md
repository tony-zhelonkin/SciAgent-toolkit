# Language Style

R and Python conventions at the level of the line.
The sibling references own structure; this one owns spelling, and the traps
that survive a competent model's defaults.

Nothing here restates PEP 8 or the tidyverse style guide.
Assume both. What follows is the house's departures from them,
and the constructs that are wrong often enough to name.

## Naming

- `snake_case` for variables and functions in both languages; verbs for functions.
- `SCREAMING_SNAKE_CASE` for constants.
- `PascalCase` for Python classes only.
- Leading underscore marks a Python internal, per `helper-family-apis`.

Two prohibitions carry weight in R:

- **No `camelCase`.** It reads as a foreign import in an R codebase.
- **No `dot.names`.** `my.func` collides with S3 dispatch, where the dot is
  syntax rather than decoration. This one is a correctness matter.

## Progress and failure in R

- `message()` for progress, so it goes to stderr and a caller can silence it.
  `print()` and `cat()` are for results a reader consumes.
- `warning()` when execution can continue on a degraded path.
- `stop()` when it cannot. A stage that carries on past a broken precondition
  produces an artifact nobody can trust.

Say what failed and what the caller should do about it. An error naming the
offending value is worth more than one naming the function.

## The traps

These are wrong in ways that pass review and fail at runtime.

```r
for (i in 1:length(vec))     # 1:0 iterates twice when vec is empty
for (i in seq_along(vec))    # correct

keep <- T                    # T is a rebindable variable
keep <- TRUE                 # a reserved literal

out <- c()
for (x in xs) out <- c(out, f(x))   # reallocates every iteration
out <- vapply(xs, f, numeric(1))    # pre-sized, and type-checked

attach(data)                 # binds an invisible scope a reader cannot see
result <<- value             # assigns somewhere the reader must guess
suppressWarnings(risky())    # discards the diagnostic with the noise
```

```python
def f(xs, acc=[]):           # one list, shared across every call
def f(xs, acc=None):         # allocate inside

except:                      # swallows KeyboardInterrupt and SystemExit
except ValueError:           # name the failure being handled
```

Suppress a specific warning by name when you have established it is benign,
and say why in the comment. Blanket suppression hides the next one.

## Paths

No absolute path in analysis code. `here()` in R, `pathlib` in Python, both
anchored on the project root. A hardcoded path is the most common reason a
stage runs for its author and nobody else — and the host and the container
disagree about the root.

Directories are created by the code that writes into them, not assumed.

## Documentation

A comment states intent in one line. A roxygen block or docstring states the
contract: what the function guarantees and what it refuses. This is the same
rule the toolkit applies to its own source; see `AGENTS.md`.

Skip the banner comment. A `# ==== SETUP ====` divider announces a section a
reader can already see, and in a narrative stage the section headings are the
named operations themselves.

## Function length

There is no line count. The test is whether the reader can hold the function
in mind, and whether extracting a piece would give that piece a name worth
having. `stage-narrative` decides what stays visible in a stage;
`helper-family-apis` decides what the extracted piece becomes.
