# Headless workflow — running review notebooks with no VS Code

The decision-gate notebook is IDE-agnostic on purpose. It needs three things, all of which
work under a bare `docker exec` + tmux + nvim session:

1. an R session (radian) started at the compartment root,
2. a way to send chunk code to that session,
3. a way to see plots.

None of it needs a Jupyter kernel. The knitr engine sends plain R to a plain R session; the
"Jupyter sidecar" experience only appears if a document opts into the jupyter engine or you
open an `.ipynb`. Keep the knitr engine and you stay kernel-free.

## The loop

```bash
docker exec -it [container] bash
cd /workspaces/.../[compartment]         # start AT the compartment root (see root-resolution.md)
tmux new -s nb                            # pane 0: editor, pane 1: R
radian                                    # in a second pane
```

- **Send code to R.** Simplest: `tmux send-keys -t nb:0.1 'CODE' Enter`. Better: an nvim R
  plugin (R.nvim or Nvim-R) that sends the line/paragraph/chunk under the cursor to the radian
  pane — the same "send to terminal" the VS Code R extension does, without VS Code.

- **See plots.** Start an **httpgd** device and open its URL in your host browser:

  ```r
  httpgd::hgd()          # prints a http://127.0.0.1:PORT/... URL
  httpgd::hgd_browse()   # or open it yourself; forward the port if remote
  ```

  Every plot appears live in the browser tab — this is the "modular right-side popup" VS Code
  gives you, decoupled from the editor. `ragg`/cairo (set in the notebook) also mean any
  `png()`/`ggsave()` writes real files headlessly.

- **Produce the committed artifacts.** Same command as under VS Code:

  ```bash
  Rscript 02_analysis/notebooks/render.R 02_analysis/notebooks/[nb]/[nb].qmd both
  ```

  With the Quarto CLI present (scdock-r-dev >= v0.5.6) you can also get live preview:

  ```bash
  quarto preview 02_analysis/notebooks/[nb]/[nb].qmd    # live-reload in the browser
  quarto render  02_analysis/notebooks/[nb]/[nb].qmd --to gfm,html
  ```

## VS Code, for contrast

Under VS Code the same notebook behaves as follows, and it is the correct non-Jupyter path:

- `r.alwaysUseActiveTerminal: false` (shipped in the devcontainer template) → the R extension
  runs a **dedicated** radian terminal and sends chunk code there;
- plots render via **httpgd** into a right-side webview (not a Jupyter kernel);
- Quarto "Run Cell" on an `{r}` chunk routes through that same R terminal.

So the choice between VS Code and bare nvim is a preference, not a constraint: both drive the
same radian session, the same httpgd, the same `render.R`.
