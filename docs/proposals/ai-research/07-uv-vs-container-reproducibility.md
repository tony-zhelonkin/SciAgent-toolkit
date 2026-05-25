# 07 — uv vs container reproducibility

ADR-009's "Personal take" section poses the question directly: working in Docker containers with pre-installed packages, does `uv` add reproducibility on top?

## Where uv and containers overlap

Both pin software state. Both reduce "works on my machine" failures. Both can be lockfile-driven. In a CI context both can produce byte-identical environments across runs.

The overlap is real but the **layer** differs:

- Container pins: OS, system libraries, Python interpreter, baseline Python packages installed at image-build time.
- uv pins: exact Python package versions per-project at resolve time; reads `pyproject.toml` + `uv.lock`; deterministic install. Resolves 10–100× faster than pip.

A container with everything pre-installed pins the *snapshot*; uv pins the *recipe*. They're not redundant — they're different layers of the same reproducibility stack.

## Where they complement

Concrete cases where uv adds something the container does not:

1. **Per-project drift.** Two repos sharing a container will each have their own dependency requirements that may conflict. uv puts each in its own `.venv` deterministically; the container provides the interpreter and base libraries.
2. **Lockfile travels with the code.** `uv.lock` checked into git means anyone (or any CI) reproduces the exact resolved tree without needing the container. The container is a fast path; the lockfile is the truth.
3. **LLM-generated code with novel imports.** A tree-search node (ADR-004) generates `import some_obscure_package`. The container probably doesn't have it. uv resolves and installs in seconds; pip would take minutes. This case is real for ERA-style search.
4. **Distroless container patterns.** Modern best practice (2026) is to build minimal containers using `uv` at image-build time so the image contains only the runtime dependencies — no apt, no /bin/sh, no pip. Containers and uv are *complementary*, not alternatives. See [Distroless Python Containers with uv — Nerd Level Tech, 2026](https://nerdleveltech.com/distroless-python-containers-with-uv-tutorial).

## The specific question: does uv matter inside a fully-pre-installed dev container?

For **interactive analysis work** in your typical Docker container — where 95% of your dependencies are already installed and you're not generating novel code with unusual imports — uv adds:

- **A lockfile.** Your container is the dev environment; collaborators not using your container can reproduce via `uv sync`. The container becomes the convenience layer, not the source of truth.
- **Per-project isolation.** If you maintain multiple analysis projects in the same container, each can have its own venv with its own pinned versions. The container's site-packages becomes a fallback / system pool, not the source of project dependencies.
- **Faster ad-hoc installs.** When you need a one-off package, `uv add` is ~30× faster than `pip install`. Minor but compounds.

For **LLM-generated code in a tree-search sandbox** (ADR-004), uv genuinely matters:

- Tree-search nodes will generate code that imports packages not in the container's site-packages. Pre-installation cannot cover this; the search is exploring novel solutions.
- Per-node venv isolation prevents one bad node from polluting subsequent nodes (e.g., a node that pins `numpy==1.20` shouldn't leak into a node expecting `numpy>=2.0`).
- Venv reuse across sibling nodes with identical dependency sets is a real performance optimisation (per skill 2.3 body).

This is the case where uv is **not optional**. The container's pre-installed packages don't help when the search is generating unforeseen imports.

## Where the spec's claim is exposed

ADR-009 says: *"The Antigravity Skills paper reports significant reproducibility loss when sandboxed code runs in user-variable Python environments. They adopted uv to pin environments per skill execution."* Per file 06, this attribution is unverified in publicly accessible ERA material. The motivation is sound regardless of whether the paper actually says this — but the spec should not lean on an unconfirmed citation.

## Recommendation

**Make uv mandatory for the `scorable` overlay (skill 2.3 `sandbox-execution`); document a "containerized user" fast-path for the base role.**

Concretely:

- In `base` role (interactive analysis): uv is not required. If your container has everything pre-installed, no friction. `sciagent doctor` can soft-warn if `uv` is missing but should not fail.
- In `scorable` overlay (tree-search node execution): uv is required. `sandbox-execution` creates a `uv venv` per node, installs the resolved dep set, runs, tears down or caches. Without uv, the overlay does not function.
- Provide a one-line installer in `sciagent doctor --fix` that runs `pipx install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh` — eliminates the friction objection for users who don't already have it.

This preserves ADR-009's intent for the place it matters (sandboxed LLM-generated code) without imposing a system dependency on users who only want the base interactive role.

## Counter-position to consider

The minority view: **make uv mandatory always.** Argument: bifurcating the dependency story creates two code paths to maintain. The cost of `pipx install uv` is genuinely trivial (~5 seconds, one command). Forcing it removes a configuration matrix.

Worth grilling. Reasonable people disagree on this — it's a values call about onboarding friction vs. design uniformity.

Sources:
- [How to use a uv lockfile for reproducible Python environments — pydevtools](https://pydevtools.com/handbook/how-to/how-to-use-a-uv-lockfile-for-reproducible-python-environments/)
- [Locking environments — uv docs](https://docs.astral.sh/uv/pip/compile/)
- [Distroless Python Containers with uv: 2026 Tutorial — Nerd Level Tech](https://nerdleveltech.com/distroless-python-containers-with-uv-tutorial)
- [Best Python Package Managers 2026 — Scopir](https://scopir.com/posts/best-python-package-managers-2026/)
- [uv on PyPI](https://pypi.org/project/uv/)
