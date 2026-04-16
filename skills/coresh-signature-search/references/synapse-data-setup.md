# Synapse data setup -- one-time download of preprocessed chunks

> **Claude cannot perform this download.** It requires the user's personal Synapse access token. Ask the user to run the steps below once, then point the skill at the resulting `preprocessed_chunks/` directory.

## What you're downloading

Synapse project **syn66227307** (<https://www.synapse.org/coresh>) hosts the preprocessed CORESH compendium:

- `preprocessed_chunks/hsa/` -- 44,253 human datasets, split into ~100 `.qs2` chunks of ~500 objects each
- `preprocessed_chunks/mmu/` -- 42,224 mouse datasets, same layout

Total size: ~20 GB. Each chunk is self-contained and can be loaded independently with `qs2::qs_read()`.

## Steps (user runs these once)

### 1. Register at Synapse

Create a free account at <https://accounts.synapse.org/authenticated/signTerms>. You must accept the data-use terms before any download.

### 2. Create a personal access token

Go to <https://accounts.synapse.org/authenticated/personalaccesstokens> and create a token with at least `view` and `download` permissions. Save the token string -- it is shown only once.

### 3. Install the CLI

```bash
pip install --upgrade synapseclient
```

### 4. Configure credentials

```bash
synapse config
# Follow the prompts; paste your personal access token when asked.
```

Credentials land in `~/.synapseConfig`.

### 5. Download the compendium

```bash
cd /workspaces/DC_hum_verse/00_data/external     # or wherever you keep large reference data
synapse get -r syn66227307                       # -r = recursive, restore directory tree
```

This creates `preprocessed_chunks/` with `hsa/` and `mmu/` subdirectories. Expect 30-90 min depending on bandwidth.

### 6. Verify presence

```bash
ls preprocessed_chunks/hsa/*.qs2 | wc -l    # expect ~90 chunks
ls preprocessed_chunks/mmu/*.qs2 | wc -l    # expect ~85 chunks
```

Run `scripts/validate_coresh_install.R` to verify from R (package loads, chunks readable, object structure matches expected).

## Path convention

Throughout this skill, `preprocessed_chunks/` is treated as a relative path from the project root. If you put the data elsewhere, set:

```r
Sys.setenv(CORESH_CHUNKS = "/absolute/path/to/preprocessed_chunks")
```

and `coresh_batch.R` will pick it up via `Sys.getenv("CORESH_CHUNKS")`.

## Disk considerations

- 20 GB is substantial. Keep it on fast storage (SSD) -- chunk loading is I/O-bound at full compendium scale.
- Chunks can be streamed from object storage, but the `bplapply` parallel pattern assumes local filesystem access; S3/GCS with FUSE works but halves throughput.
- The compendium is versioned on Synapse. Pin your download to a specific Synapse version if reproducibility matters; note the version in your methods.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Unauthorized` on `synapse get` | Not logged in or token expired | Re-run `synapse config` with a fresh token |
| Download hangs partway | Large file + flaky connection | `synapse get -r syn66227307` is resumable; just re-run |
| `qs_read` fails with "magic number mismatch" | Truncated chunk | Re-download that chunk; the rest remain valid |
| One species missing | Partial download | Re-run `synapse get -r` -- it skips existing files |
