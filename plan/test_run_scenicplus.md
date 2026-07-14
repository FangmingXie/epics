# Plan: Minimal SCENIC+ test run in the `epics` env

## Context

`epics` now contains a full SCENIC+ 1.0a2 install (cloned from `scenicplus2`, Python 3.11.4).
Before using it for real analysis we want a **minimal end-to-end test** that (1) confirms the
environment is functional and (2) exercises the actual GRN-inference pipeline on real data —
without downloading anything or running a multi-hour job.

We exploit an existing, known-good SCENIC+ run already on disk (mouse V1 astrocyte, "t2p1"):

- Prepared multiome MuData: `/home/qlyu/mydata/data/v1_astro/t2p1/spout/ACC_GEX.h5mu` (18 GB)
- Motif enrichment: `/home/qlyu/mydata/data/v1_astro/t2p1/spout/ctx_results.hdf5`
- TF list: `/data/qlyu/code/v1astro/run_scenic/scplus_pipeline/Snakemake/tf_names.txt`
- Reference outputs (for schema/sanity checks): `t2p1/spout/{tf_to_gene_adj.tsv, region_to_gene_adj.tsv, eRegulon_direct.tsv}`

The new SCENIC+ CLI exposes the pipeline as `scenicplus grn_inference {TF_to_gene, region_to_gene, eGRN, ...}`.
The regression steps (`TF_to_gene`, `region_to_gene`) need only the multiome MuData + TF names —
no motif DBs — so subsampling the MuData makes them fast and robust.

**Scope (user-approved):** smoke test + subsampled real run. All I/O under `epics/data/test_run/` (gitignored).

## Layout

```
epics/
├── scripts/
│   ├── test01_smoke.py            # Tier 1: env integrity
│   ├── test02_make_mini_mudata.py # Tier 2a: subsample ACC_GEX.h5mu -> tiny MuData
│   └── test03_run_grn.sh          # Tier 2b-d: run scenicplus grn_inference steps
└── data/test_run/                 # gitignored inputs + outputs
    ├── ACC_GEX.small.h5mu
    ├── tf_names.txt               # copied from the astro run
    ├── tf_to_gene_adj.tsv
    ├── region_to_gene_adj.tsv
    └── eRegulon.tsv
```

## Steps

### Tier 1 — Smoke test (`scripts/test01_smoke.py`)  — must pass, ~1 min
- `scenicplus --help` and import the full stack: `scenicplus, pycisTopic, pycistarget, ctxcore,
  arboreto, scanpy, anndata, mudata, ray, numba`.
- Exercise the compiled/JIT core on synthetic data: run `arboreto.grnboost2` on a tiny random
  expression matrix (~50 cells × 20 genes, 3 TFs) → returns a non-empty adjacency DataFrame.
  This confirms numba/dask/ray + GBM work in `epics`.

### Tier 2a — Build the mini MuData (`scripts/test02_make_mini_mudata.py`) — ~2–5 min
Define all paths (capitalized) at the top per CLAUDE.md.
- Open the 18 GB `ACC_GEX.h5mu` with `mudata.read_h5mu(..., backed=True)` (lazy — no full load);
  print modality names (expected `scRNA`/`scATAC`) and shapes. Fail fast if modalities are absent.
- **Subsample to guarantee signal, not just size**: seed the gene set from TFs in `tf_names.txt`
  plus target genes that appear in the reference `eRegulon_direct.tsv`; keep ~200 genes (incl.
  ~30 TFs), ~3000 regions (those referenced by the kept genes in `region_to_gene_adj.tsv`), and
  ~500 cells. Write `data/test_run/ACC_GEX.small.h5mu`. Copy `tf_names.txt` filtered to kept TFs.

### Tier 2b/2c — Regression steps (`scripts/test03_run_grn.sh`) — must pass, ~2–10 min
```
scenicplus grn_inference TF_to_gene    --multiome_mudata_fname ACC_GEX.small.h5mu \
    --tf_names tf_names.txt --temp_dir <scratch> --out_tf_to_gene_adjacencies tf_to_gene_adj.tsv --n_cpu 8
scenicplus grn_inference region_to_gene --multiome_mudata_fname ACC_GEX.small.h5mu \
    --temp_dir <scratch> --out_region_to_gene_adjacencies region_to_gene_adj.tsv --n_cpu 8
```
(Confirm exact `region_to_gene` flags via `--help` at run time.)

### Tier 2d — eGRN → eRegulons (best-effort, representative final output)
Run `scenicplus grn_inference eGRN` combining the two adjacency tables + the existing
`ctx_results.hdf5`, with **relaxed params** to guarantee non-empty output on the small set
(`--min_target_genes 3`, low `top_n`, `gsea_n_perm 100`). Write `data/test_run/eRegulon.tsv`.
- Fallback if the subsample yields zero eRegulons: rerun `eGRN` on the **full precomputed**
  `t2p1/spout/{tf_to_gene_adj.tsv, region_to_gene_adj.tsv}` + `ctx_results.hdf5` with the same
  relaxed params — still exercises the eGRN code path and produces a real eRegulon table.

## Verification

- Tier 1 prints `SMOKE OK` and a non-empty grnboost2 adjacency shape.
- `data/test_run/tf_to_gene_adj.tsv` and `region_to_gene_adj.tsv` exist, are non-empty, and their
  columns match the reference files' headers (`TF/target/importance`, `region/target/importance`).
- `eRegulon.tsv` exists and its columns match the reference `eRegulon_direct.tsv` schema; row count > 0.
- Whole run completes in minutes on CPU (no GPU); peak RAM well under the machine's 1.1 TB.

## Notes / caveats

- Read-only reuse of `/home/qlyu/mydata/data/v1_astro/**` — the test never writes there.
- `data/test_run/` is gitignored; scripts in `scripts/` are committed.
- If MuData modality keys differ from `scRNA`/`scATAC`, adjust the subsample script accordingly
  (discovered in step 2a, fail-fast otherwise).
- No new packages are installed; if any turn out to be missing, stop and ask (CLAUDE.md rule).
