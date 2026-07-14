# Plan: Install SCENIC+ into the `epics` conda env

## Context

The `epics` project needs SCENIC+ (gene regulatory network inference from scRNA + scATAC
multiome data) available in its dedicated conda env. The env was initially created empty at
**Python 3.14.6**, but SCENIC+'s dependency stack has **no builds for Python 3.14** — it
requires **Python 3.11**. A fully working SCENIC+ install already exists on this machine in the
`scenicplus2` env (Python 3.11.4, `scenicplus` 1.0a2). Rather than risk SCENIC+'s notoriously
fragile dependency resolution, we reproduce that proven env exactly.

**Decision (user-approved):** recreate `epics` as a clone of `scenicplus2`.

## What gets installed

The full `scenicplus2` stack (291 pip pkgs + conda bioconda tools). Core packages
(see `env/scenicplus_core_packages.txt` for the pinned reference list):

| Package | Version | | Package | Version |
|---|---|---|---|---|
| scenicplus | 1.0a2 | | numpy | 1.26.4 |
| pycisTopic | 2.0a0 | | pandas | 1.5.0 |
| pycistarget | 1.1 | | scipy | 1.12.0 |
| ctxcore | 0.2.0 | | numba | 0.59.0 |
| arboreto | 0.1.6 | | scikit-learn | 1.3.2 |
| loomxpy | 0.4.2 | | ray | 2.9.3 |
| scanpy | 1.8.2 | | dask | 2024.2.1 |
| anndata | 0.10.5.post1 | | polars | 0.20.13 |
| mudata | 0.2.3 | | pyranges | 0.0.111 |
| pyBigWig / pybiomart | 0.3.22 / 0.2.0 | | MACS2 | 2.2.9.1 |
| gseapy / tspex | 0.10.8 / 0.6.3 | | bedtools | 2.31.1 (bioconda) |

Python: **3.11.4**.

## Steps

1. **Snapshot current env** (CLAUDE.md rule — already done): current empty-epics list saved to
   `env/epics_before_scenicplus_conda.txt`; SCENIC+ reference list saved to
   `env/scenicplus_core_packages.txt`.
2. **Remove the empty py3.14 epics env**: `conda env remove -n epics -y`.
   (Nothing depends on it — it was just created and is Python-only.)
3. **Clone the proven env**: `conda create -n epics --clone scenicplus2 -y`.
   This copies the exact working package set, including the `scenicplus` install and
   bioconda `bedtools`/`macs2`.
4. **Record the resulting env** into `env/` for reproducibility:
   - `conda list -n epics > env/epics_conda_list.txt`
   - `conda run -n epics pip freeze > env/epics_pip_freeze.txt`

## Verification

- `conda run -n epics python --version` → `Python 3.11.4`.
- `conda run -n epics python -c "import scenicplus, pycistopic, pycistarget, scanpy, anndata; print(scenicplus.__version__)"`
  imports cleanly and prints `1.0a2`.
- `conda run -n epics which bedtools macs2` resolve inside the env.

## Notes / caveats

- This makes `epics` a functional duplicate of `scenicplus2`. That is intentional — `epics` is the
  project-dedicated env going forward; `scenicplus2` remains untouched as the source.
- No package installs beyond the clone. Any later additions follow the CLAUDE.md rule
  (ask first, save `env/` snapshot first, conda before pip).
- `scenicplus` is under an academic non-commercial license (VIB) — already accepted in the
  source env; the clone inherits it.
