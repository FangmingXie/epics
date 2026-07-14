"""Tier 2a: build a tiny multiome MuData for a minimal SCENIC+ GRN test.

Subsamples the existing prepared `ACC_GEX.h5mu` (mouse V1 astrocyte, t2p1) down to a
small set of cells / genes / regions, seeded from the genes, TFs and regions that appear
in the reference eRegulon table so the downstream GRN steps produce non-empty output.

Run: conda run --no-capture-output -n epics python -u scripts/test02_make_mini_mudata.py
"""

import os
import numpy as np
import pandas as pd
import mudata as mu

# ---------------------------------------------------------------- paths / params
MU_IN          = "/home/qlyu/mydata/data/v1_astro/t2p1/spout/ACC_GEX.h5mu"
EREGULON_REF   = "/home/qlyu/mydata/data/v1_astro/t2p1/spout/eRegulon_direct.tsv"
TF_NAMES_IN    = "/data/qlyu/code/v1astro/run_scenic/scplus_pipeline/Snakemake/tf_names.txt"

OUT_DIR        = "/data/qlyu/project/epics/data/test_run"
MU_OUT         = os.path.join(OUT_DIR, "ACC_GEX.small.h5mu")
TF_NAMES_OUT   = os.path.join(OUT_DIR, "tf_names.txt")

N_CELLS        = 600
MAX_GENES      = 250
MAX_REGIONS    = 4000
SEED           = 666

RNG = np.random.default_rng(SEED)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    # --- seed sets from the reference eRegulons (genes/TFs/regions with real signal) ---
    ereg = pd.read_csv(EREGULON_REF, sep="\t", usecols=["Region", "Gene", "TF"])
    seed_genes = set(ereg["Gene"]) | set(ereg["TF"])
    seed_regions = list(dict.fromkeys(ereg["Region"]))  # unique, preserve order
    all_tfs = set(pd.read_csv(TF_NAMES_IN, header=None)[0])
    print(f"reference eRegulons: {len(seed_genes)} genes/TFs, {len(seed_regions)} regions")

    # --- load prepared MuData (dense scATAC ~18GB; machine has ample RAM) ---
    print("loading MuData (this reads the full object) ...", flush=True)
    mdata = mu.read_h5mu(MU_IN)
    rna, atac = mdata["scRNA"], mdata["scATAC"]
    print(f"loaded scRNA {rna.shape}, scATAC {atac.shape}", flush=True)

    # --- pick genes: all seed genes present, padded with random genes up to MAX_GENES ---
    rna_vars = set(rna.var_names)
    genes = [g for g in seed_genes if g in rna_vars]
    if len(genes) < MAX_GENES:
        pad = [g for g in rna.var_names if g not in set(genes)]
        genes += list(RNG.choice(pad, size=min(MAX_GENES - len(genes), len(pad)), replace=False))
    genes = list(dict.fromkeys(genes))[:MAX_GENES]

    # --- pick regions: seed regions present, capped ---
    atac_vars = set(atac.var_names)
    regions = [r for r in seed_regions if r in atac_vars][:MAX_REGIONS]

    # --- pick cells ---
    n = mdata.n_obs
    cell_idx = np.sort(RNG.choice(n, size=min(N_CELLS, n), replace=False))

    # fail-fast on empty selections
    assert genes and regions, f"empty selection: {len(genes)} genes, {len(regions)} regions"

    rna_small = rna[cell_idx, genes].copy()
    atac_small = atac[cell_idx, regions].copy()
    small = mu.MuData({"scRNA": rna_small, "scATAC": atac_small})
    small.write(MU_OUT)

    # --- TF names present in the mini scRNA (for TF_to_gene) ---
    tfs_present = [g for g in genes if g in all_tfs]
    pd.Series(tfs_present).to_csv(TF_NAMES_OUT, index=False, header=False)

    print(f"wrote {MU_OUT}")
    print(f"  scRNA  {rna_small.shape}  ({len(tfs_present)} of them are TFs)")
    print(f"  scATAC {atac_small.shape}")
    print(f"wrote {TF_NAMES_OUT} ({len(tfs_present)} TFs)")
    print("MINI MUDATA OK")


if __name__ == "__main__":
    main()
