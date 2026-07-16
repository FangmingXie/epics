#!/usr/bin/env python
"""Stage 1 (v2) — Select ALL IT excitatory neurons (L2/3–L6) from the Wang25 RNA h5ad.

New copy of 10_select_l23.py adapted to the broadened scope (see plan/04):
subset to `Type ∈ {EN-L2_3-IT, EN-L4-IT, EN-L5-IT, EN-L6-IT}` AND `Region ∈ {PFC, V1}`
(expected 49,283 cells), guarantee raw integer counts live in `.X`, and emit a barcode list.
`Type`, `Region`, `Group` are retained in obs — Stage 6 QC reads all three (keyed by barcode).
Fail-fast: stop with a clear error if the cell count is wrong or counts are not integers.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/10.v2.select_it.py
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# --- paths (capitalized per CLAUDE.md) ---
DATA = "/data/qlyu/project/epics/data/wang25"
GEX_H5AD = f"{DATA}/GEX_with_raw_counts.h5ad"          # input: full RNA (232,328 x 35,477)
IT_GEX_H5AD = f"{DATA}/it_gex.h5ad"                    # output: all-IT subset, raw counts in .X
IT_BARCODES = f"{DATA}/it_barcodes.txt"               # output: one barcode per line

# --- selection criteria ---
TYPE_COL = "Type"
IT_TYPES = ["EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L6-IT"]   # layer-resolved IT only
REGION_COL = "Region"
REGIONS = ["PFC", "V1"]                                # drop the tiny "General" pool
EXPECTED_N = 49283
GENE_SYMBOL_COL = "feature_name"    # var column with gene symbols (CELLxGENE index is Ensembl)


def is_integer_counts(x) -> bool:
    """True if the matrix holds non-negative integer-valued counts (float dtype allowed)."""
    m = x[:2000] if x.shape[0] > 2000 else x            # sample rows; enough to catch normalized data
    data = m.data if sparse.issparse(m) else np.asarray(m).ravel()
    if data.size == 0:
        return False
    return bool((data >= 0).all() and np.allclose(data, np.rint(data)))


def pick_counts_matrix(adata):
    """Return (matrix, source_label) for the slot holding raw integer counts. Fail-fast if none."""
    candidates = [("X", adata.X)]
    if "counts" in adata.layers:
        candidates.append(('layers["counts"]', adata.layers["counts"]))
    if adata.raw is not None:
        candidates.append((".raw.X", adata.raw.X))
    for label, mat in candidates:
        if is_integer_counts(mat):
            return mat, label
    checked = ", ".join(l for l, _ in candidates)
    raise ValueError(
        f"No integer-count slot found among [{checked}]. SCENIC+ requires raw counts. "
        "Inspect the h5ad and update this script."
    )


def main():
    print(f"[Stage 1 v2] reading {GEX_H5AD}")
    adata = sc.read_h5ad(GEX_H5AD)
    print(f"  full RNA: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in (TYPE_COL, REGION_COL):
        if col not in adata.obs.columns:
            raise KeyError(f"obs column '{col}' not found; available: {list(adata.obs.columns)}")

    mask = adata.obs[TYPE_COL].isin(IT_TYPES) & adata.obs[REGION_COL].isin(REGIONS)
    it = adata[mask].copy()
    print(f"  selected {TYPE_COL} in {IT_TYPES} & {REGION_COL} in {REGIONS}: {it.n_obs:,} cells")
    print("  per-type counts:\n" + it.obs[TYPE_COL].value_counts().to_string())
    print("  per-region counts:\n" + it.obs[REGION_COL].value_counts().to_string())
    if it.n_obs != EXPECTED_N:
        raise ValueError(f"expected {EXPECTED_N} cells, got {it.n_obs}. Check labels/region values.")

    # guarantee raw integer counts in .X, drop .raw/layers so downstream reads .X unambiguously
    counts, source = pick_counts_matrix(it)
    print(f"  raw integer counts found in {source}; moving to .X")
    it.X = counts.copy() if sparse.issparse(counts) else np.asarray(counts).copy()
    it.raw = None
    it.layers.clear()

    # Use gene SYMBOLS as var_names — SCENIC+ matches RNA genes to the (symbol-based) gene annotation
    # and tf_names. CELLxGENE ships Ensembl IDs as the index. Collapse duplicate symbols by keeping the
    # highest-total-count gene per symbol so names stay unique AND match the annotation (no -1/-2 suffixes).
    if GENE_SYMBOL_COL not in it.var.columns:
        raise KeyError(f"var column '{GENE_SYMBOL_COL}' not found; have {list(it.var.columns)}")
    totals = np.asarray(it.X.sum(axis=0)).ravel()
    keep_order = (
        pd.DataFrame({"sym": it.var[GENE_SYMBOL_COL].astype(str).values, "tot": totals,
                      "pos": np.arange(it.n_vars)})
        .sort_values("tot", ascending=False).drop_duplicates("sym", keep="first")
        .sort_values("pos")["pos"].values
    )
    n_before = it.n_vars
    it = it[:, keep_order].copy()
    it.var_names = it.var[GENE_SYMBOL_COL].astype(str).values
    it.var_names_make_unique()  # no-op after dedup; guards against any residual collision
    print(f"  var_names -> gene symbols: {n_before:,} genes -> {it.n_vars:,} unique symbols")

    it.write_h5ad(IT_GEX_H5AD)
    print(f"  wrote {IT_GEX_H5AD}")

    pd.Series(it.obs_names).to_csv(IT_BARCODES, index=False, header=False)
    print(f"  wrote {IT_BARCODES} ({it.n_obs:,} barcodes)")
    print(f"  first barcode: {it.obs_names[0]}")


if __name__ == "__main__":
    main()
