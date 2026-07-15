#!/usr/bin/env python
"""Stage 1 — Select V1 L2/3 (EN-L2_3-IT) neurons from the Wang25 RNA h5ad.

Subset to `Type == "EN-L2_3-IT"` AND `Region == "V1"` (expected 5,558 cells),
guarantee raw integer counts live in `.X`, and emit a barcode list for pycisTopic.
Fail-fast: stop with a clear error if the cell count is wrong or counts are not integers.
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# --- paths (capitalized per CLAUDE.md) ---
DATA = "/data/qlyu/project/epics/data/wang25"
GEX_H5AD = f"{DATA}/GEX_with_raw_counts.h5ad"          # input: full RNA (232,328 x 35,477)
L23_GEX_H5AD = f"{DATA}/v1_l23_gex.h5ad"               # output: V1 L2/3 subset, raw counts in .X
L23_BARCODES = f"{DATA}/v1_l23_barcodes.txt"          # output: one barcode per line

# --- selection criteria ---
TYPE_COL, TYPE_VAL = "Type", "EN-L2_3-IT"
REGION_COL, REGION_VAL = "Region", "V1"
EXPECTED_N = 5558
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
    print(f"[Stage 1] reading {GEX_H5AD}")
    adata = sc.read_h5ad(GEX_H5AD)
    print(f"  full RNA: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in (TYPE_COL, REGION_COL):
        if col not in adata.obs.columns:
            raise KeyError(f"obs column '{col}' not found; available: {list(adata.obs.columns)}")

    mask = (adata.obs[TYPE_COL] == TYPE_VAL) & (adata.obs[REGION_COL] == REGION_VAL)
    l23 = adata[mask].copy()
    print(f"  selected {TYPE_COL}=={TYPE_VAL} & {REGION_COL}=={REGION_VAL}: {l23.n_obs:,} cells")
    if l23.n_obs != EXPECTED_N:
        raise ValueError(f"expected {EXPECTED_N} cells, got {l23.n_obs}. Check labels/region values.")

    # guarantee raw integer counts in .X, drop .raw/layers so downstream reads .X unambiguously
    counts, source = pick_counts_matrix(l23)
    print(f"  raw integer counts found in {source}; moving to .X")
    l23.X = counts.copy() if sparse.issparse(counts) else np.asarray(counts).copy()
    l23.raw = None
    l23.layers.clear()

    # Use gene SYMBOLS as var_names — SCENIC+ matches RNA genes to the (symbol-based) gene annotation
    # and tf_names. CELLxGENE ships Ensembl IDs as the index. Collapse duplicate symbols by keeping the
    # highest-total-count gene per symbol so names stay unique AND match the annotation (no -1/-2 suffixes).
    if GENE_SYMBOL_COL not in l23.var.columns:
        raise KeyError(f"var column '{GENE_SYMBOL_COL}' not found; have {list(l23.var.columns)}")
    totals = np.asarray(l23.X.sum(axis=0)).ravel()
    keep_order = (
        pd.DataFrame({"sym": l23.var[GENE_SYMBOL_COL].astype(str).values, "tot": totals,
                      "pos": np.arange(l23.n_vars)})
        .sort_values("tot", ascending=False).drop_duplicates("sym", keep="first")
        .sort_values("pos")["pos"].values
    )
    n_before = l23.n_vars
    l23 = l23[:, keep_order].copy()
    l23.var_names = l23.var[GENE_SYMBOL_COL].astype(str).values
    l23.var_names_make_unique()  # no-op after dedup; guards against any residual collision
    print(f"  var_names -> gene symbols: {n_before:,} genes -> {l23.n_vars:,} unique symbols")

    l23.write_h5ad(L23_GEX_H5AD)
    print(f"  wrote {L23_GEX_H5AD}")

    pd.Series(l23.obs_names).to_csv(L23_BARCODES, index=False, header=False)
    print(f"  wrote {L23_BARCODES} ({l23.n_obs:,} barcodes)")
    print(f"  first barcode: {l23.obs_names[0]}")


if __name__ == "__main__":
    main()
