#!/usr/bin/env python
"""Stage 1 (gao25 · l4) — Select cortical L4/5 IT CTX Glut + multiome barcodes (Allen DevVIS, mouse).

Fork of scripts/10.gao25.select_l23.py for the deeper IT layer L4 (plan/06). Method is unchanged; only
the subclass label, paths, and the fail-fast intersection floor differ. The atlas has no pure "L4 IT"
subclass — L4 lives in the merged `L4/5 IT CTX Glut` subclass (user-confirmed), so the "l4" run selects
that subclass. gao25's RNA `X` is already int32 raw counts and `var_names` are already mouse gene symbols.

Scope: `subclass_label == "L4/5 IT CTX Glut"` AND `roi in VIS_ROIS` (5 cortical/visual ROIs). Final
multiome cells = RNA cellNames ∩ ATAC _index (the common pairing key). We standardize the RNA obs_names
to `cellNames` so cisTopic (named by ATAC `_index`) and this RNA h5ad reconcile after the ___gao25l4
suffix in Stage 4.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/10.gao25.l4.select.py
"""
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# --- paths (capitalized per CLAUDE.md) ---
ORIG = "/home/qlyu/mydata/data/gao25/orig"
RNA_H5AD = f"{ORIG}/DevVIS_multiome_snRNA_processed.h5ad"      # 5.6G, genes×cells int32 raw counts
ATAC_H5AD = f"{ORIG}/DevVIS_multiome_snATAC_processed.h5ad"    # 45G, peaks×cells matrix (obs only read here)
DATA = "/data/qlyu/project/epics/data/gao25/l4"
L4_GEX_H5AD = f"{DATA}/vis_l4_gex.h5ad"                        # output: cortical L4/5, raw counts in .X
L4_BARCODES = f"{DATA}/vis_l4_barcodes.txt"                    # output: final RNA∩ATAC barcodes (common key)

# --- selection criteria ---
SUBCLASS_COL = "subclass_label"
L4_SUBCLASS = "L4/5 IT CTX Glut"     # atlas has no pure "L4 IT"; L4 lives in the merged L4/5 subclass
ROI_COL = "roi"
VIS_ROIS = ["Mouse Multiome VIS", "Mouse Multiome VIS - VISam-VISpm",
            "Mouse Multiome VIS - VISal-VISpl-VISl-VISli-VISpor",
            "Mouse Multiome VISp", "Mouse Multiome Isocortex"]
CELLNAMES_COL = "cellNames"     # RNA obs column = common multiome key (== ATAC _index)
MIN_INTERSECT = 15000           # fail-fast floor (plan: expect ≈ 21.6k)


def is_integer_counts(x) -> bool:
    """True if the matrix holds non-negative integer-valued counts (float dtype allowed)."""
    m = x[:2000] if x.shape[0] > 2000 else x
    data = m.data if sparse.issparse(m) else np.asarray(m).ravel()
    if data.size == 0:
        return False
    return bool((data >= 0).all() and np.allclose(data, np.rint(data)))


def read_atac_l4_barcodes():
    """Read ATAC obs only (h5py) and return the set of _index barcodes for cortical L4/5 cells."""
    def col(f, name):
        g = f["obs"][name]
        if isinstance(g, h5py.Group):                 # categorical
            cats = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
            return np.asarray(cats)[g["codes"][:]]
        arr = g[:]
        return np.asarray([x.decode() if isinstance(x, bytes) else x for x in arr])

    with h5py.File(ATAC_H5AD, "r") as f:
        idx = col(f, "_index")
        sub = col(f, SUBCLASS_COL)
        roi = col(f, ROI_COL)
    mask = (sub == L4_SUBCLASS) & np.isin(roi, VIS_ROIS)
    print(f"  ATAC cortical L4/5 cells: {int(mask.sum()):,} (of {len(idx):,})")
    return set(idx[mask])


def main():
    print(f"[Stage 1 gao25·l4] reading ATAC obs from {ATAC_H5AD}")
    atac_bc = read_atac_l4_barcodes()

    print(f"[Stage 1 gao25·l4] reading {RNA_H5AD}")
    adata = sc.read_h5ad(RNA_H5AD)
    print(f"  full RNA: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in (SUBCLASS_COL, ROI_COL, CELLNAMES_COL):
        if col not in adata.obs.columns:
            raise KeyError(f"obs column '{col}' not found; available: {list(adata.obs.columns)}")

    mask = (adata.obs[SUBCLASS_COL].astype(str) == L4_SUBCLASS) & \
           (adata.obs[ROI_COL].astype(str).isin(VIS_ROIS))
    it = adata[mask.values].copy()
    print(f"  RNA cortical L4/5 cells: {it.n_obs:,}")
    print("  per-ROI counts:\n" + it.obs[ROI_COL].astype(str).value_counts().to_string())

    # guarantee raw integer counts in .X, drop .raw/layers so downstream reads .X unambiguously
    if not is_integer_counts(it.X):
        raise ValueError("RNA .X is not integer counts. SCENIC+ requires raw counts; inspect the h5ad.")
    print("  raw integer counts confirmed in .X")
    it.raw = None
    it.layers.clear()

    # dedup duplicate gene symbols (keep highest-total-count per symbol) so var_names stay unique
    totals = np.asarray(it.X.sum(axis=0)).ravel()
    keep_order = (
        pd.DataFrame({"sym": it.var_names.astype(str).values, "tot": totals,
                      "pos": np.arange(it.n_vars)})
        .sort_values("tot", ascending=False).drop_duplicates("sym", keep="first")
        .sort_values("pos")["pos"].values
    )
    n_before = it.n_vars
    if len(keep_order) < n_before:
        it = it[:, keep_order].copy()
    it.var_names_make_unique()
    print(f"  genes: {n_before:,} -> {it.n_vars:,} unique symbols")

    # standardize obs_names to the common multiome key (RNA cellNames == ATAC _index)
    it.obs_names = it.obs[CELLNAMES_COL].astype(str).values
    if not it.obs_names.is_unique:
        raise ValueError("RNA cellNames are not unique after L4/5 subset — cannot use as multiome key.")

    # intersect RNA cellNames ∩ ATAC _index -> final multiome cells
    rna_bc = set(it.obs_names)
    inter = rna_bc & atac_bc
    print(f"  RNA-only L4/5: {len(rna_bc):,}  ATAC-only L4/5: {len(atac_bc):,}  intersection: {len(inter):,}")
    if len(inter) <= MIN_INTERSECT:
        raise ValueError(f"intersection collapsed to {len(inter)} (<= {MIN_INTERSECT}); investigate barcode keys.")

    it = it[it.obs_names.isin(inter)].copy()
    it = it[np.argsort(it.obs_names.values)].copy()   # deterministic order
    print(f"  final multiome RNA subset: {it.n_obs:,} cells x {it.n_vars:,} genes")

    it.write_h5ad(L4_GEX_H5AD)
    print(f"  wrote {L4_GEX_H5AD}")

    pd.Series(it.obs_names).to_csv(L4_BARCODES, index=False, header=False)
    print(f"  wrote {L4_BARCODES} ({it.n_obs:,} barcodes)")
    print(f"  first barcode: {it.obs_names[0]}")


if __name__ == "__main__":
    main()
