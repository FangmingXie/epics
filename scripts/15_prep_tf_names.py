#!/usr/bin/env python
"""Stage 4/5 prep — candidate TF list for GRN inference (TF_to_gene).

tf_names = (TFs annotated in the aertslab v10nr motif collection) ∩ (genes detected in the
L2/3 RNA). Only motif-backed TFs can form eRegulons downstream, and a TF must be expressed to
be a candidate regulator, so we intersect the motif-table TFs with genes detected in ≥ MIN_CELLS
L2/3 nuclei. Fail-fast if the list is implausibly small.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/15_prep_tf_names.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# --- paths (capitalized per CLAUDE.md) ---
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/wang25"
MOTIF_TBL = "/data/qlyu/data/common/aertslab_motif_collection/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"
L23_GEX_H5AD = f"{DATA}/v1_l23_gex.h5ad"        # Stage 1 (raw counts in .X, 5,558 cells)
TF_NAMES_OUT = f"{DATA}/tf_names.txt"           # OUTPUT: one TF symbol per line

MIN_CELLS = 10                                  # TF must be detected (count>0) in ≥ this many nuclei
GENE_COL = "feature_name"                       # var column holding gene symbols

# known L2/3 IT markers for a sanity report (Stage 6 positive controls)
KNOWN_L23_TFS = ["CUX2", "SATB2", "POU3F2", "RORB", "MEF2C", "FOXP1"]

for f in (MOTIF_TBL, L23_GEX_H5AD):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")

# --- 1. TF universe from the motif annotation table ---
tbl = pd.read_csv(MOTIF_TBL, sep="\t")
motif_tfs = set(tbl["gene_name"].dropna().unique()) - {"None"}
print(f"[tf_names] motif-annotated TFs (v10nr): {len(motif_tfs):,}")

# --- 2. genes present + detection counts in L2/3 RNA ---
adata = sc.read_h5ad(L23_GEX_H5AD)
if GENE_COL not in adata.var.columns:
    raise KeyError(f"var column '{GENE_COL}' not found; have {list(adata.var.columns)}")
symbols = adata.var[GENE_COL].astype(str).values
X = adata.X
n_cells_detected = np.asarray((X > 0).sum(axis=0)).ravel() if sparse.issparse(X) else (X > 0).sum(axis=0)
detect = pd.Series(n_cells_detected, index=symbols)
# collapse duplicate symbols by max detection
detect = detect.groupby(level=0).max()
print(f"[tf_names] genes in RNA panel: {adata.n_vars:,} ({detect.size:,} unique symbols)")

# --- 3. intersect + expression filter ---
in_panel = sorted(motif_tfs & set(detect.index))
expressed = sorted([tf for tf in in_panel if detect[tf] >= MIN_CELLS])
print(f"[tf_names] motif-TFs in RNA panel: {len(in_panel):,}")
print(f"[tf_names] ...detected in ≥{MIN_CELLS} nuclei (written): {len(expressed):,}")
# transparency: a couple of alternate thresholds
for thr in (1, 50, 200):
    print(f"           (≥{thr:>3} nuclei: {sum(detect[tf] >= thr for tf in in_panel):,})")

if len(expressed) < 500:
    raise ValueError(f"only {len(expressed)} TFs — implausibly small; check inputs/threshold.")

# --- 4. write ---
with open(TF_NAMES_OUT, "w") as fh:
    fh.write("\n".join(expressed) + "\n")
print(f"[tf_names] wrote {TF_NAMES_OUT} ({len(expressed):,} TFs)")

# --- 5. sanity: known L2/3 TFs ---
present = [tf for tf in KNOWN_L23_TFS if tf in expressed]
missing = [tf for tf in KNOWN_L23_TFS if tf not in expressed]
print(f"[tf_names] known L2/3 TFs present: {present}")
if missing:
    print(f"[tf_names] known L2/3 TFs NOT in list: {missing} "
          f"(detection: {{ {', '.join(f'{tf}:{int(detect.get(tf, 0))}' for tf in missing)} }})")
