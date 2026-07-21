#!/usr/bin/env python
"""Stage 3c (gao25 · l6) — pycisTopic object FROM the precomputed count matrix + LDA topic modelling.

Fork of scripts/14.gao25.run_pycistopic.py for the L6 layer (plan/06). Method unchanged; only the
npz/gex paths and the PROJECT tag (→ ___gao25l6 barcode suffix) differ. Build the cisTopic object with
`create_cistopic_object` (in-memory, from the regions×cells matrix produced in Stage 3a). No fragments,
no fragment-based QC (TSS/FRiP); we rely on the Stage-3a peak filter + the h5ad obs QC.

Steps: load the npz matrix -> create_cistopic_object -> attach cell metadata -> LDA sweep + select ->
binarize topics (ntop=3000) -> per-topic region-set BEDs for Stage 4 motif enrichment.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/14.gao25.l6.run_pycistopic.py
"""
import os
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ─────────────────────────────────────────────
# Paths (capitalized per CLAUDE.md) + params
# ─────────────────────────────────────────────
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/gao25/l6"
WORK = f"{DATA}/work"
DB = f"{DATA}/db"

L6_MATRIX_NPZ = f"{WORK}/vis_l6_matrix.npz"                   # Stage 3a: regions×cells CSR + names
L6_GEX_H5AD = f"{DATA}/vis_l6_gex.h5ad"                       # Stage 1: obs carries cell metadata
MM10_BLACKLIST = "/data/qlyu/data/common/mm10/mm10.blacklist.bed"

CISTOPIC_OBJ = f"{WORK}/cistopic_obj_model.pkl"              # output: cisTopic obj + LDA model
MODELS_DIR = f"{WORK}/lda_models"
REGION_SETS_DIR = f"{WORK}/region_sets"
TOPICS_SUBFOLDER = "topics_top_3k"

PROJECT = "gao25l6"                                           # → barcode suffix ___gao25l6 (Stage 4/6)
SPLIT_PATTERN = "___"                                         # → tagged cell names 'barcode___gao25l6'
N_TOPICS = [10, 20, 30, 40]
N_ITER = 150
N_CPU = 3                                                     # LDA: one model per worker across the sweep
RANDOM_STATE = 555
NTOP = 3000                                                  # top regions per binarized topic
META_COLS = ["subclass_label", "roi", "age_label", "sex", "donor_name"]

os.makedirs(WORK, exist_ok=True)
os.makedirs(MODELS_DIR, exist_ok=True)
for f in (L6_MATRIX_NPZ, L6_GEX_H5AD, MM10_BLACKLIST):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")

from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics

# ─────────────────────────────────────────────
# 1. Build the cisTopic object FROM the matrix
# ─────────────────────────────────────────────
npz = np.load(L6_MATRIX_NPZ, allow_pickle=True)
mat = sp.csr_matrix((npz["data"], npz["indices"], npz["indptr"]), shape=tuple(npz["shape"]))
# counts are integer-valued but stored as float64 in the source h5ad; pycisTopic's LDA (collapsed
# Gibbs sampling) requires an INTEGER matrix — cast losslessly to int32.
if not np.issubdtype(mat.dtype, np.integer):
    if not np.allclose(mat.data, np.rint(mat.data)):
        raise ValueError("matrix has non-integer values — cannot be ATAC counts.")
    mat = mat.astype(np.int32)
cell_names = [str(c) for c in npz["cell_names"]]              # ATAC _index (== common key)
region_names = [str(r) for r in npz["region_names"]]         # 'chr:start-end'
print(f"[Stage 3c gao25·l6] matrix regions×cells = {mat.shape[0]:,} x {mat.shape[1]:,}; nnz={mat.nnz:,}")

cistopic_obj = create_cistopic_object(
    fragment_matrix=mat,                 # regions×cells CSR
    cell_names=cell_names,
    region_names=region_names,
    path_to_blacklist=MM10_BLACKLIST,    # drops blacklisted regions
    project=PROJECT,                     # tag_cells → cell names get '___gao25l6'
    split_pattern=SPLIT_PATTERN,
    min_frag=1, min_cell=1,              # peak-level filtering already done in Stage 3a
)
print(f"  cisTopic binary matrix: {cistopic_obj.binary_matrix.shape} (regions × cells)")

# ─────────────────────────────────────────────
# 2. Attach cell metadata (indexed by the TAGGED cell names = barcode + '___gao25l6')
# ─────────────────────────────────────────────
import scanpy as sc
obs = sc.read_h5ad(L6_GEX_H5AD, backed="r").obs
meta = obs.loc[:, META_COLS].copy()
meta.index = [f"{bc}{SPLIT_PATTERN}{PROJECT}" for bc in meta.index]   # match cistopic cell names
meta = meta.reindex(cistopic_obj.cell_names)
if meta.isna().all(axis=0).any():
    raise ValueError(f"cell metadata has all-NaN column(s): {meta.columns[meta.isna().all()].tolist()}")
n_bad = int(meta.isna().any(axis=1).sum())
if n_bad:
    raise ValueError(f"{n_bad} cistopic cells not matched to obs metadata — barcode/tag mismatch")
cistopic_obj.add_cell_data(meta)
print(f"  attached cell metadata {META_COLS}; cell_data shape={cistopic_obj.cell_data.shape}")

# ─────────────────────────────────────────────
# 3. LDA topic modelling (sweep + auto-select)
# ─────────────────────────────────────────────
print(f"[Stage 3c gao25·l6] LDA sweep n_topics={N_TOPICS}, n_iter={N_ITER}")
models = run_cgs_models(
    cistopic_obj,
    n_topics=N_TOPICS, n_cpu=N_CPU, n_iter=N_ITER,
    random_state=RANDOM_STATE, alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False,
    save_path=MODELS_DIR,
)
model = evaluate_models(models, select_model=None, return_model=True,
                        metrics=["Minmo_2011", "loglikelihood"], plot=False,
                        save=os.path.join(WORK, "model_evaluation.pdf"))
cistopic_obj.add_LDA_model(model)
print(f"  selected model: {model.n_topic} topics")

with open(CISTOPIC_OBJ, "wb") as f:
    pickle.dump(cistopic_obj, f)
print(f"  wrote {CISTOPIC_OBJ}")

# ─────────────────────────────────────────────
# 4. Binarize topics → per-topic region-set BEDs (topics only — no DAR, per plan)
# ─────────────────────────────────────────────
print(f"[Stage 3c gao25·l6] binarizing topics (ntop={NTOP})")
region_bin = binarize_topics(cistopic_obj, method="ntop", ntop=NTOP, plot=False)
out_dir = os.path.join(REGION_SETS_DIR, TOPICS_SUBFOLDER)
os.makedirs(out_dir, exist_ok=True)
for topic, regions in region_bin.items():
    with open(os.path.join(out_dir, f"{topic}.bed"), "w") as fh:
        for r in regions.index.values:                   # 'chr:start-end'
            chrom, rest = r.split(":")
            start, end = rest.split("-")
            fh.write(f"{chrom}\t{start}\t{end}\n")
print(f"  wrote {len(region_bin)} topic BEDs -> {out_dir}")
print("[Stage 3c gao25·l6] complete.")
