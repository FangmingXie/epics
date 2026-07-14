#!/usr/bin/env python
"""Stage 3c — pycisTopic object + LDA topic modelling for V1 L2/3 ATAC.

Mirrors v1astro 11.run_pycistopic.py, adapted to: the UCSC/primary-chrom fragments (Stage 2b),
the V1 L2/3 barcodes (Stage 1), the consensus peaks (Stage 3a). Builds the cisTopic object
restricted to L2/3 barcodes, sweeps LDA topics, auto-selects the best model, binarizes topics
into per-topic region sets (BED) for motif enrichment (Stage 4).

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/14_run_pycistopic.py

NOTE: only one cell type is present (EN-L2_3-IT), so the DAR (find_diff_features) step from the
v1astro template is intentionally omitted — region sets come from topic binarization only.
"""
import os
import pickle

# ─────────────────────────────────────────────
# Paths (capitalized per CLAUDE.md) + params
# ─────────────────────────────────────────────
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/wang25"
WORK = f"{DATA}/work"
DB = f"{DATA}/db"

FRAG_UCSC_BGZ = f"{DATA}/atac_fragments.ucsc.tsv.bgz"          # Stage 2b (UCSC, primary chroms)
L23_BARCODES = f"{DATA}/v1_l23_barcodes.txt"                   # Stage 1 (5,558 barcodes)
PEAKS_FILTERED = f"{WORK}/l23_summits_501_filtered.bed"        # Stage 3a consensus peaks
HG38_BLACKLIST = f"{DB}/hg38-blacklist.v2.bed"                 # Stage 2

CISTOPIC_OBJ = f"{WORK}/cistopic_obj_model.pkl"                # output: cisTopic obj + LDA model
MODELS_DIR = f"{WORK}/lda_models"
REGION_SETS_DIR = f"{WORK}/region_sets"
TOPICS_SUBFOLDER = "topics_top_3k"

PROJECT = "l23"
N_TOPICS = [10, 20, 30, 40]      # LDA sweep; auto-select by coherence
N_ITER = 150
N_CPU = 8
RANDOM_STATE = 555
NTOP = 3000                      # top regions per binarized topic

os.makedirs(WORK, exist_ok=True)
os.makedirs(MODELS_DIR, exist_ok=True)

# fail-fast on missing inputs
for f in (FRAG_UCSC_BGZ, L23_BARCODES, PEAKS_FILTERED, HG38_BLACKLIST):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics

# ─────────────────────────────────────────────
# 1. Build the cisTopic object (restricted to L2/3 barcodes)
# ─────────────────────────────────────────────
valid_bc = open(L23_BARCODES).read().split()
print(f"[Stage 3c] building cisTopic object on {len(valid_bc):,} L2/3 barcodes")
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=FRAG_UCSC_BGZ,
    path_to_regions=PEAKS_FILTERED,
    path_to_blacklist=HG38_BLACKLIST,
    valid_bc=valid_bc,           # restrict to L2/3 cells
    n_cpu=N_CPU,
    project=PROJECT,
    # do NOT pass split_pattern — it creates duplicate cells (v1astro note)
)
print(f"  ATAC binary matrix: {cistopic_obj.binary_matrix.shape}")

# ─────────────────────────────────────────────
# 2. LDA topic modelling (sweep + auto-select)
# ─────────────────────────────────────────────
print(f"[Stage 3c] LDA sweep n_topics={N_TOPICS}, n_iter={N_ITER}")
models = run_cgs_models(
    cistopic_obj,
    n_topics=N_TOPICS,
    n_cpu=N_CPU,
    n_iter=N_ITER,
    random_state=RANDOM_STATE,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=MODELS_DIR,
)
model = evaluate_models(
    models,
    select_model=None,           # None → auto-select best coherence
    return_model=True,
    metrics=["Minmo_2011", "loglikelihood"],
    plot=False,
    save=os.path.join(WORK, "model_evaluation.pdf"),
)
cistopic_obj.add_LDA_model(model)
print(f"  selected model: {model.n_topic} topics")

with open(CISTOPIC_OBJ, "wb") as f:
    pickle.dump(cistopic_obj, f)
print(f"  wrote {CISTOPIC_OBJ}")

# ─────────────────────────────────────────────
# 3. Binarize topics → per-topic region sets (BED)
# ─────────────────────────────────────────────
print(f"[Stage 3c] binarizing topics (ntop={NTOP})")
region_bin = binarize_topics(cistopic_obj, method="ntop", ntop=NTOP, plot=False)

out_dir = os.path.join(REGION_SETS_DIR, TOPICS_SUBFOLDER)
os.makedirs(out_dir, exist_ok=True)
for topic, regions in region_bin.items():
    names = regions.index.values            # e.g. "chr1:100-200"
    with open(os.path.join(out_dir, f"{topic}.bed"), "w") as fh:
        for r in names:
            chrom, rest = r.split(":")
            start, end = rest.split("-")
            fh.write(f"{chrom}\t{start}\t{end}\n")
print(f"  wrote {len(region_bin)} topic BEDs -> {out_dir}")
print("[Stage 3c] complete.")