#!/usr/bin/env python
"""Stage 3c (v2) — pycisTopic object + LDA topic modelling for the all-IT ATAC.

New copy of 14_run_pycistopic.py for the all-IT scope (plan/04). Structurally identical: build the
cisTopic object restricted to the IT barcodes, sweep LDA topics, auto-select, binarize topics into
per-topic region sets (BED) for motif enrichment (Stage 4).

Per the approved plan, DAR region sets stay OMITTED (topics-only region sets — user decision), even
though multiple cell types are now present. Region sets come from topic binarization only.

Changes vs the L2/3 run: paths → work_it/, PROJECT="it" (barcode suffix ___it), widened topic sweep
N_TOPICS=[20,30,40,50,60] to accommodate the far richer structure of 4 layers × 2 regions × 5 stages.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/14.v2.run_pycistopic_it.py
"""
import os
import pickle

# ─────────────────────────────────────────────
# Paths (capitalized per CLAUDE.md) + params
# ─────────────────────────────────────────────
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/wang25"
WORK = f"{DATA}/work_it"                                       # NEW isolated work dir
DB = f"{DATA}/db"

FRAG_UCSC_BGZ = f"{DATA}/it_fragments.ucsc.tsv.gz"            # Stage 3c-prep (14a.v2): IT-only, .gz
IT_BARCODES = f"{DATA}/it_barcodes.txt"                       # Stage 1 v2 (49,283 barcodes)
PEAKS_FILTERED = f"{WORK}/it_summits_501_filtered.bed"        # Stage 3a v2 consensus peaks
HG38_BLACKLIST = f"{DB}/hg38-blacklist.v2.bed"                # Stage 2 (shared)

CISTOPIC_OBJ = f"{WORK}/cistopic_obj_model.pkl"               # output: cisTopic obj + LDA model
MODELS_DIR = f"{WORK}/lda_models"
REGION_SETS_DIR = f"{WORK}/region_sets"
TOPICS_SUBFOLDER = "topics_top_3k"

PROJECT = "it"                                                # → barcode suffix ___it (Stage 4/6)
N_TOPICS = [20, 30, 40]          # reduced LDA sweep (user request); auto-select by coherence
N_ITER = 150
N_CPU = 3                        # LDA (run_cgs_models): one model per worker across the sweep.
                                 # Memory is the limiter — each worker holds the ~49k×~150-200k
                                 # binary matrix; drop this (e.g. run in two batches) if RAM-constrained.
N_CPU_CISTOPIC = 1               # cisTopic build join: serial avoids the pyranges/ray parallel-join
                                 # NaN in regions.join(fragments).Score.astype(int32). ~30min-2h at scale.
RANDOM_STATE = 555
NTOP = 3000                      # top regions per binarized topic

os.makedirs(WORK, exist_ok=True)
os.makedirs(MODELS_DIR, exist_ok=True)

# fail-fast on missing inputs
for f in (FRAG_UCSC_BGZ, IT_BARCODES, PEAKS_FILTERED, HG38_BLACKLIST):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics

# ─────────────────────────────────────────────
# 1. Build the cisTopic object (restricted to IT barcodes)
# ─────────────────────────────────────────────
valid_bc = open(IT_BARCODES).read().split()
print(f"[Stage 3c v2] building cisTopic object on {len(valid_bc):,} IT barcodes")
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=FRAG_UCSC_BGZ,
    path_to_regions=PEAKS_FILTERED,
    path_to_blacklist=HG38_BLACKLIST,
    valid_bc=valid_bc,           # restrict to IT cells (file is already IT-only; belt-and-suspenders)
    n_cpu=N_CPU_CISTOPIC,        # serial join (see note above)
    project=PROJECT,
    # do NOT pass split_pattern — it creates duplicate cells (v1astro note)
)
print(f"  ATAC binary matrix: {cistopic_obj.binary_matrix.shape}")

# ─────────────────────────────────────────────
# 2. LDA topic modelling (sweep + auto-select)
# ─────────────────────────────────────────────
print(f"[Stage 3c v2] LDA sweep n_topics={N_TOPICS}, n_iter={N_ITER}")
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
# 3. Binarize topics → per-topic region sets (BED)  [topics only — no DAR, per plan/04]
# ─────────────────────────────────────────────
print(f"[Stage 3c v2] binarizing topics (ntop={NTOP})")
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
print("[Stage 3c v2] complete.")
