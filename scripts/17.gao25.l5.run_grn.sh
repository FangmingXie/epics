#!/usr/bin/env bash
# Stage 5 (gao25 · l5) — SCENIC+ GRN inference → eRegulons + per-cell activity (mouse L5).
# Fork of scripts/17.gao25.run_grn.sh for the L5 layer (plan/06). Params are species-agnostic and
# unchanged; only the work/DB paths point at the l5 namespace.
#   5a TF_to_gene      -> tf_to_gene_adj.tsv     (GBM co-expression)
#   5b region_to_gene  -> region_to_gene_adj.tsv (GBM importance + Spearman)
#   5c eGRN            -> eRegulon.tsv           (direct cistromes)
#   5d AUCell          -> AUCell.h5mu            (per-cell eRegulon activity, for Stage 6)
# Run: bash scripts/17.gao25.l5.run_grn.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/gao25/l5
WORK=${DATA}/work
DB=${DATA}/db
TMP=${WORK}/tmp

MUDATA=${WORK}/ACC_GEX.h5mu                                # Stage 4a
TF_NAMES=${WORK}/tf_names_menr.txt                        # Stage 4d (cistrome TFs)
SEARCH_SPACE=${WORK}/search_space.tsv                     # Stage 4b
CISTROMES_DIRECT=${WORK}/cistromes_direct.h5ad            # Stage 4d (direct motif→TF)
CTX_RANKINGS=${DB}/gao25_l5.regions_vs_motifs.rankings.feather   # Stage 3b (per-layer DB)

TF2G=${WORK}/tf_to_gene_adj.tsv                           # 5a output
R2G=${WORK}/region_to_gene_adj.tsv                       # 5b output
EREGULON=${WORK}/eRegulon.tsv                            # 5c output (final eRegulons)
AUCELL_OUT=${WORK}/AUCell.h5mu                           # 5d output (per-cell activity)

N_CPU=16
SEED=666
SP="conda run --no-capture-output -n epics scenicplus"

# --- fail-fast ---
for f in "${MUDATA}" "${TF_NAMES}" "${SEARCH_SPACE}" "${CISTROMES_DIRECT}" "${CTX_RANKINGS}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
mkdir -p "${TMP}"

# =============================================================================
# 5a — TF → gene (GBM co-expression)
# =============================================================================
echo "[$(date)] 5a TF_to_gene ($(wc -l < "${TF_NAMES}") TFs) -> ${TF2G}"
${SP} grn_inference TF_to_gene \
    --multiome_mudata_fname "${MUDATA}" --tf_names "${TF_NAMES}" \
    --temp_dir "${TMP}" --out_tf_to_gene_adjacencies "${TF2G}" \
    --method GBM --n_cpu ${N_CPU} --seed ${SEED}

# =============================================================================
# 5b — region → gene (GBM importance + Spearman correlation)
# =============================================================================
echo "[$(date)] 5b region_to_gene -> ${R2G}"
${SP} grn_inference region_to_gene \
    --multiome_mudata_fname "${MUDATA}" --search_space_fname "${SEARCH_SPACE}" \
    --temp_dir "${TMP}" --out_region_to_gene_adjacencies "${R2G}" \
    --importance_scoring_method GBM --correlation_scoring_method SR --n_cpu ${N_CPU}

# =============================================================================
# 5c — eGRN (direct cistromes; full-pipeline params)
# =============================================================================
echo "[$(date)] 5c eGRN -> ${EREGULON}"
${SP} grn_inference eGRN \
    --TF_to_gene_adj_fname "${TF2G}" \
    --region_to_gene_adj_fname "${R2G}" \
    --cistromes_fname "${CISTROMES_DIRECT}" --ranking_db_fname "${CTX_RANKINGS}" \
    --eRegulon_out_fname "${EREGULON}" \
    --temp_dir "${TMP}" \
    --order_regions_to_genes_by importance --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --min_regions_per_gene 0 --rho_threshold 0.05 --min_target_genes 10 \
    --n_cpu ${N_CPU} --seed ${SEED}

# =============================================================================
# 5d — AUCell (per-cell eRegulon activity)
# =============================================================================
echo "[$(date)] 5d AUCell -> ${AUCELL_OUT}"
${SP} grn_inference AUCell \
    --eRegulon_fname "${EREGULON}" \
    --multiome_mudata_fname "${MUDATA}" \
    --aucell_out_fname "${AUCELL_OUT}" \
    --n_cpu ${N_CPU}

echo "[$(date)] Stage 5 gao25·l5 complete:"
wc -l "${TF2G}" "${R2G}" "${EREGULON}"
ls -lh "${AUCELL_OUT}"
