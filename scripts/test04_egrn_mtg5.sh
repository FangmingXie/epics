#!/usr/bin/env bash
# Diagnostic — re-run ONLY the eGRN assembly with --min_target_genes 5 (default was 10),
# to check whether weak regulons (e.g. MEIS2) emerge below the 10-target threshold.
# Reuses the Stage 5a/5b adjacencies + cistromes; writes NEW _mtg5 outputs (no overwrite).
# AUCell is intentionally skipped (not needed to inspect regulon membership).
# Run: bash scripts/test04_egrn_mtg5.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work
DB=${DATA}/db
TMP=${WORK}/tmp

TF2G=${WORK}/tf_to_gene_adj.tsv                            # Stage 5a (reused)
R2G=${WORK}/region_to_gene_adj.tsv                        # Stage 5b (reused)
CISTROMES_DIRECT=${WORK}/cistromes_direct.h5ad            # Stage 4d
CISTROMES_EXTENDED=${WORK}/cistromes_extended.h5ad        # Stage 4d
CTX_RANKINGS=${DB}/wang25_l23.regions_vs_motifs.rankings.feather   # Stage 3b

EREGULON_DIRECT_MTG5=${WORK}/eRegulon_mtg5.tsv            # output (direct, min_target_genes=5)
EREGULON_EXT_MTG5=${WORK}/eRegulon_extended_mtg5.tsv      # output (extended, min_target_genes=5)

MIN_TARGET_GENES=5
N_CPU=8
SEED=666
SP="conda run --no-capture-output -n epics scenicplus"

for f in "${TF2G}" "${R2G}" "${CISTROMES_DIRECT}" "${CISTROMES_EXTENDED}" "${CTX_RANKINGS}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
mkdir -p "${TMP}"

# --- direct ---
echo "[$(date)] eGRN direct (min_target_genes=${MIN_TARGET_GENES}) -> ${EREGULON_DIRECT_MTG5}"
${SP} grn_inference eGRN \
    --TF_to_gene_adj_fname "${TF2G}" \
    --region_to_gene_adj_fname "${R2G}" \
    --cistromes_fname "${CISTROMES_DIRECT}" --ranking_db_fname "${CTX_RANKINGS}" \
    --eRegulon_out_fname "${EREGULON_DIRECT_MTG5}" \
    --temp_dir "${TMP}" \
    --order_regions_to_genes_by importance --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --min_regions_per_gene 0 --rho_threshold 0.05 --min_target_genes ${MIN_TARGET_GENES} \
    --n_cpu ${N_CPU} --seed ${SEED}

# --- extended ---
echo "[$(date)] eGRN extended (min_target_genes=${MIN_TARGET_GENES}) -> ${EREGULON_EXT_MTG5}"
${SP} grn_inference eGRN \
    --is_extended \
    --TF_to_gene_adj_fname "${TF2G}" \
    --region_to_gene_adj_fname "${R2G}" \
    --cistromes_fname "${CISTROMES_EXTENDED}" --ranking_db_fname "${CTX_RANKINGS}" \
    --eRegulon_out_fname "${EREGULON_EXT_MTG5}" \
    --temp_dir "${TMP}" \
    --order_regions_to_genes_by importance --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --min_regions_per_gene 0 --rho_threshold 0.05 --min_target_genes ${MIN_TARGET_GENES} \
    --n_cpu ${N_CPU} --seed ${SEED}

echo "[$(date)] diagnostic eGRN complete:"
wc -l "${EREGULON_DIRECT_MTG5}" "${EREGULON_EXT_MTG5}"
