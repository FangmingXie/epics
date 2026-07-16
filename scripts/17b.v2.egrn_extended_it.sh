#!/usr/bin/env bash
# Stage 5 extended variant (v2) — eGRN using the ORTHOLOGY-EXTENDED cistromes (all-IT scope, plan/04).
# New copy of 17b_egrn_extended.sh; work_it/ paths, wang25_it DB. Reuses the Stage 5a/5b adjacencies
# (tf_to_gene_adj.tsv, region_to_gene_adj.tsv); only the eGRN assembly + AUCell differ (--is_extended,
# cistromes_extended.h5ad). Broader TF coverage at lower confidence than the direct run. Separate outputs.
# Run: bash scripts/17b.v2.egrn_extended_it.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work_it                                      # NEW isolated work dir
DB=${DATA}/db
TMP=${WORK}/tmp

TF2G=${WORK}/tf_to_gene_adj.tsv                            # Stage 5a v2 (reused)
R2G=${WORK}/region_to_gene_adj.tsv                       # Stage 5b v2 (reused)
CISTROMES_EXTENDED=${WORK}/cistromes_extended.h5ad        # Stage 4d v2 (orthology-extended)
CTX_RANKINGS=${DB}/wang25_it.regions_vs_motifs.rankings.feather   # Stage 3b v2
MUDATA=${WORK}/ACC_GEX.h5mu                              # Stage 4a v2

EREGULON_EXT=${WORK}/eRegulon_extended.tsv              # output
AUCELL_EXT=${WORK}/AUCell_extended.h5mu                 # output

N_CPU=16
SEED=666
SP="conda run --no-capture-output -n epics scenicplus"

for f in "${TF2G}" "${R2G}" "${CISTROMES_EXTENDED}" "${CTX_RANKINGS}" "${MUDATA}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
mkdir -p "${TMP}"

echo "[$(date)] eGRN (extended) -> ${EREGULON_EXT}"
${SP} grn_inference eGRN \
    --is_extended \
    --TF_to_gene_adj_fname "${TF2G}" \
    --region_to_gene_adj_fname "${R2G}" \
    --cistromes_fname "${CISTROMES_EXTENDED}" --ranking_db_fname "${CTX_RANKINGS}" \
    --eRegulon_out_fname "${EREGULON_EXT}" \
    --temp_dir "${TMP}" \
    --order_regions_to_genes_by importance --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --min_regions_per_gene 0 --rho_threshold 0.05 --min_target_genes 10 \
    --n_cpu ${N_CPU} --seed ${SEED}

echo "[$(date)] AUCell (extended) -> ${AUCELL_EXT}"
${SP} grn_inference AUCell \
    --eRegulon_fname "${EREGULON_EXT}" \
    --multiome_mudata_fname "${MUDATA}" \
    --aucell_out_fname "${AUCELL_EXT}" \
    --n_cpu ${N_CPU}

echo "[$(date)] extended eGRN (v2) complete:"
wc -l "${EREGULON_EXT}"
ls -lh "${AUCELL_EXT}"
