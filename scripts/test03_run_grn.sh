#!/usr/bin/env bash
# Tier 2b/2c/2d: run the SCENIC+ GRN-inference steps on the mini MuData.
#   2b TF_to_gene      -> tf_to_gene_adj.tsv
#   2c region_to_gene  -> region_to_gene_adj.tsv   (reuses the reference search space)
#   2d eGRN            -> eRegulon.tsv             (reuses reference cistromes + ranking db)
# Run: bash scripts/test03_run_grn.sh
set -euo pipefail

# ---------------------------------------------------------------- paths / params
ENV=epics
RUN="conda run --no-capture-output -n ${ENV} python -u"   # unbuffered (CLAUDE.md)
SP="conda run --no-capture-output -n ${ENV} scenicplus"

OUT_DIR=/data/qlyu/project/epics/data/test_run
MU=${OUT_DIR}/ACC_GEX.small.h5mu
TF_NAMES=${OUT_DIR}/tf_names.txt
TMP=${OUT_DIR}/tmp
mkdir -p "${TMP}"

# reference inputs reused read-only from the prior astro run
REF=/home/qlyu/mydata/data/v1_astro/t2p1/spout
SEARCH_SPACE=${REF}/search_space.tsv
CISTROMES=${REF}/cistromes_direct.h5ad
RANKING_DB=/home/qlyu/mydata/data/v1_astro/t2/spin/t2.regions_vs_motifs.rankings.feather

TF2G=${OUT_DIR}/tf_to_gene_adj.tsv
R2G=${OUT_DIR}/region_to_gene_adj.tsv
EREGULON=${OUT_DIR}/eRegulon.tsv

N_CPU=8

echo "=== 2b TF_to_gene ==="
${SP} grn_inference TF_to_gene \
    --multiome_mudata_fname "${MU}" \
    --tf_names "${TF_NAMES}" \
    --temp_dir "${TMP}" \
    --out_tf_to_gene_adjacencies "${TF2G}" \
    --method GBM --n_cpu ${N_CPU} --seed 666

echo "=== 2c region_to_gene ==="
${SP} grn_inference region_to_gene \
    --multiome_mudata_fname "${MU}" \
    --search_space_fname "${SEARCH_SPACE}" \
    --temp_dir "${TMP}" \
    --out_region_to_gene_adjacencies "${R2G}" \
    --importance_scoring_method GBM --correlation_scoring_method SR --n_cpu ${N_CPU}

echo "=== 2d eGRN (relaxed params to guarantee signal on the small set) ==="
${SP} grn_inference eGRN \
    --TF_to_gene_adj_fname "${TF2G}" \
    --region_to_gene_adj_fname "${R2G}" \
    --cistromes_fname "${CISTROMES}" \
    --ranking_db_fname "${RANKING_DB}" \
    --eRegulon_out_fname "${EREGULON}" \
    --temp_dir "${TMP}" \
    --gsea_n_perm 100 \
    --min_target_genes 3 \
    --min_regions_per_gene 0 \
    --n_cpu ${N_CPU} --seed 666

echo "=== outputs ==="
wc -l "${TF2G}" "${R2G}" "${EREGULON}"
echo "GRN RUN OK"
