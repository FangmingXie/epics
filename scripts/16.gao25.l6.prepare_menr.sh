#!/usr/bin/env bash
# Stage 4 (gao25 · l6) — SCENIC+ prepare_data + motif enrichment → MuData, search space, cistromes.
# Fork of scripts/16.gao25.prepare_menr.sh for the L6 layer (plan/06). Method unchanged; only the gex,
# cisTarget DB prefix, and the ___gao25l6 barcode suffix differ. The mm10 GRCm38 gene annotation and
# scenicplus chromsizes are subclass-independent and REUSED from the shared L2/3 refs (data/gao25/db).
#   4a prepare_GEX_ACC            -> ACC_GEX.h5mu   (multiome MuData; per-cell, barcode-matched)
#   4b search_spance              -> search_space.tsv
#   4c motif_enrichment_cistarget -> ctx_results.hdf5  (topic region-sets vs custom mm10 cisTarget DB)
#   4d prepare_menr               -> cistromes_{direct,extended}.h5ad (+ tf_names); marks TFs in MuData
# Run: bash scripts/16.gao25.l6.prepare_menr.sh   (all steps use the epics env internally)
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/gao25/l6
WORK=${DATA}/work
DB=${DATA}/db
SHARED_DB=${PROJ}/data/gao25/db                           # shared subclass-independent refs (Stage 2)
TMP=${WORK}/tmp

CISTOPIC_OBJ=${WORK}/cistopic_obj_model.pkl                # Stage 3c
L6_GEX_H5AD=${DATA}/vis_l6_gex.h5ad                       # Stage 1 (raw counts in .X, .raw cleared)
GENE_ANNOT=${SHARED_DB}/mm10_gene_annotation.tsv         # shared (GRCm38, UCSC chr names)
MM10_CHROMSIZES_SP=${SHARED_DB}/mm10_chromsizes_scenicplus.tsv  # shared (Chromosome/Start/End header)
REGION_SETS=${WORK}/region_sets                          # parent; motif_enrichment iterates its SUBFOLDERS
CTX_RANKINGS=${DB}/gao25_l6.regions_vs_motifs.rankings.feather   # Stage 3b (per-layer DB)
MOTIF_TBL_MGI=/data/qlyu/data/common/aertslab_motif_collection/v10nr_clust_public/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl

MUDATA=${WORK}/ACC_GEX.h5mu                               # 4a output
SEARCH_SPACE=${WORK}/search_space.tsv                     # 4b output
CTX_RESULT=${WORK}/ctx_results.hdf5                       # 4c output (pycistarget hdf5)
CISTROMES_DIRECT=${WORK}/cistromes_direct.h5ad            # 4d output (direct motif→TF)
CISTROMES_EXTENDED=${WORK}/cistromes_extended.h5ad        # 4d output (orthology-extended)
TF_NAMES_MENR=${WORK}/tf_names_menr.txt                   # 4d output (cistrome TFs; Stage 5 uses this)

PROJECT_SUFFIX=___gao25l6                                 # cisTopic appends this to ATAC barcodes
N_CPU=8
SP="conda run --no-capture-output -n epics scenicplus"

# --- fail-fast ---
for f in "${CISTOPIC_OBJ}" "${L6_GEX_H5AD}" "${GENE_ANNOT}" "${MM10_CHROMSIZES_SP}" "${CTX_RANKINGS}" "${MOTIF_TBL_MGI}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -d "${REGION_SETS}" ]] || { echo "ERROR: region-set dir missing: ${REGION_SETS}" >&2; exit 1; }
mkdir -p "${TMP}"

# =============================================================================
# 4a — MuData (multiome, per-cell). RNA barcodes get the ___gao25l6 suffix to match ATAC.
#      --do_not_use_raw_for_GEX_anndata: Stage 1 put raw counts in .X and cleared .raw.
# =============================================================================
echo "[$(date)] 4a prepare_GEX_ACC -> ${MUDATA}"
${SP} prepare_data prepare_GEX_ACC \
    --cisTopic_obj_fname "${CISTOPIC_OBJ}" \
    --GEX_anndata_fname "${L6_GEX_H5AD}" \
    --out_file "${MUDATA}" \
    --do_not_use_raw_for_GEX_anndata \
    --bc_transform_func "lambda x: x + '${PROJECT_SUFFIX}'"

# =============================================================================
# 4b — region↔gene search space (explicit ±1kb–150kb). GRCm38 annotation guards open risk #1:
#      a non-trivial row count confirms peaks (mm10) and genes (GRCm38) share coordinates.
# =============================================================================
echo "[$(date)] 4b search_spance -> ${SEARCH_SPACE}"
${SP} prepare_data search_spance \
    --multiome_mudata_fname "${MUDATA}" \
    --gene_annotation_fname "${GENE_ANNOT}" \
    --chromsizes_fname "${MM10_CHROMSIZES_SP}" \
    --out_fname "${SEARCH_SPACE}" \
    --upstream 1000 150000 \
    --downstream 1000 150000

# =============================================================================
# 4c — motif enrichment: topic region-sets vs the custom mm10 cisTarget DB (MOUSE: mus_musculus + MGI)
# =============================================================================
echo "[$(date)] 4c motif_enrichment_cistarget -> ${CTX_RESULT}"
${SP} grn_inference motif_enrichment_cistarget \
    --region_set_folder "${REGION_SETS}" \
    --cistarget_db_fname "${CTX_RANKINGS}" \
    --output_fname_cistarget_result "${CTX_RESULT}" \
    --temp_dir "${TMP}" --species mus_musculus \
    --path_to_motif_annotations "${MOTIF_TBL_MGI}" \
    --annotation_version v10nr_clust --n_cpu ${N_CPU}

# =============================================================================
# 4d — cistromes (direct + extended) + TF names; marks TFs in the MuData
# =============================================================================
echo "[$(date)] 4d prepare_menr -> ${CISTROMES_DIRECT} / ${CISTROMES_EXTENDED}"
${SP} prepare_data prepare_menr \
    --paths_to_motif_enrichment_results "${CTX_RESULT}" \
    --multiome_mudata_fname "${MUDATA}" \
    --out_file_tf_names "${TF_NAMES_MENR}" \
    --out_file_direct_annotation "${CISTROMES_DIRECT}" \
    --out_file_extended_annotation "${CISTROMES_EXTENDED}" \
    --direct_annotation Direct_annot \
    --extended_annotation Orthology_annot

echo "[$(date)] Stage 4 gao25·l6 complete:"
ls -lh "${MUDATA}" "${SEARCH_SPACE}" "${CTX_RESULT}" "${CISTROMES_DIRECT}" "${CISTROMES_EXTENDED}" "${TF_NAMES_MENR}"
echo "  search-space rows: $(wc -l < "${SEARCH_SPACE}")"
echo "  cistrome TFs (menr): $(wc -l < "${TF_NAMES_MENR}")"
