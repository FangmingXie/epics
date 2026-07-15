#!/usr/bin/env bash
# Stage 4 — SCENIC+ prepare_data + motif enrichment → MuData, search space, cistromes.
#   4a prepare_GEX_ACC          -> ACC_GEX.h5mu   (multiome MuData; per-cell, barcode-matched)
#   4b search_spance            -> search_space.tsv
#   4c motif_enrichment_cistarget -> ctx_results.hdf5  (topic region-sets vs custom cisTarget DB)
#   4d prepare_menr             -> cistromes_{direct,extended}.h5ad (+ tf_names); marks TFs in MuData
# Run: bash scripts/16_prepare_menr.sh   (all steps use the epics env internally)
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
WORK=${DATA}/work
DB=${DATA}/db
TMP=${WORK}/tmp

CISTOPIC_OBJ=${WORK}/cistopic_obj_model.pkl                # Stage 3c
L23_GEX_H5AD=${DATA}/v1_l23_gex.h5ad                       # Stage 1 (raw counts in .X, .raw cleared)
GENE_ANNOT=${DB}/hg38_gene_annotation.tsv                  # Stage 2 (UCSC, primary chroms)
HG38_CHROMSIZES=${DB}/hg38.chrom.sizes                     # Stage 2 (2-col: chrom<TAB>size)
HG38_CHROMSIZES_SP=${DB}/hg38_chromsizes_scenicplus.tsv    # search_spance format: Chromosome/Start/End header
REGION_SETS=${WORK}/region_sets                           # parent: motif_enrichment iterates its SUBFOLDERS
                                                          #   (topics_top_3k/*.bed → keys "topics_top_3k_TopicN")
CTX_RANKINGS=${DB}/wang25_l23.regions_vs_motifs.rankings.feather   # Stage 3b
MOTIF_TBL=/data/qlyu/data/common/aertslab_motif_collection/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl

MUDATA=${WORK}/ACC_GEX.h5mu                                # 4a output
SEARCH_SPACE=${WORK}/search_space.tsv                      # 4b output
CTX_RESULT=${WORK}/ctx_results.hdf5                        # 4c output (pycistarget hdf5)
CISTROMES_DIRECT=${WORK}/cistromes_direct.h5ad             # 4d output (direct motif→TF)
CISTROMES_EXTENDED=${WORK}/cistromes_extended.h5ad         # 4d output (orthology-extended)
TF_NAMES_MENR=${WORK}/tf_names_menr.txt                    # 4d output (cistrome TFs; Stage 5 uses data/tf_names.txt)

PROJECT_SUFFIX=___l23                                      # cisTopic appends this to ATAC barcodes
N_CPU=8
SP="conda run --no-capture-output -n epics scenicplus"

# --- fail-fast ---
for f in "${CISTOPIC_OBJ}" "${L23_GEX_H5AD}" "${GENE_ANNOT}" "${HG38_CHROMSIZES}" "${CTX_RANKINGS}" "${MOTIF_TBL}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -d "${REGION_SETS}" ]] || { echo "ERROR: region-set dir missing: ${REGION_SETS}" >&2; exit 1; }
mkdir -p "${TMP}"

# =============================================================================
# 4a — MuData (multiome, per-cell). RNA barcodes get the ___l23 suffix to match ATAC.
#      --do_not_use_raw_for_GEX_anndata: Stage 1 put raw counts in .X and cleared .raw.
# =============================================================================
echo "[$(date)] 4a prepare_GEX_ACC -> ${MUDATA}"
${SP} prepare_data prepare_GEX_ACC \
    --cisTopic_obj_fname "${CISTOPIC_OBJ}" \
    --GEX_anndata_fname "${L23_GEX_H5AD}" \
    --out_file "${MUDATA}" \
    --do_not_use_raw_for_GEX_anndata \
    --bc_transform_func "lambda x: x + '${PROJECT_SUFFIX}'"

# =============================================================================
# 4b — region↔gene search space (explicit ±1kb–150kb per open risk #4)
# search_spance reads chromsizes via pd.read_table expecting Chromosome/Start/End
# columns (with header) — build that from the 2-col chrom.sizes (primary chroms).
# =============================================================================
printf 'Chromosome\tStart\tEnd\n' > "${HG38_CHROMSIZES_SP}"
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XYM])$/ {print $1, 0, $2}' "${HG38_CHROMSIZES}" >> "${HG38_CHROMSIZES_SP}"
echo "[$(date)] 4b search_spance -> ${SEARCH_SPACE}  (chromsizes: $(($(wc -l < "${HG38_CHROMSIZES_SP}") - 1)) chroms)"
${SP} prepare_data search_spance \
    --multiome_mudata_fname "${MUDATA}" \
    --gene_annotation_fname "${GENE_ANNOT}" \
    --chromsizes_fname "${HG38_CHROMSIZES_SP}" \
    --out_fname "${SEARCH_SPACE}" \
    --upstream 1000 150000 \
    --downstream 1000 150000

# =============================================================================
# 4c — motif enrichment: topic region-sets vs the custom hg38 cisTarget DB
# =============================================================================
echo "[$(date)] 4c motif_enrichment_cistarget -> ${CTX_RESULT}"
${SP} grn_inference motif_enrichment_cistarget \
    --region_set_folder "${REGION_SETS}" \
    --cistarget_db_fname "${CTX_RANKINGS}" \
    --output_fname_cistarget_result "${CTX_RESULT}" \
    --temp_dir "${TMP}" --species homo_sapiens \
    --path_to_motif_annotations "${MOTIF_TBL}" \
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

echo "[$(date)] Stage 4 complete:"
ls -lh "${MUDATA}" "${SEARCH_SPACE}" "${CTX_RESULT}" "${CISTROMES_DIRECT}" "${CISTROMES_EXTENDED}" "${TF_NAMES_MENR}"
echo "  search-space rows: $(wc -l < "${SEARCH_SPACE}")"
echo "  cistrome TFs (menr): $(wc -l < "${TF_NAMES_MENR}")"