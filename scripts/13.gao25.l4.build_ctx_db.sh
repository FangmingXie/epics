#!/usr/bin/env bash
# Stage 3b (gao25 · l4) — Build a custom cisTarget motif database on the L4/5-accessible mouse peaks.
# Fork of scripts/13.gao25.build_ctx_db.sh for the L4 layer (plan/06). The build is species-neutral;
# only the peak BED and the output DB prefix change. Built ON the L4/5 peak BED (Stage 3a) with a NEW
# prefix gao25_l4 in the per-layer db dir (do NOT share the L2/3 DB — each layer builds on its own peaks).
# MUST run inside the `epics` env (bedtools + cbust + pyarrow):
#   conda run --no-capture-output -n epics bash scripts/13.gao25.l4.build_ctx_db.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/gao25/l4
WORK=${DATA}/work
DB=${DATA}/db

DB_BUILDER=/data/qlyu/code/trial/create_cisTarget_databases
MM10_FASTA=/data/qlyu/data/common/mm10/mm10.fa
MM10_CHROMSIZES=/data/qlyu/data/common/mm10/mm10.chrom.sizes
L4_PEAKS_BED=${WORK}/vis_l4_peaks.bed                  # Stage 3a output (chr\tstart\tend)
MOTIF_DIR=/data/qlyu/data/common/aertslab_motif_collection/v10nr_clust_public/singletons

OUT_MOTIF_LIST=${DB}/motifs.txt                        # list of motif .cb filenames
OUT_FASTA=${DB}/l4_regions_1kb_bg.fa                   # padded-background region FASTA
CTX_DB_PREFIX=${DB}/gao25_l4                           # → *.regions_vs_motifs.{rankings,scores}.feather

BG_PADDING=1000
N_CPU=48                                                # cbust scoring cores (machine has 128); scales ~linearly

# --- fail-fast checks (PATH tools come from the epics env) ---
command -v bedtools >/dev/null || { echo "ERROR: bedtools not on PATH — run under: conda run -n epics bash $0" >&2; exit 1; }
command -v cbust    >/dev/null || { echo "ERROR: cbust (Cluster-Buster) not on PATH — run under the epics env." >&2; exit 1; }
for f in "${MM10_FASTA}" "${MM10_CHROMSIZES}" "${L4_PEAKS_BED}"; do
    [[ -s "${f}" ]] || { echo "ERROR: required input missing/empty: ${f}" >&2; exit 1; }
done
[[ -d "${MOTIF_DIR}" ]] || { echo "ERROR: motif dir missing: ${MOTIF_DIR}" >&2; exit 1; }
mkdir -p "${DB}"

# --- motif list (one .cb filename per line) ---
ls "${MOTIF_DIR}" > "${OUT_MOTIF_LIST}"
echo "[$(date)] motifs listed: $(wc -l < "${OUT_MOTIF_LIST}") -> ${OUT_MOTIF_LIST}"

# --- padded-background FASTA for the L4/5 regions ---
echo "[$(date)] building padded-bg FASTA (${BG_PADDING} bp) for $(wc -l < "${L4_PEAKS_BED}") regions -> ${OUT_FASTA}"
"${DB_BUILDER}/create_fasta_with_padded_bg_from_bed.sh" \
    "${MM10_FASTA}" "${MM10_CHROMSIZES}" "${L4_PEAKS_BED}" \
    "${OUT_FASTA}" "${BG_PADDING}" yes

# --- score regions × motifs → rankings/scores feathers ---
echo "[$(date)] scoring regions × motifs (t=${N_CPU}) -> ${CTX_DB_PREFIX}.*"
"${DB_BUILDER}/create_cistarget_motif_databases.py" \
    -f "${OUT_FASTA}" \
    -M "${MOTIF_DIR}" -m "${OUT_MOTIF_LIST}" \
    -o "${CTX_DB_PREFIX}" --bgpadding "${BG_PADDING}" -t "${N_CPU}"

echo "[$(date)] Stage 3b gao25·l4 complete:"
ls -lh "${CTX_DB_PREFIX}".regions_vs_motifs.* 2>/dev/null
