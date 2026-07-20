#!/usr/bin/env bash
# Driver (gao25) — run the mouse DevVIS L2/3 SCENIC+ eGRN pipeline (Stages 1–7) in order.
# Stage 0 (both multiome h5ads) is done. Each stage streams to its own log; set -e halts on failure.
# Run in the background:  bash scripts/run_all.gao25.sh
set -euo pipefail

PROJ=/data/qlyu/project/epics
cd "${PROJ}"
LOGS=logs
mkdir -p "${LOGS}"
EPICS="conda run --no-capture-output -n epics"

step () {  # step <logname> <cmd...>
    local name="$1"; shift
    echo "===== [$(date)] START ${name} ====="
    "$@" 2>&1 | tee "${LOGS}/${name}.log"
    echo "===== [$(date)] DONE  ${name} ====="
}

# Stage 1 — select cortical L2/3 + multiome barcodes
step 10.gao25  ${EPICS} python -u scripts/10.gao25.select_l23.py

# Stage 2 — mouse GRCm38 references (gene annotation + chromsizes_sp)
step 11.gao25  bash scripts/11.gao25.refs.sh

# Stage 3a — L2/3-accessible peak matrix (BED + regions×cells npz)
step 12.gao25  ${EPICS} python -u scripts/12.gao25.prep_peaks_matrix.py

# Stage 3b — custom mm10 cisTarget DB on the L2/3 peaks (epics env for bedtools/cbust)
step 13.gao25  ${EPICS} bash scripts/13.gao25.build_ctx_db.sh

# Stage 3c — cisTopic FROM matrix + LDA + binarize
step 14.gao25  ${EPICS} python -u scripts/14.gao25.run_pycistopic.py

# Stage 4 — prepare_data + motif enrichment (mus_musculus + MGI)
step 16.gao25  bash scripts/16.gao25.prepare_menr.sh

# Stage 5 — GRN inference (direct + extended)
step 17.gao25   bash scripts/17.gao25.run_grn.sh
step 17b.gao25  bash scripts/17b.gao25.egrn_extended.sh

# Stage 6/7 — QC + organize regulons
step 18.gao25  ${EPICS} python -u scripts/18.gao25.qc.py
step 19.gao25  ${EPICS} python -u scripts/19.gao25.organize_regulons.py

echo "===== [$(date)] ALL STAGES COMPLETE ====="
