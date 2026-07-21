#!/usr/bin/env bash
# Driver (gao25 · l6) — run the mouse DevVIS L6 SCENIC+ eGRN pipeline (Stages 1, 3–7) in order.
# Fork of scripts/run_all.gao25.sh for the L6 layer (plan/06). Stage 2 (mouse GRCm38 refs) is DROPPED —
# the gene annotation + chromsizes are subclass-independent and reused from the shared L2/3 refs
# (data/gao25/db). Each stage streams to its own log; set -e halts on failure.
# Run in the background:  bash scripts/run_all.gao25.l6.sh
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

# Stage 1 — select cortical L6 + multiome barcodes
step 10.gao25.l6  ${EPICS} python -u scripts/10.gao25.l6.select.py

# Stage 2 — SKIPPED: shared mm10 GRCm38 refs reused from the L2/3 run (data/gao25/db)

# Stage 3a — L6-accessible peak matrix (BED + regions×cells npz)
step 12.gao25.l6  ${EPICS} python -u scripts/12.gao25.l6.prep_peaks.py

# Stage 3b — custom mm10 cisTarget DB on the L6 peaks (epics env for bedtools/cbust)
step 13.gao25.l6  ${EPICS} bash scripts/13.gao25.l6.build_ctx_db.sh

# Stage 3c — cisTopic FROM matrix + LDA + binarize
step 14.gao25.l6  ${EPICS} python -u scripts/14.gao25.l6.run_pycistopic.py

# Stage 4 — prepare_data + motif enrichment (mus_musculus + MGI)
step 16.gao25.l6  bash scripts/16.gao25.l6.prepare_menr.sh

# Stage 5 — GRN inference (direct + extended)
step 17.gao25.l6   bash scripts/17.gao25.l6.run_grn.sh
step 17b.gao25.l6  bash scripts/17b.gao25.l6.egrn_extended.sh

# Stage 6/7 — QC + organize regulons
step 18.gao25.l6  ${EPICS} python -u scripts/18.gao25.l6.qc.py
step 19.gao25.l6  ${EPICS} python -u scripts/19.gao25.l6.organize_regulons.py

echo "===== [$(date)] ALL STAGES COMPLETE ====="
