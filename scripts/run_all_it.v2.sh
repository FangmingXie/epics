#!/usr/bin/env bash
# Driver (v2) — run the all-IT SCENIC+ eGRN pipeline (Stages 3a–7) in order.
# Stage 1 (10.v2.select_it.py) is assumed already done (it_gex.h5ad + it_barcodes.txt exist).
# Stages 0/2/2b are shared and skipped. Each stage streams to its own log; set -e halts on failure.
# Run in the background:  bash scripts/run_all_it.v2.sh
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

# Stage 3a — pooled peak calling (macs3 env inside) + consensus
step 12.v2   bash scripts/12.v2.call_peaks_it.sh
step 12b.v2  ${EPICS} python -u scripts/12b.v2.consensus_peaks_it.py

# Stage 3b — custom cisTarget DB (needs epics env for bedtools/cbust)
step 13.v2   ${EPICS} bash scripts/13.v2.build_ctx_db_it.sh

# Stage 3c — IT-only fragments, then cisTopic + LDA
step 14a.v2  bash scripts/14a.v2.extract_it_fragments.sh
step 14.v2   ${EPICS} python -u scripts/14.v2.run_pycistopic_it.py

# Stage 4 — prepare_data + motif enrichment (spawns epics env internally)
step 16.v2   bash scripts/16.v2.prepare_menr_it.sh

# Stage 5 — GRN inference (direct + extended)
step 17.v2   bash scripts/17.v2.run_grn_it.sh
step 17b.v2  bash scripts/17b.v2.egrn_extended_it.sh

# Stage 6/7 — QC + organize regulons
step 18.v2   ${EPICS} python -u scripts/18.v2.qc_it.py
step 19.v2   ${EPICS} python -u scripts/19.v2.organize_regulons_it.py

echo "===== [$(date)] ALL STAGES COMPLETE ====="
