#!/usr/bin/env python
"""Stage 3a consensus (v2) — non-overlapping consensus peaks from the all-IT MACS3 narrowPeak.

New copy of 12b_consensus_peaks.py for the all-IT scope (plan/04). Same pycisTopic
get_consensus_peaks logic (re-center each summit to a fixed 501 bp, iteratively drop less
significant overlapping peaks, apply blacklist); single-entry narrow_peaks_dict {LABEL: PyRanges}.
Only paths and LABEL change; output goes to the isolated work_it/ dir.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/12b.v2.consensus_peaks_it.py
"""
import os
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks

# --- paths (capitalized per CLAUDE.md) ---
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/wang25"
WORK = f"{DATA}/work_it"                                    # NEW isolated work dir
DB = f"{DATA}/db"

NARROWPEAK = f"{WORK}/it_peaks.narrowPeak"                  # Stage 3a MACS3 output (12.v2)
HG38_CHROMSIZES = f"{DB}/hg38.chrom.sizes"                  # Stage 2 (UCSC, shared)
HG38_BLACKLIST = f"{DB}/hg38-blacklist.v2.bed"             # Stage 2 (shared)
PEAKS_FILTERED = f"{WORK}/it_summits_501_filtered.bed"     # OUTPUT: consensus peaks (consumed by 3b/3c)

PEAK_HALF_WIDTH = 250                                       # → 501 bp
LABEL = "it"

# narrowPeak columns (MACS3); avoid naming col6 "Strand" so PyRanges stays unstranded.
NP_COLS = ["Chromosome", "Start", "End", "Name", "Score", "dot",
           "FC_summit", "neglog10_pval", "neglog10_qval", "Summit"]

for f in (NARROWPEAK, HG38_CHROMSIZES, HG38_BLACKLIST):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")

print(f"[consensus] reading {NARROWPEAK}")
np_df = pd.read_csv(NARROWPEAK, sep="\t", header=None, names=NP_COLS)
print(f"  raw MACS3 peaks: {len(np_df):,}")
narrow_peaks_dict = {LABEL: pr.PyRanges(np_df)}

# chromsizes as a PyRanges (Chromosome, Start, End)
cs = pd.read_csv(HG38_CHROMSIZES, sep="\t", header=None, names=["Chromosome", "End"])
cs["Start"] = 0
chromsizes = pr.PyRanges(cs[["Chromosome", "Start", "End"]])

print(f"[consensus] get_consensus_peaks (peak_half_width={PEAK_HALF_WIDTH}, +blacklist)")
consensus = get_consensus_peaks(
    narrow_peaks_dict,
    peak_half_width=PEAK_HALF_WIDTH,
    chromsizes=chromsizes,
    path_to_blacklist=HG38_BLACKLIST,
)
cons_df = consensus.df.sort_values(["Chromosome", "Start"]).reset_index(drop=True)
print(f"  consensus peaks: {len(cons_df):,}")

widths = (cons_df["End"] - cons_df["Start"]).value_counts()
print(f"  width distribution: {dict(widths.head())}")

cons_df[["Chromosome", "Start", "End", "Name"]].to_csv(
    PEAKS_FILTERED, sep="\t", header=False, index=False
)
print(f"  wrote {PEAKS_FILTERED}")
print("[consensus] done.")
