#!/usr/bin/env python
"""Stage 6 (gao25 · l5) — QC / developmental dynamics for the mouse L5 eGRN.

Fork of scripts/18.gao25.qc.py for the L5 layer (plan/06). A SINGLE subclass (L5 IT CTX Glut), so the
analysis axes are DEVELOPMENT (`age_label`) and REGION (`roi`, a covariate). Only the paths, barcode
suffix, and the layer-specific positive-control markers differ from L2/3.
Outputs:
  1. eRegulon summary table (TF, direction, #target genes/regions/triplets) -> CSV + top-sizes plot.
  2. Mouse positive-control markers (mixed-case) that surfaced as eRegulons — headline Fezf2 / Bcl11b,
     the canonical deep-layer (L5) identity TFs.
  3. AUCell activity z-scored across `age_label` (developmental order) -> developmental cascade CSV+PNG.
     Caveat: near-zero embryonic N; the real signal is the postnatal P0→P58 trajectory.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/18.gao25.l5.qc.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import mudata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- paths (capitalized per CLAUDE.md) ---
PROJ = "/data/qlyu/project/epics"
DATA = f"{PROJ}/data/gao25/l5"
WORK = f"{DATA}/work"
QC = f"{WORK}/qc"

EREGULON = f"{WORK}/eRegulon.tsv"                 # Stage 5c
AUCELL = f"{WORK}/AUCell.h5mu"                    # Stage 5d
L5_GEX_H5AD = f"{DATA}/vis_l5_gex.h5ad"          # Stage 1 (has roi / age_label / sex / donor)

SUMMARY_CSV = f"{QC}/eRegulon_summary.csv"
STAGE_CSV = f"{QC}/eRegulon_activity_by_age.csv"
REGION_CSV = f"{QC}/eRegulon_activity_by_roi.csv"
FIG_SIZES = f"{QC}/eRegulon_target_sizes.png"
FIG_CASCADE = f"{QC}/eRegulon_developmental_cascade.png"

STAGE_COL, REGION_COL = "age_label", "roi"
BC_SUFFIX = "___gao25l5"                          # AUCell cell names carry the cisTopic project suffix
# developmental order (plan): near-zero embryonic; signal is the postnatal P0→P58 trajectory
AGE_ORDER = ["E15.5", "E16", "E17", "E18", "P0", "P2", "P4", "P5", "P8", "P9", "P11", "P14", "P56", "P58"]
# mouse L5 positive-control marker TFs (mixed-case); headline Fezf2 / Bcl11b (canonical L5 identity TFs)
KNOWN_TFS = ["Fezf2", "Bcl11b", "Etv1", "Foxp2", "Sox5", "Mef2c"]

for f in (EREGULON, AUCELL, L5_GEX_H5AD):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")
os.makedirs(QC, exist_ok=True)

# ─────────────────────────────────────────────
# 1. eRegulon summary table
# ─────────────────────────────────────────────
er = pd.read_csv(EREGULON, sep="\t")
summ = (er.groupby(["eRegulon_name", "TF", "is_extended"])
          .agg(n_genes=("Gene", "nunique"), n_regions=("Region", "nunique"),
               n_triplets=("Gene", "size"))
          .reset_index().sort_values("n_genes", ascending=False))
summ.to_csv(SUMMARY_CSV, index=False)
print(f"[QC] {summ.shape[0]} eRegulons ({er['TF'].nunique()} unique TFs) -> {SUMMARY_CSV}")
print(f"[QC] target genes per eRegulon: median={int(summ.n_genes.median())} "
      f"range={summ.n_genes.min()}-{summ.n_genes.max()}")

# ─────────────────────────────────────────────
# 2. Mouse positive-control markers
# ─────────────────────────────────────────────
er_tfs = set(er["TF"].unique())
hits = [f"{tf}{'✓' if tf in er_tfs else '✗'}" for tf in KNOWN_TFS]
print("[QC] positive-control marker TFs as eRegulons: " + " ".join(hits))
for tf in ("Fezf2", "Bcl11b"):
    print(f"[QC] HEADLINE — {tf} eRegulon present: {'YES' if tf in er_tfs else 'no'} "
          f"(canonical L5 deep-layer identity TF)")

# ─────────────────────────────────────────────
# AUCell matrix + metadata (age_label / roi), keyed by barcode
# ─────────────────────────────────────────────
m = mudata.read(AUCELL)
gauc = m.mod["Gene_based"]                        # cells × eRegulons
auc = pd.DataFrame(gauc.X, index=gauc.obs_names, columns=gauc.var_names)
auc.index = [c[: -len(BC_SUFFIX)] if c.endswith(BC_SUFFIX) else c for c in auc.index]

obs = sc.read_h5ad(L5_GEX_H5AD, backed="r").obs
meta = obs.loc[:, [STAGE_COL, REGION_COL]].reindex(auc.index)
if meta.isna().any().any():
    n_bad = int(meta.isna().any(axis=1).sum())
    raise ValueError(f"{n_bad} AUCell cells not matched to obs metadata — barcode mismatch")
ages_present = [s for s in AGE_ORDER if s in set(meta[STAGE_COL].astype(str))]
meta[STAGE_COL] = pd.Categorical(meta[STAGE_COL].astype(str), categories=ages_present, ordered=True)
print(f"[QC] AUCell cells: {auc.shape[0]:,}; eRegulons: {auc.shape[1]}")
print("[QC] cells per age_label:\n" + meta[STAGE_COL].value_counts().reindex(ages_present).to_string())

# ─────────────────────────────────────────────
# 3a. mean activity by age_label (developmental trajectory) + per-age N
# ─────────────────────────────────────────────
age_means = auc.groupby(meta[STAGE_COL], observed=True).mean().reindex(ages_present)   # age × eRegulon
n_per_age = meta[STAGE_COL].value_counts().reindex(ages_present).rename("n_cells")
age_out = age_means.T.copy()                       # eRegulon × age
age_out.to_csv(STAGE_CSV)
print(f"[QC] wrote per-age activity -> {STAGE_CSV}")

# 3b. mean activity by roi (region covariate)
roi_means = auc.groupby(meta[REGION_COL].astype(str), observed=True).mean()
roi_means.T.to_csv(REGION_CSV)
print(f"[QC] wrote per-roi activity -> {REGION_CSV}")

# ─────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────
# (a) top eRegulons by target-gene count (red = mouse marker positive control)
top = summ.head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 7))
colors = ["#c0392b" if tf in KNOWN_TFS else "#4472c4" for tf in top["TF"]]
ax.barh(top["eRegulon_name"], top["n_genes"], color=colors)
ax.set_xlabel("# target genes")
ax.set_title("Top 20 eRegulons (red = mouse L5-marker positive control)")
plt.tight_layout(); plt.savefig(FIG_SIZES, dpi=130); plt.close()
print(f"[QC] wrote {FIG_SIZES}")

# (b) developmental cascade: eRegulon activity z-scored across age_label, ordered by peak age.
#     Restrict columns to ages with a reasonable N so near-empty embryonic bins don't dominate.
z = (age_means - age_means.mean(axis=0)) / (age_means.std(axis=0) + 1e-9)   # age × eRegulon
z = z.T                                                                       # eRegulon × age
z = z.iloc[np.argsort(np.nan_to_num(z.values).argmax(axis=1))]               # order by peak age
fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(ages_present)), 12))
im = ax.imshow(np.nan_to_num(z.values), aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax.set_xticks(range(len(ages_present)))
ax.set_xticklabels([f"{a}\n(n={int(n_per_age[a])})" for a in ages_present], rotation=45, ha="right", fontsize=6)
ax.set_yticks(range(z.shape[0])); ax.set_yticklabels(z.index, fontsize=5)
ax.set_title("eRegulon developmental cascade\n(z-scored across age_label, AUCell gene-based)")
fig.colorbar(im, ax=ax, shrink=0.4, label="z")
plt.tight_layout(); plt.savefig(FIG_CASCADE, dpi=130); plt.close()
print(f"[QC] wrote {FIG_CASCADE}")

print("[QC] Stage 6 gao25·l5 complete.")
