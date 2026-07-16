#!/usr/bin/env python
"""Stage 6 — QC / sanity for the V1 L2/3 eGRN.

1. eRegulon summary table (TF, direction, #target genes/regions/triplets) -> CSV.
2. Positive-control check: known L2/3 IT TFs that surfaced as eRegulons.
3. AUCell developmental dynamics: mean eRegulon activity across the 4 stages V1 L2/3 spans
   (Second_trimester -> Third_trimester -> Infancy -> Adolescence) — the paper's theme.
4. Plots: top eRegulons by target size; z-scored eRegulon activity heatmap across stages.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/18_qc.py
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
DATA = f"{PROJ}/data/wang25"
WORK = f"{DATA}/work"
QC = f"{WORK}/qc"

EREGULON = f"{WORK}/eRegulon.tsv"                 # Stage 5c
AUCELL = f"{WORK}/AUCell.h5mu"                    # Stage 5d
L23_GEX_H5AD = f"{DATA}/v1_l23_gex.h5ad"          # Stage 1 (has Group = developmental stage)

SUMMARY_CSV = f"{QC}/eRegulon_summary.csv"
STAGE_CSV = f"{QC}/eRegulon_activity_by_stage.csv"
FIG_SIZES = f"{QC}/eRegulon_target_sizes.png"
FIG_STAGE = f"{QC}/eRegulon_activity_by_stage.png"

STAGE_COL = "Group"
STAGE_ORDER = ["Second_trimester", "Third_trimester", "Infancy", "Adolescence"]
BC_SUFFIX = "___l23"                             # AUCell cell names carry the cisTopic project suffix
KNOWN_L23_TFS = ["CUX2", "SATB2", "POU3F2", "RORB", "MEF2C", "FOXP1"]

for f in (EREGULON, AUCELL, L23_GEX_H5AD):
    if not os.path.exists(f):
        raise FileNotFoundError(f"required input missing: {f}")
os.makedirs(QC, exist_ok=True)

# ─────────────────────────────────────────────
# 1. eRegulon summary table
# ─────────────────────────────────────────────
er = pd.read_csv(EREGULON, sep="\t")
summ = (er.groupby(["eRegulon_name", "TF", "is_extended"])
          .agg(n_genes=("Gene", "nunique"),
               n_regions=("Region", "nunique"),
               n_triplets=("Gene", "size"))
          .reset_index()
          .sort_values("n_genes", ascending=False))
summ.to_csv(SUMMARY_CSV, index=False)
print(f"[QC] {summ.shape[0]} eRegulons ({er['TF'].nunique()} unique TFs) -> {SUMMARY_CSV}")
print(f"[QC] target genes per eRegulon: median={int(summ.n_genes.median())} "
      f"range={summ.n_genes.min()}-{summ.n_genes.max()}")

# ─────────────────────────────────────────────
# 2. Positive controls
# ─────────────────────────────────────────────
er_tfs = set(er["TF"].unique())
print("[QC] positive-control L2/3 TFs as eRegulons:")
for tf in KNOWN_L23_TFS:
    print(f"      {tf}: {'YES' if tf in er_tfs else 'no'}")

# ─────────────────────────────────────────────
# 3. AUCell developmental dynamics (Gene-based)
# ─────────────────────────────────────────────
m = mudata.read(AUCELL)
gauc = m.mod["Gene_based"]                        # cells × eRegulons
auc = pd.DataFrame(gauc.X, index=gauc.obs_names, columns=gauc.var_names)
auc.index = [c[: -len(BC_SUFFIX)] if c.endswith(BC_SUFFIX) else c for c in auc.index]

obs = sc.read_h5ad(L23_GEX_H5AD, backed="r").obs
stage = obs[STAGE_COL].reindex(auc.index)
if stage.isna().any():
    raise ValueError(f"{stage.isna().sum()} AUCell cells not matched to {STAGE_COL} — barcode mismatch")

stages_present = [s for s in STAGE_ORDER if s in set(stage)]
stage_means = auc.groupby(stage).mean().reindex(stages_present)      # stages × eRegulons
stage_means.T.to_csv(STAGE_CSV)
print(f"[QC] stages present ({stage.value_counts().reindex(stages_present).to_dict()})")
print(f"[QC] wrote per-stage activity -> {STAGE_CSV}")

# ─────────────────────────────────────────────
# 4a. Plot: top eRegulons by target-gene count
# ─────────────────────────────────────────────
top = summ.head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 7))
colors = ["#c0392b" if tf in KNOWN_L23_TFS else "#4472c4" for tf in top["TF"]]
ax.barh(top["eRegulon_name"], top["n_genes"], color=colors)
ax.set_xlabel("# target genes"); ax.set_title("Top 20 eRegulons (red = L2/3 positive control)")
plt.tight_layout(); plt.savefig(FIG_SIZES, dpi=130); plt.close()
print(f"[QC] wrote {FIG_SIZES}")

# ─────────────────────────────────────────────
# 4b. Plot: z-scored eRegulon activity across developmental stages
# ─────────────────────────────────────────────
z = (stage_means - stage_means.mean(axis=0)) / (stage_means.std(axis=0) + 1e-9)   # z per eRegulon
z = z.T                                                                            # eRegulons × stages
z = z.iloc[np.argsort(z.values.argmax(axis=1))]                                    # order by peak stage
fig, ax = plt.subplots(figsize=(5, 11))
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax.set_xticks(range(len(stages_present))); ax.set_xticklabels(stages_present, rotation=45, ha="right")
ax.set_yticks(range(z.shape[0])); ax.set_yticklabels(z.index, fontsize=6)
ax.set_title("eRegulon activity across development\n(z-scored per eRegulon, AUCell gene-based)")
fig.colorbar(im, ax=ax, shrink=0.4, label="z")
plt.tight_layout(); plt.savefig(FIG_STAGE, dpi=130); plt.close()
print(f"[QC] wrote {FIG_STAGE}")

# eRegulons with strongest developmental dynamics (range of stage-mean activity)
dyn = (stage_means.max(axis=0) - stage_means.min(axis=0)).sort_values(ascending=False)
print("[QC] most dynamic eRegulons across development (top 8):")
for name, rng in dyn.head(8).items():
    peak = stage_means[name].idxmax()
    print(f"      {name}: range={rng:.3f}, peaks at {peak}")
print("[QC] Stage 6 complete.")
