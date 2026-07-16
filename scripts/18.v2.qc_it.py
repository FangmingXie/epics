#!/usr/bin/env python
"""Stage 6 (v2) — QC / sanity for the all-IT eGRN (L2/3–L6, PFC+V1, all stages; plan/04).

Reworked from the single-cell-type 18_qc.py into a 3-way Type × Region × Group analysis:
  1. eRegulon summary table (TF, direction, #target genes/regions/triplets) -> CSV.
  2. Layer-specific positive controls: known marker TFs per IT layer that surfaced as eRegulons.
     Headline check — did CUX2 / RORB (absent in the single-cell-type L2/3 run, for lack of contrast)
     now form eRegulons?
  3. Regulon × Type activity (layer identity): mean AUCell per IT layer.
  4. Region effect (PFC vs V1) WITHIN matched (Type, Group) strata — the paper's central theme;
     controls for layer + developmental stage given the strong stage×region imbalance.
  5. Developmental trajectories faceted by Type × Region (all 5 stages on x-axis).
  Per-stratum N is reported alongside every activity estimate so effects aren't confounded by
  composition.

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/18.v2.qc_it.py
"""
import os
import itertools
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
WORK = f"{DATA}/work_it"                          # NEW isolated work dir
QC = f"{WORK}/qc"

EREGULON = f"{WORK}/eRegulon.tsv"                 # Stage 5c v2
AUCELL = f"{WORK}/AUCell.h5mu"                    # Stage 5d v2
IT_GEX_H5AD = f"{DATA}/it_gex.h5ad"              # Stage 1 v2 (has Type / Region / Group)

SUMMARY_CSV = f"{QC}/eRegulon_summary.csv"
TYPE_CSV = f"{QC}/eRegulon_activity_by_type.csv"
STAGE_CSV = f"{QC}/eRegulon_activity_by_type_region_stage.csv"
REGION_DIFF_CSV = f"{QC}/eRegulon_region_diff_PFCvsV1.csv"
FIG_SIZES = f"{QC}/eRegulon_target_sizes.png"
FIG_TYPE = f"{QC}/eRegulon_activity_by_type.png"
FIG_REGION = f"{QC}/eRegulon_region_diff_PFCvsV1.png"
FIG_TRAJ = f"{QC}/eRegulon_trajectories_by_type_region.png"

TYPE_COL, REGION_COL, STAGE_COL = "Type", "Region", "Group"
IT_TYPES = ["EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L6-IT"]
REGIONS = ["PFC", "V1"]
STAGE_ORDER = ["First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence"]
BC_SUFFIX = "___it"                              # AUCell cell names carry the cisTopic project suffix

# known marker TFs per IT layer (positive controls). CUX2/RORB are the headline recovery check.
KNOWN_TFS_BY_LAYER = {
    "EN-L2_3-IT": ["CUX2", "POU3F2", "POU3F3", "SATB2", "MEF2C"],
    "EN-L4-IT":   ["RORB", "MEF2C", "RFX3"],
    "EN-L5-IT":   ["FEZF2", "BCL11B", "ETV1", "FOXP2", "SOX5"],
    "EN-L6-IT":   ["TBR1", "FOXP2", "TLE4", "SOX5", "FOXP1", "NR4A2"],
}
ALL_KNOWN_TFS = sorted({t for v in KNOWN_TFS_BY_LAYER.values() for t in v})

for f in (EREGULON, AUCELL, IT_GEX_H5AD):
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
# 2. Layer-specific positive controls
# ─────────────────────────────────────────────
er_tfs = set(er["TF"].unique())
print("[QC] positive-control marker TFs as eRegulons (per layer):")
for layer, tfs in KNOWN_TFS_BY_LAYER.items():
    hits = [f"{tf}{'✓' if tf in er_tfs else '✗'}" for tf in tfs]
    print(f"      {layer}: {' '.join(hits)}")
for tf in ("CUX2", "RORB"):
    print(f"[QC] HEADLINE — {tf} eRegulon present: {'YES' if tf in er_tfs else 'no'} "
          f"(absent in the single-cell-type L2/3 run)")

# ─────────────────────────────────────────────
# AUCell matrix + metadata (Type / Region / Group), keyed by barcode
# ─────────────────────────────────────────────
m = mudata.read(AUCELL)
gauc = m.mod["Gene_based"]                        # cells × eRegulons
auc = pd.DataFrame(gauc.X, index=gauc.obs_names, columns=gauc.var_names)
auc.index = [c[: -len(BC_SUFFIX)] if c.endswith(BC_SUFFIX) else c for c in auc.index]

obs = sc.read_h5ad(IT_GEX_H5AD, backed="r").obs
meta = obs.loc[:, [TYPE_COL, REGION_COL, STAGE_COL]].reindex(auc.index)
if meta.isna().any().any():
    n_bad = int(meta.isna().any(axis=1).sum())
    raise ValueError(f"{n_bad} AUCell cells not matched to obs metadata — barcode mismatch")
meta[TYPE_COL] = pd.Categorical(meta[TYPE_COL], categories=IT_TYPES, ordered=True)
meta[REGION_COL] = pd.Categorical(meta[REGION_COL], categories=REGIONS, ordered=True)
stages_present = [s for s in STAGE_ORDER if s in set(meta[STAGE_COL])]
meta[STAGE_COL] = pd.Categorical(meta[STAGE_COL], categories=stages_present, ordered=True)
print(f"[QC] AUCell cells: {auc.shape[0]:,}; eRegulons: {auc.shape[1]}")
print("[QC] cells per Type×Region:\n" +
      pd.crosstab(meta[TYPE_COL], meta[REGION_COL]).to_string())

# ─────────────────────────────────────────────
# 3. Regulon × Type activity (layer identity)
# ─────────────────────────────────────────────
type_means = auc.groupby(meta[TYPE_COL], observed=True).mean().reindex(IT_TYPES)   # Type × eRegulon
type_means.T.to_csv(TYPE_CSV)
print(f"[QC] wrote per-Type activity -> {TYPE_CSV}")

# ─────────────────────────────────────────────
# 4. Region effect (PFC vs V1) within matched (Type, Group) strata
# ─────────────────────────────────────────────
rows = []
for t, g in itertools.product(IT_TYPES, stages_present):
    sel = (meta[TYPE_COL] == t) & (meta[STAGE_COL] == g)
    sub_pfc = auc[sel & (meta[REGION_COL] == "PFC").values]
    sub_v1 = auc[sel & (meta[REGION_COL] == "V1").values]
    n_pfc, n_v1 = len(sub_pfc), len(sub_v1)
    if n_pfc < 20 or n_v1 < 20:        # skip strata too small to compare reliably
        continue
    diff = sub_v1.mean() - sub_pfc.mean()          # >0 → V1-biased, <0 → PFC-biased
    for reg, d in diff.items():
        rows.append({"Type": t, "Group": g, "eRegulon_name": reg,
                     "n_PFC": n_pfc, "n_V1": n_v1,
                     "mean_PFC": sub_pfc[reg].mean(), "mean_V1": sub_v1[reg].mean(),
                     "V1_minus_PFC": d})
region_diff = pd.DataFrame(rows)
if region_diff.empty:
    print("[QC] WARNING: no (Type,Group) stratum had ≥20 cells in both regions — region diff skipped")
else:
    region_diff.to_csv(REGION_DIFF_CSV, index=False)
    # average the region effect across strata (equal weight per stratum) → consensus PFC/V1 bias
    consensus = (region_diff.groupby("eRegulon_name")["V1_minus_PFC"].mean()
                 .sort_values())
    print(f"[QC] wrote region diff (PFC vs V1, matched Type×Group) -> {REGION_DIFF_CSV} "
          f"({region_diff[['Type','Group']].drop_duplicates().shape[0]} strata)")
    print("[QC] most PFC-biased eRegulons (mean V1−PFC, top 5):")
    for reg, v in consensus.head(5).items():
        print(f"      {reg}: {v:+.4f}")
    print("[QC] most V1-biased eRegulons (mean V1−PFC, top 5):")
    for reg, v in consensus.tail(5).iloc[::-1].items():
        print(f"      {reg}: {v:+.4f}")

# ─────────────────────────────────────────────
# 5. Full Type × Region × Group mean-activity table
# ─────────────────────────────────────────────
grp = meta.groupby([TYPE_COL, REGION_COL, STAGE_COL], observed=True)
stage_tbl = auc.groupby([meta[TYPE_COL], meta[REGION_COL], meta[STAGE_COL]], observed=True).mean()
n_per = grp.size().rename("n_cells")
stage_tbl = stage_tbl.join(n_per)
stage_tbl.to_csv(STAGE_CSV)
print(f"[QC] wrote Type×Region×Group activity (+ per-stratum N) -> {STAGE_CSV}")

# ─────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────
# (a) top eRegulons by target-gene count
top = summ.head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 7))
colors = ["#c0392b" if tf in ALL_KNOWN_TFS else "#4472c4" for tf in top["TF"]]
ax.barh(top["eRegulon_name"], top["n_genes"], color=colors)
ax.set_xlabel("# target genes"); ax.set_title("Top 20 eRegulons (red = layer-marker positive control)")
plt.tight_layout(); plt.savefig(FIG_SIZES, dpi=130); plt.close()
print(f"[QC] wrote {FIG_SIZES}")

# (b) regulon × Type heatmap (z-scored per eRegulon), ordered by peak layer
z = (type_means - type_means.mean(axis=0)) / (type_means.std(axis=0) + 1e-9)   # z per eRegulon
z = z.T                                                                          # eRegulon × Type
z = z.iloc[np.argsort(z.values.argmax(axis=1))]
fig, ax = plt.subplots(figsize=(5, 12))
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax.set_xticks(range(len(IT_TYPES))); ax.set_xticklabels(IT_TYPES, rotation=45, ha="right")
ax.set_yticks(range(z.shape[0])); ax.set_yticklabels(z.index, fontsize=6)
ax.set_title("eRegulon activity by IT layer\n(z-scored per eRegulon, AUCell gene-based)")
fig.colorbar(im, ax=ax, shrink=0.4, label="z")
plt.tight_layout(); plt.savefig(FIG_TYPE, dpi=130); plt.close()
print(f"[QC] wrote {FIG_TYPE}")

# (c) region-diff bar (consensus V1−PFC across matched strata) — top |effect| eRegulons
if not region_diff.empty:
    consensus = region_diff.groupby("eRegulon_name")["V1_minus_PFC"].mean()
    topn = consensus.reindex(consensus.abs().sort_values(ascending=False).index).head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.barh(topn.index, topn.values,
            color=["#2c7fb8" if v > 0 else "#d95f0e" for v in topn.values])
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("mean(V1 − PFC) activity  (blue = V1-biased, orange = PFC-biased)")
    ax.set_title("Region-differential eRegulon activity\n(matched Type × Group strata)")
    ax.tick_params(axis="y", labelsize=6)
    plt.tight_layout(); plt.savefig(FIG_REGION, dpi=130); plt.close()
    print(f"[QC] wrote {FIG_REGION}")

# (d) developmental trajectories of the most-dynamic eRegulons, faceted by Type × Region
overall_stage = auc.groupby(meta[STAGE_COL], observed=True).mean().reindex(stages_present)
dyn = (overall_stage.max(axis=0) - overall_stage.min(axis=0)).sort_values(ascending=False)
top_dyn = list(dyn.head(6).index)
fig, axes = plt.subplots(len(IT_TYPES), len(REGIONS),
                         figsize=(4 * len(REGIONS), 2.4 * len(IT_TYPES)),
                         sharex=True, sharey=True, squeeze=False)
x = range(len(stages_present))
for i, t in enumerate(IT_TYPES):
    for j, rgn in enumerate(REGIONS):
        ax = axes[i][j]
        sel = (meta[TYPE_COL] == t) & (meta[REGION_COL] == rgn)
        sub = auc[sel.values]
        sub_stage = meta.loc[sel.values, STAGE_COL]
        mns = sub.groupby(sub_stage, observed=True).mean().reindex(stages_present)
        ncell = sub_stage.value_counts().reindex(stages_present).fillna(0).astype(int)
        for reg in top_dyn:
            ax.plot(x, mns[reg].values, marker="o", ms=3, lw=1, label=reg)
        ax.set_title(f"{t} / {rgn}  (n={int(ncell.sum())})", fontsize=7)
        if i == len(IT_TYPES) - 1:
            ax.set_xticks(list(x)); ax.set_xticklabels(stages_present, rotation=45, ha="right", fontsize=6)
axes[0][-1].legend(fontsize=5, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.suptitle("Top-6 dynamic eRegulons across development, by IT layer × region", y=1.0)
plt.tight_layout(); plt.savefig(FIG_TRAJ, dpi=130, bbox_inches="tight"); plt.close()
print(f"[QC] wrote {FIG_TRAJ}")

print("[QC] Stage 6 v2 complete.")
