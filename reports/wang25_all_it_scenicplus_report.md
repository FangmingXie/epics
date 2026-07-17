# Enhancer-driven GRN (SCENIC+) for all IT neurons (L2/3–L6), PFC + V1 — Wang et al. 2025

**Date:** 2026-07-17 · **Env:** `epics` (SCENIC+ 1.0a2, pycisTopic 2.0a0) + `chrombpnet` (htslib/bedtools) + `macs3`
**Plan:** [`plan/04.wang25_all_it_scenicplus_executable.md`](../plan/04.wang25_all_it_scenicplus_executable.md) · **Scripts:** `scripts/*.v2.*` · **Driver:** `scripts/run_all_it.v2.sh`

> Companion to the single-region run [`reports/wang25_l23_scenicplus_report.md`](wang25_l23_scenicplus_report.md)
> (V1 `EN-L2_3-IT`, 5,558 cells). This run scales the same pipeline to **all four layer-resolved IT excitatory
> types across both cortical regions and all developmental stages**. Outputs are isolated under
> `data/wang25/work_it/` + DB prefix `wang25_it`; the L2/3 run is untouched.

---

## 1. Summary

We inferred an enhancer-driven gene-regulatory network (eGRN / eRegulons) for **all intratelencephalic (IT)
excitatory neurons spanning cortical layers L2/3–L6**, in **both prefrontal (PFC) and primary visual (V1)
cortex**, across **all five developmental stages**, from the paired 10x snMultiome atlas of
**Wang L. et al., *Nature* 647:169–178 (2025)**, "Molecular and cellular dynamics of the developing human
neocortex" (hg38).

**Headline result:** a **147-eRegulon network (113 TFs)** from direct cistromes (+ a 93-eRegulon extended
variant). Adding cross-layer contrast delivered its central expected payoff — **CUX2 now forms an eRegulon**
(it could not in the single-cell-type L2/3 run) — and the network resolves **layer-identity** programs
(CUX2/L2-3, CUX1·RORA/L4, BCL11A/L5, PBX3·NFIA/L6) and **areal (PFC↔V1) regulatory bias**, the paper's central
theme, which the single-region run could not address at all.

---

## 2. Data & scope

| | |
| --- | --- |
| Assay | Paired 10x snMultiome (RNA + ATAC, same nuclei), hg38 |
| Source | CELLxGENE (RNA h5ad + ATAC fragments) — shared with the L2/3 run |
| Cell types | `Type ∈ {EN-L2_3-IT, EN-L4-IT, EN-L5-IT, EN-L6-IT}` (layer-resolved IT only) |
| Regions | `Region ∈ {PFC, V1}` (the tiny "General" fetal pool and unlayered `EN-IT-Immature` excluded) |
| Stages | all 5 `Group` levels: First / Second / Third_trimester, Infancy, Adolescence |
| Cells | **49,283** (~9× the L2/3 run) |

Population sizes (`obs`):

| Type | PFC | V1 | total |
| --- | ---: | ---: | ---: |
| EN-L2_3-IT | 6,292 | 5,558 | 11,850 |
| EN-L4-IT | 8,239 | 11,834 | 20,073 |
| EN-L5-IT | 6,069 | 3,552 | 9,621 |
| EN-L6-IT | 5,335 | 2,404 | 7,739 |
| **total** | **25,935** | **23,348** | **49,283** |

**Scope rationale.** The L2/3 run deliberately restricted to V1 to *avoid* the PFC↔V1 confound. Here we do the
opposite on purpose: pool layers and regions so the model has the **cell-type and areal contrast** needed to
build layer-identity and region-biased regulons. Strong stage×region imbalance persists (V1 postnatal-dominated;
PFC carries more prenatal cells) — carried as a caveat into all Stage-6 stratified comparisons.

---

## 3. Pipeline (Stages 1–7)

Design choices for this run (user-confirmed): **pooled single MACS3 peak call** (not per-population
pseudobulk), **LDA topics only** for region sets (no DAR / `find_diff_features`), and a **reduced topic sweep
`[20,30,40]`**. Stages 0/2/2b (download, refs, UCSC-fragment normalization) are population-independent and were
reused from the L2/3 run.

| Stage | Script (`.v2.`) | What | Key result | Wall |
| --- | --- | --- | --- | --- |
| 1 | `10.v2.select_it.py` | Subset 4 IT types ∩ PFC/V1; raw counts → `.X`; gene symbols as var_names | **49,283 × 35,468 genes** | — |
| 3a | `12.v2.call_peaks_it.sh` + `12b.v2…` | MACS3 `-g hs` on **1.01 B** pooled IT fragments → consensus | 306,132 raw → **196,132** non-overlapping (501 bp) | ~1.6 h |
| 3b | `13.v2.build_ctx_db_it.sh` | Custom cisTarget motif DB on the consensus peaks | **196,132 regions × 10,249 motifs** (rankings 5.4 GB, scores 3.6 GB) | ~1.7 h |
| 3c | `14a.v2…` + `14.v2.run_pycistopic_it.py` | IT-only fragments (14 GB) → cisTopic + LDA sweep | 196,132 × 49,283; **40 topics** (auto-selected) → 40 region-sets | ~8.7 h |
| 4 | `16.v2.prepare_menr_it.sh` | MuData, search space, motif enrichment, cistromes | MuData 49,283 × (35,468 + 196,132); search space 593,670; **833 cistrome TFs** | ~7 min |
| 5 | `17.v2.run_grn_it.sh` (+ `17b.v2…`) | TF→gene, region→gene, eGRN, AUCell (direct + extended) | **147 eRegulons / 113 TFs / 36,372 triplets** | ~4.1 h |
| 6 | `18.v2.qc_it.py` | 3-way Type×Region×Group QC + per-layer controls | layer-identity + PFC↔V1 + developmental | <1 min |
| 7 | `19.v2.organize_regulons_it.py` | Merge direct+extended → clean regulon→gene table | **104 regulons (88 direct + 16 extended), 8,889 (regulon,gene) rows** | <1 min |

Total wall-clock ~16 h (2026-07-16 15:02 → 2026-07-17 07:16), single run, no failures. The **LDA sweep
dominated** (~7.6 h; the 40-topic model alone ~7.2 h single-threaded CGS at 49k×196k). Everything else, incl.
Stage-5 GRN, was faster than the pre-run estimates (TF2G ~1.5 h, R2G ~2 h, eGRN+AUCell ~30 min).

Scale vs the L2/3 run: consensus peaks **141,169 → 196,132**; cistrome TFs **756 → 833**; direct eRegulons
**39 → 147**; triplets **9,766 → 36,372** — the extra layers/regions expose substantially more regulatory
structure, as intended.

---

## 4. Key results

### 4.1 eGRN — direct vs extended cistromes

| | Direct (`eRegulon.tsv`) | Extended (`eRegulon_extended.tsv`) |
| --- | --- | --- |
| eRegulons | **147** | 93 |
| unique TFs | 113 | 71 |
| triplets | 36,372 | 18,438 |
| target genes / eRegulon | median 58 (10–1,782) | median 48 (10–1,720) |

As in the L2/3 run, the two are **complementary, not nested** (orthology annotation reassigns region↔TF links).
The largest direct eRegulons are broad/ubiquitous factors (KLF9, ZNF148, EGR1, ETV5, ATF4, YY1, NRF1, MAZ) —
the biologically informative signal is in the mid-size, layer- and region-selective regulons below.

### 4.2 Positive controls — CUX2 recovered

| TF | Direct | Extended | Expected layer | Note |
| --- | --- | --- | --- | --- |
| **CUX2** | ✅ | – | L2/3–L4 | **absent in the L2/3-only run — now recovered** (the key benefit of cross-layer contrast) |
| POU3F2 (BRN2) | ✅ | ✅ | L2/3 | upper-layer |
| POU3F3 (BRN1) | ✅ | ✅ | L2/3 | upper-layer |
| MEF2C | ✅ | ✅ | L2/3–L4 | |
| RFX3 | ✅ | ✅ | L2/3–L4 | |
| SATB2 | – | ✅ | pan-IT | direct→extended only (same as L2/3 run) |
| TBR1 | ✅ | ✅ | L6 | deep-layer |
| RORB | ✗ | ✗ | L4 | **no cistrome** even with L4 present (motif not enriched) |
| FEZF2, BCL11B, ETV1, FOXP2, SOX5, TLE4, FOXP1, NR4A2 | ✗ | ✗ | L5/L6 | deep-layer markers did not form eRegulons |

So cross-layer contrast recovered **CUX2** (upper-layer) and **TBR1** (deep-layer) that the single-type run
lacked, but the **L5/L6 identity TFs (FEZF2, BCL11B, FOXP2, SOX5, TLE4) still fail to form eRegulons**, and
**RORB** remains absent — the main biological gaps (see §7).

### 4.3 Layer identity (AUCell by `Type`)

Each direct eRegulon's peak layer (z-scored activity across the four types); counts: **L2/3 45, L4 43, L5 33,
L6 26**. Representative layer-peaking regulons (by activity range) are biologically coherent:

```
L2/3  CUX2, ZNF704, SMAD3, GLIS1, MEIS1, RFX3        (upper-layer identity — CUX2 peaks in its own layer)
L4    CUX1, RORA, JDP2, TEF, NFIX, MEF2D             (L4/upper — CUX1 & RORA mark L4)
L5    BCL11A, DRAP1, CREB3, TEAD1, NR1D1, CHURC1
L6    PBX3, NFIA, TCF4, KLF6, NR3C1, SMARCC1          (deep-layer)
```

Figure: `work_it/qc/eRegulon_activity_by_type.png` (z-scored regulon × layer heatmap).

### 4.4 Areal regulation — PFC vs V1 (the paper's theme)

Differential AUCell computed **within matched (Type, Group) strata** (16 strata with ≥20 cells per region),
so region effects are not confounded by layer or developmental composition:

| V1-biased (mean V1−PFC > 0) | PFC-biased (mean V1−PFC < 0) |
| --- | --- |
| **NFIX, POU3F2, PBX3, RORA, CUX1** | **RARB, PDLIM5, KLF6, TCF4, NR3C1** |

The V1-biased set includes recognized areal/upper-layer factors (POU3F2, PBX3, NFIX); the PFC-biased set
includes RARB (retinoic-acid signalling, a known anterior-cortex patterning input) and the deep-layer/mature
TCF4·NR3C1·KLF6 module. Full per-stratum table: `work_it/qc/eRegulon_region_diff_PFCvsV1.csv`; figure:
`eRegulon_region_diff_PFCvsV1.png`.

### 4.5 Developmental dynamics (Type × Region × Group)

Per-cell activity means for every eRegulon across all five stages, faceted by layer × region, are written to
`work_it/qc/eRegulon_activity_by_type_region_stage.csv` (with per-stratum N) and plotted in
`eRegulon_trajectories_by_type_region.png` (top-6 most dynamic regulons per facet). *Caveat:* prenatal V1
strata are very low-N (e.g. V1 3rd-trimester 935 cells spread across layers; some Type×Region×stage cells <100),
so early-stage trajectories in V1 are noisy — read alongside the reported N.

---

## 5. Notes / design differences vs the L2/3 run

- **Pooled peak calling** (one MACS3 call over all 49,283 IT barcodes) rather than per-population pseudobulk.
  Simpler and faithful to the L2/3 recipe; a per-layer×region pseudobulk call (SCENIC+-canonical) is the main
  lever if L5/L6-specific enhancers prove under-recovered (see §7).
- **Topics-only region sets** (no DAR). With four layers now present, adding `find_diff_features` DARs by Type
  and by Region is the natural next enrichment — it would supply cell-type-contrastive region sets that could
  recover identity-TF cistromes (RORB, FEZF2, …) the topic-only sets miss.
- **Reduced LDA sweep `[20,30,40]`**; auto-selected **40** (top of the range → coherence still rising). As in
  the L2/3 run this is accepted (topics are only candidate region-sets, re-filtered by the eGRN), but a
  higher-topic sweep is a cheap refinement.
- All the L2/3 runtime fixes carried over unchanged (UCSC/primary chroms, `.gz` fragment naming, `n_cpu=1`
  cisTopic serial join, gene-symbol var_names, `search_spance` chromsizes header, `___it` barcode suffix,
  `prepare_menr` cistromes). No new gotchas surfaced at 9× scale.

---

## 6. Reproducibility

- **Scripts:** `scripts/10.v2.select_it.py` … `scripts/19.v2.organize_regulons_it.py` (+ driver
  `run_all_it.v2.sh`); all paths defined at the top of each script per repo convention. The finished L2/3
  scripts are untouched.
- **Envs:** `epics` (SCENIC+/pycisTopic/pycistarget), `chrombpnet` (htslib/bedtools), `macs3`. Nothing
  installed into `epics`.
- **Commit:** `60fe757` ("Add v2 all-IT SCENIC+ pipeline …") on `main`.
- **Per-stage logs:** `logs/*.v2.log` (+ driver `logs/run_all_it.v2.log`).

### Output manifest (`data/wang25/`, gitignored)

| File | Contents |
| --- | --- |
| `work_it/cistopic_obj_model.pkl` | cisTopic object + 40-topic LDA (4.2 GB) |
| `db/wang25_it.regions_vs_motifs.{rankings,scores}.feather` | custom cisTarget DB (5.4 / 3.6 GB) |
| `work_it/ACC_GEX.h5mu` | multiome MuData 49,283 × (35,468 + 196,132) |
| `work_it/cistromes_{direct,extended}.h5ad` | motif→TF cistromes (833 TFs) |
| `work_it/eRegulon.tsv` / `eRegulon_extended.tsv` | **final eGRNs** (8.0 / ~4 MB) |
| `work_it/AUCell.h5mu` / `AUCell_extended.h5mu` | per-cell eRegulon activity (Gene + Region based; 122 / 81 MB) |
| `work_it/regulon_gene_table.tsv` | merged direct+extended regulon→gene table (104 regulons, 8,889 rows) |
| `work_it/qc/` | summary CSV, per-Type & Type×Region×Group activity CSVs, PFC-vs-V1 diff CSV, 4 figures |

---

## 7. Suggested next steps

1. **Add DAR region sets** (`find_diff_features` by Type and by Region) alongside topics, then re-run motif
   enrichment — the most direct route to recovering the missing layer-identity cistromes (RORB, FEZF2, BCL11B,
   FOXP2, SOX5). This is the biggest methodological lever left on the table.
2. **Per-population pseudobulk peak calling** (Type×Region, then `get_consensus_peaks` over the 8 sets) if (1)
   suggests deep-layer enhancers are under-represented in the pooled peak set.
3. **Higher-topic LDA sweep** (e.g. `[40,60,80]`) since coherence was still rising at 40.
4. **Formalize the areal analysis**: per-regulon linear model `AUCell ~ Type + Region + Group + Type:Region`
   to partition variance and rank region-identity regulons with significance, beyond the current matched-stratum
   means.
5. **Cross-check layer/areal regulons** against the paper's reported cell-type and PFC↔V1 signatures (Wang25),
   and against the V1-L2/3 single-region eGRN for consistency on the shared L2/3 program.
