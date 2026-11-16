# Reproducing the Results

This document describes how to reproduce the figures and tables reported in:

> **"Improving Target–Adverse Event Association Prediction via Inverse Node Frequency
> Transformation to Mitigate Knowledge Graph Topological Bias"**

File, test and notebook-section names follow the manuscript numbering (`figN_*` =
Figure N; `figS*` / `tableS*` = Supplementary Information).

---

## What is already included in this package

| Location | Contents |
|---|---|
| `data/fig4_top20_targets.csv` | Pre-computed data for Figure 4 (top 20 targets by AE count) |
| `data/hetionet-v1.0-nodes_version2.tsv` | Hetionet v1.0 node mapping (Entrez ID → gene name) |
| `data/soldatos_2018/` | Soldatos et al. (2018) supplementary files 2, 3, and 4 |
| `data/pickles/` | Small pre-computed pickles (CV results, slice metrics, feature importance) |
| `data/figure_data/` | Plotted values of Figures 8, 11 and S4, whose inputs are restricted |
| `metapaths_src/` | The `metapaths` package (INF toolbox + utilities) |
| `requirements.txt` | Pinned Python dependencies |

The large pickles (~9.5 GB: trained models, prediction DataFrames, full feature matrices)
contain the proprietary Evotec target–AE labels and are available on request — see the
next section. They are only needed to recompute Figures 8, 11 and S4 from the predictions.

---

## Step 1 — (Optional) obtain the restricted pickle files

The large prediction pickles contain the Evotec ground truth and are not
distributed publicly; they are available from the corresponding author on request.
Place the files in `data/pickles/` as described in `data/pickles/DOWNLOAD.md`.

Without them: Figures 4, 6, 7, 9 and 10 and Supplementary Tables S1–S2 / Figure S3
are reproduced from the bundled pickles; Figure 8 is re-plotted, and the values of
Figures 8, 11 and S4 are provided, from `data/figure_data/`. With them,
`supplement.ipynb` recomputes Figures 7–9, 11 and S4 from the predictions.

---

## Step 2 — Set up the Python environment

```bash
python3.10 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install ./metapaths_src --use-pep517
```

No path configuration is required — `config.py` uses paths relative to the project
root by default.

---

## Step 3 — Run the reproducibility notebook

`supplement.ipynb` is the canonical notebook for reproducing the paper outputs
from the underlying predictions (requires the restricted files from Step 1).
It contains only the cells needed to produce Figures 7–9 and 11 and Supplementary
Figure S4, with all
intermediate results loaded from pre-computed pickles.

```bash
.venv/bin/jupyter nbconvert --to notebook --execute --inplace supplement.ipynb
```

Expected runtime: ~20 minutes (dominated by the permutation-importance cell).

| Output | Notebook section | Requires restricted files? |
|---|---|---|
| Figure 11 | `### Figure 11 — Permutation importance on Soldatos external test set` | Yes |
| Figure S4 | `## Figure S4 — Synthesis: Holdout vs. External Test Performance` | Yes |
| Figures 7, 8, 9 | `## Figures 7, 8, 9 — Per-decile performance` | Yes |

> `INF_notebook_clean.ipynb` is the full exploratory notebook. Its "Top features"
> section produces the holdout permutation-importance results shown in Figure 10
> (values bundled in `data/pickles/config_top_feats_by_slice.pkl`); everything
> else reported is reproduced by `supplement.ipynb`.

---

## Step 4 — Verify with automated tests

```bash
source .venv/bin/activate
pytest tests/ -v
```

All 36 tests should pass using only the bundled data (no restricted files needed):

| Test file | Paper claim verified |
|---|---|
| `tests/test_fig4.py` | Figure 4: top target is CXCL8 (333 AEs); top-20 sum = 2,997 (~45%) |
| `tests/test_fig6.py` | Figure 6: CV TPR at deciles 1, 2, 10 for L3 and L4 (tolerance ±0.001) |
| `tests/test_fig10.py` | Figure 10: rank-1 INF feature is Gene_Compound_Gene_AE; full top-5 rankings for INF, raw counts and DWPC |
| `tests/test_figure_data.py` | Figures 8, 11, S4: bundled values match the manuscript (holdout 1st-decile TPR 0.48 / 0.26 / 0.20; overall TPR 0.90 / 0.87 / 0.85 holdout, 0.40 / 0.22 / 0.22 external; rank-1 external INF feature Gene_Compound_Gene_AE) |
| `tests/test_si_inf_grid.py` | §3.2 / Tables S1–S2: decile-1 best config 12 (TPR 0.424); top-3 at decile 1 use deflator exponent 0; best global config 20 |

---

## Expected outputs

### Figure 4 — Top 20 targets by AE count
- **Source file:** `data/fig4_top20_targets.csv` (no notebook cell required)
- **Key values:** CXCL8 (Gene::3576) at rank 1 with 333 AEs; top-20 total = 2,997 (≈45% of 6,657)

### Figure 6 — Cross-validation TPR by decile
- **Source:** `figures/fig6_cv_tpr.png` (pre-generated; all values hardcoded)
- **Key values:**

  | Feature set | Decile 1 TPR | Decile 2 TPR | Decile 10 TPR |
  |---|---|---|---|
  | L3 (raw counts) | 0.202 | 0.560 | 0.959 |
  | L4 (raw counts) | 0.223 | 0.575 | 0.966 |

### Figure 8 — Holdout test set, 1st decile
- **Script:** `generate_figs7_8_9.py` (uses `data/figure_data/fig8_holdout_d1.csv` if the restricted predictions are absent)
- **Key values (TPR):** Best INF Config 0.476; DWPC 0.256; Original 0.201. Hits@10: 0.376 / 0.376 / 0.321

### Figure S4 — Holdout vs. external test performance (all targets)
- **Notebook:** `supplement.ipynb` → `## Figure S4 — Synthesis`
- **Output file:** `figures/figS4_holdout_vs_external.png`; values in `data/figure_data/figS4_holdout_vs_external.csv`
- **Key values:**

  | Metric | Test set | INF Ensemble | DWPC | Original |
  |---|---|---|---|---|
  | Hits@10 | Holdout (n=998) | 0.353 | 0.362 | 0.361 |
  | Hits@10 | External (n=7,804) | 0.064 | 0.068 | 0.064 |
  | TPR | Holdout (n=998) | 0.901 | 0.868 | 0.847 |
  | TPR | External (n=7,804) | 0.402 | 0.215 | 0.217 |

  The INF ensemble is retrained and the bootstrap is unseeded, so recomputed values can
  differ in the third decimal (see `data/figure_data/README.md`).

### Figure 10 — Holdout feature importance, 1st decile
- **Source:** `data/pickles/config_top_feats_by_slice.pkl` (bundled), produced by the
  "Top features" section of `INF_notebook_clean.ipynb`
- **Key values (INF model):** Gene–Compound–Gene–AE 0.130; Gene–Pathway–Gene–AE 0.081;
  Gene–MolecularFunction–Gene–AE 0.079; Gene–BiologicalProcess–Gene–AE 0.079; Gene–Gene–AE 0.059

### Figure 11 — External-set feature importance, 1st decile
- **Notebook:** `supplement.ipynb` → `### Figure 11 — Permutation importance on Soldatos external test set`
- **Values:** `data/figure_data/fig11_external_feat_importance_d1.csv`
- **Key result:** Rank-1 feature for the INF model is `Gene_Compound_Gene_AE`
  (drug-target co-occurrence path), consistent with the pharmacovigilance premise

---

## Regenerating pickles from scratch

The pickles can be regenerated end-to-end, but require a running Neo4j instance loaded
with the TSA knowledge graph. Contact the authors for access.

If the graph is available, set `NEO4J_URI_READ` / `NEO4J_URI_WRITE` in `config.py` and
run `INF_notebook_clean.ipynb` from the top. Estimated runtimes: L3 CV grid search ~3 hours;
L4 ~10 hours; holdout test evaluation ~2–4 hours.

Feature extraction (metapath count and INF computation) is not included in this notebook.
Contact the authors for the extraction scripts.

---

## Regenerating Figure 4 data from the graph

`generate_fig4_data.py` re-queries Neo4j for the top-20 targets and writes
`data/fig4_top20_targets.csv`. It includes built-in assertions that verify the
paper's cited values (sum = 2,997, share ≈ 45%):

```bash
source .venv/bin/activate
python generate_fig4_data.py
```

The CSV is already bundled so this step is optional.
