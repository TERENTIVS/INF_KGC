# INF — Inverse Node Frequency for biomedical KG completion

Code and data accompanying the manuscript:

> **Improving Target–Adverse Event Association Prediction via Inverse
> Node Frequency Transformation to Mitigate Knowledge Graph
> Topological Bias.**
> Terence Egbelo, J Charles G Jeynes, Zeyneb Kurt, Rahul V Veettil, Valerie J Gillet.
> *Bioinformatics* (under review).

Submitted version archived at Software Heritage:
[`swh:1:rev:09ff3525acb878c37a58d744b32ff11e45b14fad`](https://archive.softwareheritage.org/swh:1:rev:09ff3525acb878c37a58d744b32ff11e45b14fad).

INF is a feature transformation for KG completion that re-weights
metapath instances by the inverse frequency of each intermediate node
under a specific edge type. It mitigates the topological / degree bias
that systematically degrades KG-completion predictions for sparsely
annotated entities — in this work, understudied protein targets for
which AE inference is most needed.

## What's in this repo

| Path | Purpose |
|---|---|
| `supplement.ipynb` | Canonical reproduction notebook for manuscript Figures 7–9 and 11 and Supplementary Figure S4 (requires the restricted data). |
| `generate_figs7_8_9.py` | Standalone script for Figures 7, 8, 9 (CV / holdout / external D1 bars). |
| `generate_si_inf_grid.py` | Generates Supplementary Tables S1–S2 and Figure S3 (per-decile CV TPR across the 35-point INF parameter grid). |
| `generate_fig4_data.py` | Re-queries source for Figure 4 (top-20 targets by AE count). |
| `INF_notebook_clean.ipynb` | Full exploratory notebook; source of the holdout feature-importance results in Figure 10. |
| `metapaths_src/` | The `metapaths` Python package — INF implementation and utilities. |
| `tests/` | Anchored reproduction tests (`pytest tests/` ≈ 10 s). |
| `data/` | Bundled data (Figure 4 CSV, Hetionet node mapping, Soldatos 2018 supplementary files, small pickles, `figure_data/` values for Figures 8, 11, S4). |
| `figures/` | Pre-rendered PNG figures from the paper and SI. |
| `REPRODUCING_RESULTS.md` | Detailed reproduction recipe — per-figure source mapping, runtimes, dependencies. |

## Quickstart

```bash
# 1. Set up the environment (Python 3.10)
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install ./metapaths_src --use-pep517

# 2. Verify every reported value against the bundled data (~10 s)
.venv/bin/pytest tests/ -v

# 3. Regenerate figures from the bundled data
python generate_figs7_8_9.py      # Figures 7, 8, 9
python generate_si_inf_grid.py    # Tables S1-S2, Figure S3

# 4. (Optional, restricted data) Recompute Figures 7-9, 11 and S4 from
#    the underlying predictions (~25 minutes). Requires the restricted
#    pickles, available from the first author on request:
#    see data/pickles/DOWNLOAD.md
.venv/bin/jupyter nbconvert --to notebook --execute --inplace supplement.ipynb
```

`REPRODUCING_RESULTS.md` covers what each step produces, runtimes per
output, and which restricted files are needed for which figure. The
`pytest` suite (36 tests, ~10 s) runs against the bundled data alone
and anchors every numerical claim in the paper to its source pickle.

## Outputs by manuscript item

File, test and notebook-section names follow the manuscript numbering
(`figN_*` = Figure N; `figS*` / `tableS*` = Supplementary Information).

| Manuscript item | Source | Test |
|---|---|---|
| Figure 4 — top 20 targets by AE count | `data/fig4_top20_targets.csv` (`generate_fig4_data.py`) | `tests/test_fig4.py` |
| Figure 6 — CV TPR, deciles 1, 2, 10 | `figures/fig6_cv_tpr.png` | `tests/test_fig6.py` |
| Figures 7, 8, 9 — 1st-decile CV / holdout / external performance | `generate_figs7_8_9.py`; `supplement.ipynb`; Figure 8 values in `data/figure_data/` | `tests/test_figure_data.py` (Figure 8) |
| Figure 10 — holdout feature importance, 1st decile | `INF_notebook_clean.ipynb` ("Top features"); values in `data/pickles/config_top_feats_by_slice.pkl` | `tests/test_fig10.py` |
| Figure 11 — external-set feature importance, 1st decile | `supplement.ipynb` ("Figure 11" section); values in `data/figure_data/` | `tests/test_figure_data.py` |
| Table S1 — CV TPR, all 35 INF configurations | `figures/tableS1_inf_grid_heatmap.png` (`generate_si_inf_grid.py`) | `tests/test_si_inf_grid.py` |
| Table S2 — top-3 INF configurations per decile | `data/tableS2_inf_grid_top3.{csv,md}` (`generate_si_inf_grid.py`) | `tests/test_si_inf_grid.py` |
| Figure S3 — decile-best vs global-best vs raw | `figures/figS3_decile_best_vs_global_vs_raw.png` (`generate_si_inf_grid.py`) | — |
| Figure S4 — holdout vs. external performance, all targets | `figures/figS4_holdout_vs_external.png` (`supplement.ipynb`); values in `data/figure_data/` | `tests/test_figure_data.py` |

## Notebooks vs scripts

- For a from-scratch reproduction of the paper's headline outputs (Figures 7–9, 11 and S4) from the predictions, run `supplement.ipynb` end-to-end (requires the restricted data).
- For just the SI parameter-grid artefacts (Tables S1–S2, Figure S3), `generate_si_inf_grid.py` is faster (~10 s) and reads bundled pickles only — no external download required.
- `INF_notebook_clean.ipynb` is the full exploratory notebook with the original grid searches; provided for reference only.

## The `metapaths` package

The `metapaths` Python package implementing the INF transformation is
bundled here as `metapaths_src/`. This bundled copy is a frozen
publication snapshot pinned to the version used in the paper. The
canonical public source of the package is:

> **https://github.com/TERENTIVS/metapaths_publ**

(active development happens in a private repository; the
`metapaths_publ` repository is the publication-aligned mirror.)

The Quickstart above installs from the bundled snapshot — recommended
for exact reproduction of the paper's numbers. To track the
publication snapshot from GitHub instead, replace step 1's package
install with:

```bash
pip install git+https://github.com/TERENTIVS/metapaths_publ@v1.0.0
```

(adjust the tag to the desired release).

> **Note.** The `metapaths` package is **not currently on PyPI** — it was withdrawn pending publication of this paper. Any `pip install metapaths` instructions a reader may encounter elsewhere (including any older copies of the upstream README) are stale; install from the bundled snapshot or from the GitHub URL above instead.

## Data sources

- **TSA knowledge graph.** Built on the Drug Repurposing Knowledge Graph (Ioannidis et al., 2020) and augmented in-house at Evotec; node mapping (`data/hetionet-v1.0-nodes_version2.tsv`) bundled.
- **Ground-truth target–AE associations.** Curated at Evotec (6,657 associations) and proprietary; not included in this repository. A 200-association sample illustrating the ground-truth schema, and all pickles containing labelled target–AE pairs, are available from the first author on request — see `data/pickles/DOWNLOAD.md`.
- **External validation set.** Soldatos et al. (2018) supplementary files 2, 3, 4 bundled in `data/soldatos_2018/` under their original CC-BY 4.0 licence — credit and the original publication are noted in `data/soldatos_2018/README.md`.

## Licence

Code and documentation are released under the MIT License; see
[`LICENSE`](./LICENSE). The bundled Soldatos et al. (2018) files keep
their original CC BY 4.0 licence.

## Citation

```bibtex
@article{Egbelo2026INF,
  title  = {Improving Target--Adverse Event Association Prediction via
            Inverse Node Frequency Transformation to Mitigate Knowledge
            Graph Topological Bias},
  author = {Egbelo, Terence and Jeynes, J Charles G and Kurt, Zeyneb and
            Veettil, Rahul V and Gillet, Valerie J},
  journal = {Bioinformatics},
  year    = {2026},
  note    = {Forthcoming}
}
```

## Authors

- Terence Egbelo¹·² — University of Sheffield · Evotec UK · terence.egbelo@gmail.com
- J Charles G Jeynes² — Evotec UK
- Zeyneb Kurt¹ — University of Sheffield
- Rahul V Veettil² — Evotec UK
- Valerie J Gillet¹·* — University of Sheffield · v.gillet@sheffield.ac.uk

¹ School of Information, Journalism and Communication, University of Sheffield.
² Evotec UK, Abingdon.
\* Corresponding author.
