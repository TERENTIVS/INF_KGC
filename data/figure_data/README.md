# Figure data — aggregate values for figures with restricted inputs

Figures 8, 11 and Supplementary Figure S4 are computed from prediction
pickles that contain the proprietary Evotec target–AE labels (restricted
access; see `data/pickles/DOWNLOAD.md`). The values plotted in those
figures are provided here so they can be inspected and re-plotted
without the restricted files. No target or AE identifiers are included.

| File | Manuscript item | Contents |
|---|---|---|
| `fig8_holdout_d1.csv` | Figure 8 | Holdout test set, 1st decile (0–5 AEs): Hits@10 and TPR with 95% CIs for Best INF Config, DWPC, Original |
| `fig11_external_feat_importance_d1.csv` | Figure 11 | Soldatos external test set, 1st decile: top-5 metapath features by permutation importance (mean decrease in TPR) for INF config 12, Original (config 0) and DWPC |
| `figS4_holdout_vs_external.csv` | Supplementary Figure S4 | All targets: Hits@10 and TPR with 95% CIs, holdout (n = 998) and external (n = 7,804) test sets |

**Provenance.** Exported from a full execution of `supplement.ipynb`
against the restricted pickles (16 Nov 2026). The Figure 11 values are
identical to the stored `config_top_feats_by_slice_soldatos.pkl`.

**Run-to-run variation.** The notebook's bootstrap (`calculate_hits_at_k`,
`calculate_recall`) is unseeded, and the INF ensemble in Figure S4 is
retrained on each run. Values recomputed from the restricted files may
therefore differ from these in the third decimal place (CI bounds by up
to ~0.01); all figures in the manuscript are reported to two decimals.

`generate_figs7_8_9.py` uses `fig8_holdout_d1.csv` automatically when
the restricted holdout predictions are absent.
