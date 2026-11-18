# Restricted data — available on request

The files listed below contain the proprietary target–adverse event (AE)
associations curated by Evotec (UK) Ltd, as labelled target–AE pairs,
per-target AE counts or the list of annotated targets. They are not
redistributed publicly and are available from the first author
(Terence Egbelo, terence.egbelo@gmail.com) on reasonable request. Please
state your name, affiliation and intended use. The files are provided
for verification and reproduction of the published results and may not
be redistributed.

## What works without these files

Everything in the public repository runs on the bundled data alone:

- all reproduction tests (`pytest tests/`)
- Figures 4, 6, 7, 9 and 10, Supplementary Tables S1–S2 and Figure S3
  (bundled pickles; `generate_figs7_8_9.py`, `generate_si_inf_grid.py`)
- Figure 8 (`generate_figs7_8_9.py` falls back to the bundled values
  in `data/figure_data/`)
- the plotted values of Figures 8, 11 and Supplementary Figure S4
  (`data/figure_data/`)

The restricted files are needed to recompute Figures 8, 11 and S4 from
the underlying predictions, i.e. to run `supplement.ipynb` end to end.

## Files available on request

Place the `.pkl` files in this directory (`data/pickles/`) and unzip the
two archives here, which creates `dwpc_features/` and `inf_features/`.

| File | Size | Contents / used for |
|---|---|---|
| `cv_inf_feats_l3_all_combos.pkl` | ~372 MB | L3 INF feature matrices for the CV grid search |
| `cv_inf_feats_l34_all_combos.pkl` | ~856 MB | L3+L4 INF feature matrices for the CV grid search |
| `cv_dwpc_feats_ready.pkl` | ~11 MB | DWPC feature matrices for the CV folds |
| `hitk_dfs.pkl` | ~1.4 GB | Trained RF models + holdout predictions, all INF configs |
| `hitk_partitions_ae_deg.pkl` | ~1.6 GB | INF feature matrices (train/test partitions, all 35 configs) |
| `inf_configs_slice_recall.pkl` | ~1.6 GB | Slice-wise holdout recall, all INF configs |
| `raw_preds.pkl` | ~48 MB | Holdout predictions — original metapath counts |
| `dwpc_preds.pkl` | ~43 MB | Holdout predictions — DWPC |
| `single_best_hitk.pkl` | ~48 MB | INF ensemble holdout predictions (Hits@k optimised) |
| `single_best_recall.pkl` | ~48 MB | INF ensemble holdout predictions (recall optimised) |
| `ideal_inf_res.pkl` | ~90 MB | INF ensemble holdout predictions and metrics |
| `negs_to_test_l3_feats_ids.pkl` | ~39 MB | DWPC feature matrix for holdout evaluation |
| `test_pos_pairs.pkl` | ~56 KB | Holdout positive target–AE pairs |
| `all_genes_with_aes.pkl` | ~8 KB | Targets with ≥1 curated AE |
| `soldatos_drkg_intersects.pkl` | ~37 MB | Soldatos–DRKG intersection with Evotec AE counts |
| `new_test.pkl` | ~36 MB | External (Soldatos) test set with Evotec AE counts |
| `preds_mapped_incl_ensemble.pkl` | ~658 MB | All predictions mapped to the Soldatos test set |
| `hitk_dfs_mapped_newtest_soldatos.pkl` | ~640 MB | Soldatos evaluation DataFrames |
| `dwpc_features.zip` | ~41 MB | DWPC train/test feature files (unzip here) |
| `inf_features.zip` | ~1.3 GB | INF train/test feature files (unzip here) |
| `ground_truth_sample_200.csv` | <1 MB | 200 randomly sampled curated associations, illustrating the ground-truth schema (Entrez Gene, MedDRA PT, common names) |

A `SHA256SUMS` file listing checksums for every file is supplied with them.
