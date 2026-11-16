"""
config.py — Project-level configuration for INF_notebook.ipynb
===============================================================
All paths are relative to this file's location by default, so the
package works out of the box; the restricted pickle files (available from
the corresponding author on request) go in data/pickles/.  See
data/pickles/DOWNLOAD.md for the file list.

The only settings that ever need changing are:
  - PICKLES_DIR  — if you keep the large pickles elsewhere
  - NEO4J_URI_*  — only needed to regenerate pickles from scratch
"""

from pathlib import Path

# Project root (directory containing this file)
_HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------

# Saved figures (filenames follow manuscript numbering: figN_* = Figure N, figS*/tableS* = supplement)
FIGS_DIR = _HERE / "figures"

# ---------------------------------------------------------------------------
# Neo4j connection URIs
# ---------------------------------------------------------------------------
# Required only to regenerate pickles from scratch (not for figure reproduction).
# Update the port to match your local TSA graph instance.
NEO4J_URI_READ  = "bolt://localhost:7687"
NEO4J_URI_WRITE = "bolt://localhost:7687"

# ---------------------------------------------------------------------------
# Root data directories  (relative paths — no editing required)
# ---------------------------------------------------------------------------

# All pickle files live here.  Small pickles are bundled; large ones are
# available on request — see data/pickles/DOWNLOAD.md.
PICKLES_DIR = _HERE / "data" / "pickles"

# Soldatos et al. (2018) supplementary files — bundled in data/soldatos_2018/
SOLDATOS_DIR = _HERE / "data" / "soldatos_2018"

# Aggregate values behind figures whose inputs are restricted (Figures 8, 11, S4)
FIGURE_DATA_DIR = _HERE / "data" / "figure_data"

# Hetionet v1.0 node mapping (Entrez ID → common name) — bundled in data/
MAPPING_FILE = _HERE / "data" / "hetionet-v1.0-nodes_version2.tsv"

# ---------------------------------------------------------------------------
# Pickle file paths  (all derived from PICKLES_DIR — no editing required)
# ---------------------------------------------------------------------------

# --- Shared / global ---
PARAM_COMBOS_PKL        = PICKLES_DIR / "param_combos.pkl"
ALL_GENES_AES_PKL       = PICKLES_DIR / "all_genes_with_aes.pkl"
TEST_POS_PAIRS_PKL      = PICKLES_DIR / "test_pos_pairs.pkl"

# --- CV grid-search results (Section I) ---
CV_GS_RESULTS_L3_PKL        = PICKLES_DIR / "cv_gs_results_l3.pkl"
CV_GS_RESULTS_L3_BOTTOM_PKL = PICKLES_DIR / "cv_gs_results_l3_bottom.pkl"
CV_GS_RESULTS_L4_PKL        = PICKLES_DIR / "cv_gs_results_l4.pkl"
CV_GS_RESULTS_L3_DWPC_PKL   = PICKLES_DIR / "cv_gs_results_l3_dwpc_recall.pkl"
CV_GS_RESULTS_L4_F1_PKL     = PICKLES_DIR / "cv_gs_results_l4_f1.pkl"

# --- CV INF transformation feature matrices (Section I grid search) ---
CV_INF_FEATS_L3_PKL         = PICKLES_DIR / "cv_inf_feats_l3_all_combos.pkl"
CV_INF_FEATS_L34_PKL        = PICKLES_DIR / "cv_inf_feats_l34_all_combos.pkl"
CV_DWPC_FEATS_PKL           = PICKLES_DIR / "cv_dwpc_feats_ready.pkl"
DRKG_FEATS_L3_PKL           = PICKLES_DIR / "drkg_feats_l3.pkl"

# --- CV slice results (Section I) ---
BEST_CONFIGS_L3_RECALL_PKL = PICKLES_DIR / "best_configs_l3_recall.pkl"
BEST_CONFIGS_L4_RECALL_PKL = PICKLES_DIR / "best_configs_l4_recall.pkl"
RAW_L3_RECALL_PKL          = PICKLES_DIR / "raw_l3_recall_copy.pkl"
RAW_L4_RECALL_PKL          = PICKLES_DIR / "raw_l4_recall_copy.pkl"

# --- INF / DWPC pre-computed features (Section II) ---
INF_FEATS_DIR         = PICKLES_DIR / "inf_features"
DWPC_FEATS_DIR        = PICKLES_DIR / "dwpc_features"
DWPC_NEGS_TO_TEST_PKL = PICKLES_DIR / "negs_to_test_l3_feats_ids.pkl"
HITK_PARTITIONS_PKL   = PICKLES_DIR / "hitk_partitions_ae_deg.pkl"

# --- Feature importance results (Section II) ---
FEAT_IMPORTANCE_DIR = PICKLES_DIR
FIG10_FEAT_IMP_PKL = PICKLES_DIR / "config_top_feats_by_slice.pkl"   # Figure 10 (holdout, 1st decile)

# --- Holdout test results (Section II) ---
HITK_DFS_PKL             = PICKLES_DIR / "hitk_dfs.pkl"
SINGLE_BEST_HITK_PKL     = PICKLES_DIR / "single_best_hitk.pkl"
SINGLE_BEST_RECALL_PKL   = PICKLES_DIR / "single_best_recall.pkl"
RAW_PREDS_PKL            = PICKLES_DIR / "raw_preds.pkl"
DWPC_PREDS_PKL           = PICKLES_DIR / "dwpc_preds.pkl"
IDEAL_INF_RES_PKL        = PICKLES_DIR / "ideal_inf_res.pkl"          # output only
IDEAL_INF_RES_RECALL_PKL = PICKLES_DIR / "ideal_inf_res_recall.pkl"
INF_CONFIGS_SLICE_RECALL_PKL = PICKLES_DIR / "inf_configs_slice_recall.pkl"

# --- Soldatos external test (Section III) ---
SOLDATOS_INTERSECTS_PKL   = PICKLES_DIR / "soldatos_drkg_intersects.pkl"
SOLDATOS_NEW_TEST_PKL     = PICKLES_DIR / "new_test.pkl"
SOLDATOS_PREDS_MAPPED_PKL = PICKLES_DIR / "preds_mapped_incl_ensemble.pkl"
SOLDATOS_RESULTS_PKL      = PICKLES_DIR / "res_soldatos_incl_recall.pkl"
SOLDATOS_HITK_MAPPED_PKL  = PICKLES_DIR / "hitk_dfs_mapped_newtest_soldatos.pkl"

# --- Soldatos supplementary data files ---
SOLDATOS_FILE_2 = SOLDATOS_DIR / "Supplementary_File_2.txt"
SOLDATOS_FILE_3 = SOLDATOS_DIR / "Supplementary_File_3.txt"
SOLDATOS_FILE_4 = SOLDATOS_DIR / "Supplementary_File_4.txt"
