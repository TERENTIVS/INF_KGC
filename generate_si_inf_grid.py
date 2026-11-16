"""
Generate SI artefacts addressing reviewer comments on §3.2.

Produces three components driven entirely by bundled pickles
(no Neo4j, no model retraining):

1. ``figures/tableS1_inf_grid_heatmap.png`` — 35 INF configurations
   x 10 deciles of CV TPR, with an "Overall" right-hand column. The
   per-decile maximum is marked with a white asterisk. This is the
   "different INF results in supporting information" the reviewer
   asked for in §3.2.

2. ``data/tableS2_inf_grid_top3.csv`` and ``.md`` — top 3 configurations
   per decile, with full parameter triple (deflator exponent, INF
   aggregator, pooling) and 95% bootstrap CI. Used both to support
   the deflator-exponent-0 worked example at decile 1 and to make the
   ensemble's per-decile winners auditable.

3. ``figures/figS3_decile_best_vs_global_vs_raw.png`` — grouped bar
   chart contrasting, for every decile, the decile-best INF config
   against the best single global INF config (the one chosen as the
   default "best INF Config" in the main text) and the raw-count
   baseline. This is what supports the §3.2 ensemble paragraph the
   reviewer flagged as "expand and show data".

All three artefacts use the **same** evaluator that produced the §3.2
headline number 0.424 — ``analyse_parameter_slices_quantile`` from
``tests/helpers.py`` (verbatim of the notebook function in
``INF_notebook_clean.ipynb``). This is the CV-fold-mean estimator:
per fold, deciles are recomputed via ``pd.qcut`` on the positives and
sklearn ``recall_score`` is computed on the positives plus their
degree-matched negatives; across folds, the point estimate is the mean
of the 5 per-fold recall values and the CI is the t-based half-width
(``ci_margin_of_error``).

This is the same evaluator used to produce the headline §3.2 values
(0.424 / 0.202 / 0.263 at D1) and Figure 7 of the main text, so the SI
sits on a single, internally consistent footing.

Sources
-------
- ``cv_gs_results_l3.pkl`` — dict[int config_idx] -> (summary tuple,
  fold_results, model_params). ``fold_results`` is the input to
  ``analyse_parameter_slices_quantile``.
- ``best_configs_l3_recall.pkl`` — DataFrame of the decile-best INF
  configs as reported in §3.2 (recall_mean ± recall_ci half-width).
- ``raw_l3_recall_copy.pkl`` — pre-computed config-0 (raw count) per-
  decile recall via the same fold-mean estimator. Used directly for the
  raw-baseline bars; for the heatmap we re-derive the same values from
  ``cv_gs_results_l3.pkl[0]``.
- ``param_combos.pkl`` — list of the 35 parameter dicts.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

import config
from metapaths_src.metapaths.starterpack import load_pickle
from tests.helpers import analyse_parameter_slices_quantile

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams


# ── Plot styling (matches generate_figs7_8_9.py) ───────────────────────────
rcParams.update({
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

OKABE_ITO = {
    "decile_best": "#0072B2",   # blue
    "best_global": "#009E73",   # green
    "raw":         "#999999",   # grey
}

FIGS_DIR = Path(config.FIGS_DIR)
DATA_DIR = Path(config.PICKLES_DIR).parent  # project_root/data
FIGS_DIR.mkdir(parents=True, exist_ok=True)


# ── Helpers ────────────────────────────────────────────────────────────────
N_CONFIGS = 35
N_DECILES = 10


def _format_params(params: dict) -> str:
    """One-line label for a parameter triple."""
    defl = params.get("path_deflator_exp")
    infl = params.get("inf_inflator")
    pool = params.get("inf_pooling")
    if infl is None and defl is None:
        return "raw counts (baseline)"
    if infl is None:
        return f"defl={defl}, agg=—, pool=—  (deflated count only)"
    defl_str = "—" if defl is None else f"{defl}"
    pool_str = "—" if pool is None else pool
    return f"defl={defl_str}, agg={infl}, pool={pool_str}"


def _load_grid():
    """Run analyse_parameter_slices_quantile per config to populate the grid.

    Returns (matrix, ci_radius, overall, param_combos):
      matrix      : (35, 10) per-decile mean recall
      ci_radius   : (35, 10) per-decile t-based half-width (symmetric)
      overall     : (35,)    mean across deciles per config
      param_combos: list of 35 parameter dicts
    """
    cv = load_pickle(config.CV_GS_RESULTS_L3_PKL)
    pc = load_pickle(config.PARAM_COMBOS_PKL)
    assert len(cv) == N_CONFIGS, f"expected {N_CONFIGS} configs, got {len(cv)}"

    matrix    = np.full((N_CONFIGS, N_DECILES), np.nan)
    ci_radius = np.full((N_CONFIGS, N_DECILES), np.nan)

    for cfg_idx in range(N_CONFIGS):
        df = analyse_parameter_slices_quantile(cv[cfg_idx][1])
        for _, row in df.iterrows():
            d = int(row["slice_num"]) - 1
            matrix[cfg_idx, d]    = row["recall_mean"]
            ci_radius[cfg_idx, d] = row["recall_ci"]

    overall = np.nanmean(matrix, axis=1)
    return matrix, ci_radius, overall, pc


# ── Artefact 1: heatmap ────────────────────────────────────────────────────
def make_heatmap(out_path: Path = None) -> Path:
    matrix, _, overall, pc = _load_grid()

    # Compose (configs x (deciles + overall)) display matrix
    display = np.column_stack([matrix, overall])
    col_labels = [f"D{d}" for d in range(1, N_DECILES + 1)] + ["Overall"]
    row_labels = [f"{i:>2} | {_format_params(pc[i])}" for i in range(N_CONFIGS)]

    fig, ax = plt.subplots(figsize=(13, 14))
    im = ax.imshow(display, aspect="auto", cmap="viridis", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(N_CONFIGS))
    ax.set_yticklabels(row_labels, fontsize=8, family="monospace")

    # Per-cell 2-decimal value annotation. Light text on dark cells, dark on light.
    for r in range(N_CONFIGS):
        for c in range(display.shape[1]):
            val = display[r, c]
            if np.isnan(val):
                continue
            ax.text(c, r, f"{val:.2f}",
                    ha="center", va="center",
                    color="white" if val < 0.55 else "black",
                    fontsize=6.5)

    # Asterisk the per-decile maximum within the INF configs (idx >= 1) only,
    # so the raw baseline is shown but not artificially flagged as "winner".
    for d in range(N_DECILES):
        best_cfg = 1 + int(np.argmax(matrix[1:, d]))
        ax.text(d + 0.30, best_cfg - 0.30, "*",
                ha="center", va="center",
                color="red", fontsize=14, fontweight="bold")
    best_overall = 1 + int(np.argmax(overall[1:]))
    ax.text(N_DECILES + 0.30, best_overall - 0.30, "*",
            ha="center", va="center",
            color="red", fontsize=14, fontweight="bold")

    # Separator between deciles and overall column
    ax.axvline(N_DECILES - 0.5, color="white", linewidth=1.5)

    ax.set_xlabel("Decile of training-set AE count (D1 = 0–5 AEs, D10 = 163–256 AEs)")
    ax.set_title("CV TPR (recall) for all 35 INF configurations × 10 deciles\n"
                 "Fold-mean estimator (analyse_parameter_slices_quantile); "
                 "* marks the best INF config in each column", pad=10)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Mean CV TPR")

    fig.tight_layout()
    out_path = out_path or (FIGS_DIR / "tableS1_inf_grid_heatmap.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── Artefact 2: top-3 per decile table ─────────────────────────────────────
def make_top3_table(csv_path: Path = None, md_path: Path = None):
    matrix, ci_radius, _, pc = _load_grid()

    rows = []
    for d in range(1, N_DECILES + 1):
        col = matrix[:, d - 1]
        order = np.argsort(-col)  # descending
        for rank, cfg_idx in enumerate(order[:3], start=1):
            mean = float(col[cfg_idx])
            err  = float(ci_radius[cfg_idx, d - 1])
            rows.append({
                "decile": d,
                "rank": rank,
                "config_idx": int(cfg_idx),
                "deflator_exp": pc[cfg_idx].get("path_deflator_exp"),
                "inf_inflator": pc[cfg_idx].get("inf_inflator"),
                "inf_pooling":  pc[cfg_idx].get("inf_pooling"),
                "cv_tpr": round(mean, 4),
                "ci_half_width": round(err, 4),
            })
    df = pd.DataFrame(rows)

    csv_path = csv_path or (DATA_DIR / "tableS2_inf_grid_top3.csv")
    md_path  = md_path  or (DATA_DIR / "tableS2_inf_grid_top3.md")
    df.to_csv(csv_path, index=False)

    # Markdown table
    md_lines = [
        "# Top-3 INF configurations per decile (CV TPR)",
        "",
        "Decile 1 corresponds to targets with 0–5 known AEs (sparsest); "
        "decile 10 corresponds to 163–256 AEs (densest). 'config_idx' is the "
        "position of the parameter triple in `param_combos.pkl`. "
        "Source: ``analyse_parameter_slices_quantile`` applied per fold to "
        "``cv_gs_results_l3.pkl`` (fold-mean estimator, t-based 95% half-width).",
        "",
        "| decile | rank | config | deflator exp | INF aggregator | pooling | CV TPR | ± 95% CI |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        defl = "—" if r["deflator_exp"] is None else r["deflator_exp"]
        infl = "—" if r["inf_inflator"] is None else r["inf_inflator"]
        pool = "—" if r["inf_pooling"]  is None else r["inf_pooling"]
        md_lines.append(
            f"| {r['decile']} | {r['rank']} | {r['config_idx']} | "
            f"{defl} | {infl} | {pool} | {r['cv_tpr']:.2f} | ± {r['ci_half_width']:.2f} |"
        )
    md_path.write_text("\n".join(md_lines) + "\n")
    return csv_path, md_path, df


# ── Artefact 3: decile-best vs best-global vs raw bar chart ────────────────
def make_decile_comparison_bar(out_path: Path = None) -> Path:
    """Grouped bar chart per decile.

    All three series use the same fold-mean estimator
    (``analyse_parameter_slices_quantile``):

    * Decile-best from ``best_configs_l3_recall.pkl`` (matches §3.2's
      headline 0.424 for D1).
    * Best-global from running the same estimator on
      ``cv_gs_results_l3.pkl[global_best_cfg]``, where global_best_cfg
      is selected by mean-across-deciles CV TPR (excluding config 0).
    * Raw baseline from ``raw_l3_recall_copy.pkl``, which is the pre-
      stored output of the same estimator on config 0.
    """
    best_configs = load_pickle(config.BEST_CONFIGS_L3_RECALL_PKL)
    raw_l3 = (
        load_pickle(config.RAW_L3_RECALL_PKL)
        .reset_index(drop=True)
    )
    matrix, ci_radius, overall, pc = _load_grid()

    # Identify global best by mean-across-deciles CV TPR, excluding raw baseline
    global_best_cfg = 1 + int(np.argmax(overall[1:]))

    # Decile-best (from main-text pickle)
    bc = best_configs.sort_values("slice_num").reset_index(drop=True)
    db_means = bc["recall_mean"].to_numpy(dtype=float)
    db_err = bc["recall_ci"].to_numpy(dtype=float)  # symmetric half-width
    db_cfg_idx = bc["param_combo"].to_numpy(dtype=int)

    # Best-global from fold-mean grid
    bg_means = matrix[global_best_cfg]
    bg_err   = ci_radius[global_best_cfg]

    # Raw baseline from raw_l3_recall_copy.pkl (same estimator, pre-stored)
    raw_means = raw_l3["recall_mean"].to_numpy(dtype=float)
    raw_err   = raw_l3["recall_ci"].to_numpy(dtype=float)

    x = np.arange(N_DECILES)
    width = 0.27

    fig, ax = plt.subplots(figsize=(12, 6))

    bars_db = ax.bar(x - width, db_means, width, yerr=db_err,
                     color=OKABE_ITO["decile_best"], capsize=3,
                     label="Decile-best INF config")
    bars_bg = ax.bar(x, bg_means, width, yerr=bg_err,
                     color=OKABE_ITO["best_global"], capsize=3,
                     label=f"Best single global INF config "
                           f"(#{global_best_cfg}: {_format_params(pc[global_best_cfg])})")
    bars_rw = ax.bar(x + width, raw_means, width, yerr=raw_err,
                     color=OKABE_ITO["raw"], capsize=3,
                     label="Raw-count baseline (config 0)")

    def _annotate(bars, means, errs, color):
        for bar, mean, err in zip(bars, means, errs):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    mean + err + 0.012,
                    f"{mean:.2f}",
                    ha="center", va="bottom",
                    fontsize=8, fontweight="bold", color=color)

    _annotate(bars_db, db_means, db_err, OKABE_ITO["decile_best"])
    _annotate(bars_bg, bg_means, bg_err, OKABE_ITO["best_global"])
    _annotate(bars_rw, raw_means, raw_err, OKABE_ITO["raw"])

    # Annotate decile-best bars with the underlying config index, to make
    # the ensemble composition auditable from the chart alone.
    for i, cfg_idx in enumerate(db_cfg_idx):
        ax.text(x[i] - width, db_means[i] - 0.06, f"#{cfg_idx}",
                ha="center", va="top", fontsize=7,
                color="white", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"D{d}\n[{int(bc.iloc[d-1]['min_param'])}-"
                        f"{int(bc.iloc[d-1]['max_param'])}]"
                        for d in range(1, N_DECILES + 1)])
    ax.set_xlabel("Decile (target AE count range shown in brackets)")
    ax.set_ylabel("CV TPR (recall)")
    ax.set_ylim(0, 1.15)
    ax.set_title("Per-decile CV TPR: decile-best INF ensemble components "
                 "vs best single global INF config vs raw-count baseline")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    out_path = out_path or (FIGS_DIR / "figS3_decile_best_vs_global_vs_raw.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── Entry point ────────────────────────────────────────────────────────────
def main():
    heatmap = make_heatmap()
    csv_path, md_path, _ = make_top3_table()
    bar = make_decile_comparison_bar()
    print("Wrote:")
    for p in (heatmap, csv_path, md_path, bar):
        print(f"  {p}")


if __name__ == "__main__":
    main()
