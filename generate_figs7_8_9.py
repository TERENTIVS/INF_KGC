"""
Generate manuscript Figures 7, 8 and 9 (both vertical and horizontal bar versions).

Figure 7: CV TPR for 1st decile (0-5 AEs) — Original vs DWPC vs Best INF config
Figure 8: Holdout test set, 1st decile — Hits@10 and TPR side-by-side, 3 methods
Figure 9: External Soldatos test set — Hits@10 and TPR side-by-side, 3 methods

For each figure: one vertical bar chart and one horizontal bar chart saved as PNG.
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import config
from metapaths_src.metapaths.starterpack import load_pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from typing import Tuple, List
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import recall_score as sk_recall_score

# ── Accessibility-friendly palette (Okabe-Ito) ──────────────────────────────
COLORS = {
    'Best INF Config': '#0072B2',   # blue
    'DWPC':            '#E69F00',   # orange/amber
    'Original':        '#999999',   # grey
}
# Fig 8/9 use this order; Fig 7 has its own method list
METHOD_ORDER = ['Best INF Config', 'DWPC', 'Original']

rcParams.update({
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
})

FIGS_DIR = config.FIGS_DIR
os.makedirs(str(FIGS_DIR), exist_ok=True)


# ── Metric implementations matching INF_notebook.ipynb exactly ──────────────
#
# calculate_hits_at_k: filtered Hits@k (standard KGE evaluation protocol).
#   - Per test pair (not per gene).
#   - Filtered rank: rank of true tail after removing other known positives
#     ranked above it, so they don't inflate the rank.
#   - Genes with no true positives in this slice are skipped (not counted).
#   - Bootstrap resamples individual pair ranks.
#
# calculate_recall: cluster-bootstrap TPR.
#   - Per gene: records (TP, FN) for each gene that has ≥1 true positive.
#   - Bootstrap resamples gene-level records, then computes aggregate recall.

def calculate_hits_at_k(df: pd.DataFrame, k: int, head_col: str, tail_col: str,
                         n_bootstrap: int = 1000) -> Tuple[float, Tuple[float, float]]:
    all_filtered_ranks = []
    if df.empty:
        return 0.0, (0.0, 0.0)
    for head, group in df.groupby(head_col):
        sorted_group = group.sort_values('model_score', ascending=False)
        true_tails = set(sorted_group[sorted_group['true_label'] == True][tail_col])
        if not true_tails:
            continue
        ranked_tails = sorted_group[tail_col].tolist()
        for true_t in true_tails:
            try:
                raw_rank_idx = ranked_tails.index(true_t)
            except ValueError:
                continue
            n_other_true_higher = len([t for t in ranked_tails[:raw_rank_idx]
                                       if t in true_tails])
            filtered_rank = (raw_rank_idx + 1) - n_other_true_higher
            all_filtered_ranks.append(filtered_rank)
    if not all_filtered_ranks:
        return 0.0, (0.0, 0.0)
    filtered_ranks = np.array(all_filtered_ranks)
    hits_score = np.mean(filtered_ranks <= k)
    boot = np.random.default_rng(42).choice(filtered_ranks,
                                             (n_bootstrap, len(filtered_ranks)),
                                             replace=True)
    boot_means = np.mean(boot <= k, axis=1)
    return hits_score, (float(np.percentile(boot_means, 2.5)),
                        float(np.percentile(boot_means, 97.5)))


def calculate_recall(df: pd.DataFrame, head_col: str, tail_col: str,
                     n_bootstrap: int = 1000) -> Tuple[float, Tuple[float, float]]:
    if df.empty:
        return 0.0, (0.0, 0.0)
    per_gene_records = []
    for head, group in df.groupby(head_col):
        tp_rows = group[group['true_label'] == True]
        if tp_rows.empty:
            continue
        pred_pos = group[group['model_score'] >= 0.5]
        tp = len(pred_pos.merge(tp_rows, on=[head_col, tail_col], how='inner'))
        fn = len(tp_rows) - tp
        per_gene_records.append({'tp': tp, 'fn': fn})
    if not per_gene_records:
        return 0.0, (0.0, 0.0)
    rng = np.random.default_rng(42)
    recall_scores = []
    idx = np.arange(len(per_gene_records))
    for _ in range(n_bootstrap):
        samp = [per_gene_records[i] for i in rng.choice(idx, size=len(idx), replace=True)]
        tot_tp = sum(r['tp'] for r in samp)
        tot_fn = sum(r['fn'] for r in samp)
        recall_scores.append(tot_tp / (tot_tp + tot_fn) if (tot_tp + tot_fn) > 0 else 0.0)
    return (float(np.mean(recall_scores)),
            (float(np.percentile(recall_scores, 2.5)),
             float(np.percentile(recall_scores, 97.5))))


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9 DATA
# Source: CV results pickles; 1st decile (0–5 AEs per target)
# ─────────────────────────────────────────────────────────────────────────────
print("Loading Figure 7 data (CV 1st decile TPR)...")

best_configs_l3 = load_pickle(str(config.BEST_CONFIGS_L3_RECALL_PKL))
raw_l3_recall   = load_pickle(str(config.RAW_L3_RECALL_PKL))
cv_gs_dwpc      = load_pickle(str(config.CV_GS_RESULTS_L3_DWPC_PKL))

# Best INF: 1st decile row from best_configs_l3
inf_d1 = best_configs_l3[best_configs_l3['slice_num'] == 1].iloc[0]
inf_d1_mean = inf_d1['recall_mean']
inf_d1_ci   = inf_d1['recall_ci']        # ± half-width (already symmetric)

# Original: 1st decile row from raw_l3_recall (index 0 = decile 1)
orig_d1_mean = raw_l3_recall.iloc[0]['recall_mean']
orig_d1_ci   = raw_l3_recall.iloc[0]['recall_ci']

# DWPC: use analyse_parameter_slices_quantile (quantile-based per-fold boundaries),
# matching the source notebook's approach exactly.
def _analyse_parameter_slices_quantile(fold_results, num_slices=10):
    """Replica of the source notebook's analyse_parameter_slices_quantile."""
    slice_data = []
    for fold_id, fold_data in fold_results.items():
        fold_df = pd.DataFrame({
            'y_val':         fold_data['y_val'],
            'y_pred':        fold_data['y_pred'],
            'param_values':  fold_data['parameter_values'],
        })
        pos_df = fold_df[fold_df['y_val'] == 1].copy()
        if pos_df.empty:
            continue
        try:
            pos_df['slice_num'] = pd.qcut(pos_df['param_values'], q=num_slices,
                                          labels=False, duplicates='drop')
        except ValueError:
            continue
        for sid in sorted(pos_df['slice_num'].unique()):
            sl = pos_df[pos_df['slice_num'] == sid]
            mn, mx = sl['param_values'].min(), sl['param_values'].max()
            negs   = fold_df[(fold_df['y_val'] == 0) &
                              fold_df['param_values'].between(mn, mx)]
            final  = pd.concat([sl, negs])
            recall = sk_recall_score(final['y_val'], final['y_pred'], zero_division=0)
            slice_data.append({'slice_num': sid + 1, 'recall': recall,
                                'min_param': mn, 'max_param': mx})
    df = pd.DataFrame(slice_data)
    if df.empty:
        return pd.DataFrame()
    # Aggregate by slice_num only (min/max vary per fold)
    agg_recall = df.groupby('slice_num')['recall'].agg(['mean', 'sem', 'count']).reset_index()
    agg_recall.rename(columns={'mean': 'recall_mean', 'sem': 'recall_sem'}, inplace=True)
    agg_recall['recall_ci'] = agg_recall.apply(
        lambda r: r['recall_sem'] * stats.t.ppf(0.975, df=r['count'] - 1)
                  if r['count'] > 1 else 0.0, axis=1)
    # Representative param boundaries from the first fold
    bounds = (df.groupby('slice_num')[['min_param', 'max_param']].first().reset_index())
    agg = agg_recall.merge(bounds, on='slice_num')
    return agg[['min_param', 'max_param', 'recall_mean', 'recall_ci']]

dwpc_slice_res = _analyse_parameter_slices_quantile(cv_gs_dwpc[1])
dwpc_d1_mean = float(dwpc_slice_res.iloc[0]['recall_mean'])
dwpc_d1_ci   = float(dwpc_slice_res.iloc[0]['recall_ci'])

fig7_data = {
    'Method':   ['Original', 'DWPC', 'Best INF Config'],
    'TPR_mean': [orig_d1_mean, dwpc_d1_mean, inf_d1_mean],
    'TPR_ci':   [orig_d1_ci,   dwpc_d1_ci,   inf_d1_ci],
    'color':    [COLORS['Original'], COLORS['DWPC'], COLORS['Best INF Config']],
}

print(f"  Original:  {orig_d1_mean:.3f} ± {orig_d1_ci:.3f}")
print(f"  DWPC:      {dwpc_d1_mean:.3f} ± {dwpc_d1_ci:.3f}")
print(f"  Best INF:  {inf_d1_mean:.3f} ± {inf_d1_ci:.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8 DATA
# Source: holdout test set predictions; 1st decile (0–5 AEs)
# ─────────────────────────────────────────────────────────────────────────────
print("\nLoading Figure 8 data (holdout test, 1st decile)...")

# INF: use pre-computed config-12 predictions (stored in inf_configs_slice_recall).
#   TPR and CI come from the stored evaluation; Hits@10 is computed from predictions.
#   This avoids re-running the 35-second ensemble loop and gives a deterministic result.
# DWPC / Original: computed from pre-computed holdout prediction pickles filtered to d1.
#   Note: DWPC gives ~0.255 (not the source notebook's 0.252, which used an older 16-feature
#   DWPC matrix). The difference is 0.003; 0.255 is kept here for reproducibility.

FIG8_CSV = config.FIGURE_DATA_DIR / "fig8_holdout_d1.csv"
_fig8_restricted = [config.INF_CONFIGS_SLICE_RECALL_PKL, config.RAW_PREDS_PKL,
                    config.DWPC_PREDS_PKL]

if all(p.exists() for p in _fig8_restricted):
    # Recompute from the labelled holdout predictions (restricted-access files)
    inf_slice_recall = load_pickle(str(config.INF_CONFIGS_SLICE_RECALL_PKL))
    ideal_recall     = load_pickle(str(config.IDEAL_INF_RES_RECALL_PKL))
    raw_preds        = load_pickle(str(config.RAW_PREDS_PKL))
    dwpc_preds       = load_pickle(str(config.DWPC_PREDS_PKL))

    # Config 12 1st-decile stored TPR (deterministic, from inf_configs_slice_recall)
    _c12_m, c12_preds = inf_slice_recall[12]
    _c0_m,  _        = inf_slice_recall[0]

    # For the INF TPR/CI, use ideal_inf_res_recall (same values, confirmed consistent)
    inf6_rec    = ideal_recall['best_slice_scores']['slice_1.0_Recall']['score']
    inf6_rec_ci = ideal_recall['best_slice_scores']['slice_1.0_Recall']['ci']

    # For Original TPR/CI, use stored config-0 slice-1 values
    orig6_rec    = _c0_m['slice_1.0_Recall']
    orig6_rec_ci = _c0_m['slice_1.0_Recall_ci']

    # Hits@10 — computed from stored predictions (deterministic point estimate)
    inf6_d1  = c12_preds[c12_preds['Gene_ADR_deg'] <= 5].copy()
    raw_d1   = raw_preds[raw_preds['Gene_ADR_deg']   <= 5].copy()
    dwpc_d1  = dwpc_preds[dwpc_preds['Gene_ADR_deg'] <= 5].copy()

    inf6_h10,  inf6_h10_ci  = calculate_hits_at_k(inf6_d1,  10, 'Gene', 'ADR')
    orig6_h10, orig6_h10_ci = calculate_hits_at_k(raw_d1,   10, 'Gene', 'ADR')
    dwpc6_h10, dwpc6_h10_ci = calculate_hits_at_k(dwpc_d1,  10, 'Gene', 'ADR')
    dwpc6_rec, dwpc6_rec_ci = calculate_recall(   dwpc_d1,      'Gene', 'ADR')

    print(f"  Best INF Config: Hits@10={inf6_h10:.3f} ({inf6_h10_ci[0]:.3f},{inf6_h10_ci[1]:.3f}), "
          f"TPR={inf6_rec:.3f} ({inf6_rec_ci[0]:.3f},{inf6_rec_ci[1]:.3f})")
    print(f"  DWPC:            Hits@10={dwpc6_h10:.3f} ({dwpc6_h10_ci[0]:.3f},{dwpc6_h10_ci[1]:.3f}), "
          f"TPR={dwpc6_rec:.3f} ({dwpc6_rec_ci[0]:.3f},{dwpc6_rec_ci[1]:.3f})")
    print(f"  Original:        Hits@10={orig6_h10:.3f} ({orig6_h10_ci[0]:.3f},{orig6_h10_ci[1]:.3f}), "
          f"TPR={orig6_rec:.3f} ({orig6_rec_ci[0]:.3f},{orig6_rec_ci[1]:.3f})")

    fig8_data = {
        'Method': METHOD_ORDER,
        'Hits10_mean': [inf6_h10,  dwpc6_h10,  orig6_h10],
        'Hits10_lo':   [inf6_h10  - inf6_h10_ci[0],   dwpc6_h10  - dwpc6_h10_ci[0],   orig6_h10  - orig6_h10_ci[0]],
        'Hits10_hi':   [inf6_h10_ci[1]  - inf6_h10,   dwpc6_h10_ci[1]  - dwpc6_h10,   orig6_h10_ci[1]  - orig6_h10],
        'TPR_mean':    [inf6_rec,  dwpc6_rec,  orig6_rec],
        'TPR_lo':      [inf6_rec  - inf6_rec_ci[0],    dwpc6_rec  - dwpc6_rec_ci[0],   orig6_rec  - orig6_rec_ci[0]],
        'TPR_hi':      [inf6_rec_ci[1]  - inf6_rec,    dwpc6_rec_ci[1]  - dwpc6_rec,   orig6_rec_ci[1]  - orig6_rec],
        'color': [COLORS[m] for m in METHOD_ORDER],
    }
else:
    # Restricted files absent: use the bundled values exported from the
    # predictions (data/figure_data/fig8_holdout_d1.csv)
    print(f"  Restricted holdout predictions not found; using {FIG8_CSV.name}")
    _f8 = pd.read_csv(FIG8_CSV).set_index(["method", "metric"])
    fig8_data = {"Method": METHOD_ORDER, "color": [COLORS[m] for m in METHOD_ORDER]}
    for key, metric in [("Hits10", "Hits@10"), ("TPR", "TPR")]:
        rows = _f8.xs(metric, level="metric").loc[METHOD_ORDER]
        fig8_data[f"{key}_mean"] = rows["value"].tolist()
        fig8_data[f"{key}_lo"] = (rows["value"] - rows["ci_low"]).tolist()
        fig8_data[f"{key}_hi"] = (rows["ci_high"] - rows["value"]).tolist()
    for m in METHOD_ORDER:
        print(f"  {m}: Hits@10={_f8.loc[(m, 'Hits@10'), 'value']:.3f}, "
              f"TPR={_f8.loc[(m, 'TPR'), 'value']:.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9 DATA
# Source: Soldatos external test set (res_soldatos_incl_recall.pkl)
# ─────────────────────────────────────────────────────────────────────────────
print("\nLoading Figure 9 data (Soldatos external test set)...")

soldatos_res = load_pickle(str(config.SOLDATOS_RESULTS_PKL))

# keys: 'inf_ensemble', 'dwpc', 0 (original)
fig9_data = {
    'Method': METHOD_ORDER,
    'Hits10_mean': [], 'Hits10_lo': [], 'Hits10_hi': [],
    'TPR_mean':    [], 'TPR_lo':    [], 'TPR_hi':    [],
    'color': [COLORS[m] for m in METHOD_ORDER],
}

for key, label in [('inf_ensemble', 'INF Ensemble'), ('dwpc', 'DWPC'), (0, 'Original')]:
    h_val, h_ci = soldatos_res[key]['Hits@10']
    r_val, r_ci = soldatos_res[key]['Recall']
    fig9_data['Hits10_mean'].append(h_val)
    fig9_data['Hits10_lo'].append(h_val - h_ci[0])
    fig9_data['Hits10_hi'].append(h_ci[1] - h_val)
    fig9_data['TPR_mean'].append(r_val)
    fig9_data['TPR_lo'].append(r_val - r_ci[0])
    fig9_data['TPR_hi'].append(r_ci[1] - r_val)
    print(f"  {label}: Hits@10={h_val:.4f} ({h_ci[0]:.4f},{h_ci[1]:.4f}), TPR={r_val:.3f} ({r_ci[0]:.3f},{r_ci[1]:.3f})")


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def annotate_bar_vertical(ax, bar, lo_err, hi_err, fmt="{:.2f}", offset_frac=0.015,
                           fontsize=9.5):
    """Place label above the error bar cap, avoiding overlap."""
    bar_top = bar.get_height()
    error_top = bar_top + hi_err
    y_limit = ax.get_ylim()[1]
    # place text slightly above the error cap
    y_text = error_top + offset_frac * y_limit
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        y_text,
        fmt.format(bar_top),
        ha='center', va='bottom',
        fontsize=fontsize, fontweight='bold',
    )

def annotate_bar_horizontal(ax, bar, hi_err, fmt="{:.2f}", offset_frac=0.015,
                             fontsize=9.5):
    """Place label to the right of the error bar cap, avoiding overlap."""
    bar_right = bar.get_width()
    x_limit = ax.get_xlim()[1]
    x_text = bar_right + hi_err + offset_frac * x_limit
    ax.text(
        x_text,
        bar.get_y() + bar.get_height() / 2,
        fmt.format(bar_right),
        ha='left', va='center',
        fontsize=fontsize, fontweight='bold',
    )

def style_ax(ax, orientation='vertical'):
    ax.spines[['top', 'right']].set_visible(False)
    if orientation == 'vertical':
        ax.yaxis.grid(True, linestyle='--', alpha=0.45, zorder=0)
        ax.set_axisbelow(True)
    else:
        ax.xaxis.grid(True, linestyle='--', alpha=0.45, zorder=0)
        ax.set_axisbelow(True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 7  — Vertical
# ─────────────────────────────────────────────────────────────────────────────
def make_fig7_vertical():
    methods   = fig7_data['Method']
    means     = fig7_data['TPR_mean']
    cis       = fig7_data['TPR_ci']
    colors    = fig7_data['color']
    n         = len(methods)
    x         = np.arange(n)
    err_lo    = cis   # already half-widths
    err_hi    = cis

    fig, ax = plt.subplots(figsize=(6.5, 5))
    bars = ax.bar(x, means, width=0.55, color=colors, alpha=0.88, zorder=3,
                  yerr=[err_lo, err_hi], capsize=6,
                  ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})

    ax.set_ylim(0, max(means) * 1.35)
    for bar, lo, hi in zip(bars, err_lo, err_hi):
        annotate_bar_vertical(ax, bar, lo, hi, offset_frac=0.02)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=10.5)
    ax.set_ylabel('Mean CV TPR', fontsize=11)
    ax.set_xlabel('Method', fontsize=11)
    ax.set_title(
        'Mean CV TPR — 1st Decile of Target–AE Associations\n'
        'by No. of Known AEs for Target (0–5 AEs), L3 Metapaths',
        fontsize=11.5, fontweight='bold', pad=10
    )
    style_ax(ax, 'vertical')
    plt.tight_layout()
    out = str(FIGS_DIR / 'fig7_cv_d1_tpr_vertical.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig7_vertical()

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 7  — Horizontal
# ─────────────────────────────────────────────────────────────────────────────
def make_fig7_horizontal():
    methods   = fig7_data['Method']
    means     = fig7_data['TPR_mean']
    cis       = fig7_data['TPR_ci']
    colors    = fig7_data['color']
    n         = len(methods)
    y         = np.arange(n)
    err_lo    = cis
    err_hi    = cis

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(y, means, height=0.55, color=colors, alpha=0.88, zorder=3,
                   xerr=[err_lo, err_hi], capsize=6,
                   ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})

    ax.set_xlim(0, max(means) * 1.35)
    for bar, hi in zip(bars, err_hi):
        annotate_bar_horizontal(ax, bar, hi, offset_frac=0.025)

    ax.set_yticks(y)
    ax.set_yticklabels(methods, fontsize=10.5)
    ax.set_xlabel('Mean CV TPR', fontsize=11)
    ax.set_ylabel('Method', fontsize=11)
    ax.set_title(
        'Mean CV TPR — 1st Decile of Target–AE Associations\n'
        'by No. of Known AEs for Target (0–5 AEs), L3 Metapaths',
        fontsize=11.5, fontweight='bold', pad=10
    )
    style_ax(ax, 'horizontal')
    plt.tight_layout()
    out = str(FIGS_DIR / 'fig7_cv_d1_tpr_horizontal.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig7_horizontal()


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8 helper — dual panel (Hits@10 | TPR)
# ─────────────────────────────────────────────────────────────────────────────
def make_fig8_vertical():
    methods = fig8_data['Method']
    colors  = fig8_data['color']
    n = len(methods)
    x = np.arange(n)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.5), sharey=False)

    metrics = [
        ('Hits@10',
         fig8_data['Hits10_mean'], fig8_data['Hits10_lo'], fig8_data['Hits10_hi']),
        ('TPR',
         fig8_data['TPR_mean'],    fig8_data['TPR_lo'],    fig8_data['TPR_hi']),
    ]

    all_means = fig8_data['Hits10_mean'] + fig8_data['TPR_mean']
    shared_ylim = max(all_means) * 1.40

    for ax, (metric, means, lo, hi) in zip(axes, metrics):
        bars = ax.bar(x, means, width=0.55, color=colors, alpha=0.88, zorder=3,
                      yerr=[lo, hi], capsize=6,
                      ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})
        ax.set_ylim(0, shared_ylim)
        for bar, l, h in zip(bars, lo, hi):
            annotate_bar_vertical(ax, bar, l, h, offset_frac=0.025)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontsize=10)
        ax.set_ylabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        style_ax(ax, 'vertical')

    fig.suptitle(
        'Holdout Test Set — 1st Decile of Target–AE Associations\n'
        'by No. of Known AEs for Target (0–5 AEs); Best INF Config doubles TPR, Hits@10 is similar',
        fontsize=11, fontweight='bold', y=1.03
    )

    # Legend
    legend_patches = [
        mpatches.Patch(color=COLORS[m], alpha=0.88, label=m) for m in METHOD_ORDER
    ]
    fig.legend(handles=legend_patches, loc='lower center', ncol=3,
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.07))

    plt.tight_layout()
    out = str(FIGS_DIR / 'fig8_holdout_d1_vertical.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig8_vertical()


def make_fig8_horizontal():
    methods = fig8_data['Method']
    colors  = fig8_data['color']
    n = len(methods)
    y = np.arange(n)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    metrics = [
        ('Hits@10',
         fig8_data['Hits10_mean'], fig8_data['Hits10_lo'], fig8_data['Hits10_hi']),
        ('TPR',
         fig8_data['TPR_mean'],    fig8_data['TPR_lo'],    fig8_data['TPR_hi']),
    ]

    all_means = fig8_data['Hits10_mean'] + fig8_data['TPR_mean']
    shared_xlim = max(all_means) * 1.40

    for ax, (metric, means, lo, hi) in zip(axes, metrics):
        bars = ax.barh(y, means, height=0.55, color=colors, alpha=0.88, zorder=3,
                       xerr=[lo, hi], capsize=6,
                       ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})
        ax.set_xlim(0, shared_xlim)
        for bar, h in zip(bars, hi):
            annotate_bar_horizontal(ax, bar, h, offset_frac=0.03)
        ax.set_yticks(y)
        ax.set_yticklabels(methods, fontsize=10)
        ax.set_xlabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        style_ax(ax, 'horizontal')

    fig.suptitle(
        'Holdout Test Set — 1st Decile of Target–AE Associations\n'
        'by No. of Known AEs for Target (0–5 AEs); Best INF Config doubles TPR, Hits@10 is similar',
        fontsize=11, fontweight='bold', y=1.03
    )

    plt.tight_layout()
    out = str(FIGS_DIR / 'fig8_holdout_d1_horizontal.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig8_horizontal()


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9 — Soldatos external test set  (dual panel Hits@10 | TPR)
# ─────────────────────────────────────────────────────────────────────────────
def make_fig9_vertical():
    methods = fig9_data['Method']
    colors  = fig9_data['color']
    n = len(methods)
    x = np.arange(n)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.5), sharey=False)

    metrics = [
        ('Hits@10',
         fig9_data['Hits10_mean'], fig9_data['Hits10_lo'], fig9_data['Hits10_hi']),
        ('TPR',
         fig9_data['TPR_mean'],    fig9_data['TPR_lo'],    fig9_data['TPR_hi']),
    ]

    # Both panels share the same axis limit (0–0.5) so the Hits@10 collapse
    # is visually immediate against the INF TPR bar in the adjacent panel.
    SHARED_YLIM = 0.5

    for ax, (metric, means, lo, hi) in zip(axes, metrics):
        bars = ax.bar(x, means, width=0.55, color=colors, alpha=0.88, zorder=3,
                      yerr=[lo, hi], capsize=6,
                      ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})
        ax.set_ylim(0, SHARED_YLIM)
        fmt = "{:.2f}"
        for bar, l, h in zip(bars, lo, hi):
            annotate_bar_vertical(ax, bar, l, h, fmt=fmt, offset_frac=0.025)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontsize=10)
        ax.set_ylabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        style_ax(ax, 'vertical')

    fig.suptitle(
        'External Validation — Soldatos et al. Test Set (n = 7,804)\n'
        'Hits@10 collapses for all methods; INF TPR remains resilient',
        fontsize=11, fontweight='bold', y=1.03
    )

    legend_patches = [
        mpatches.Patch(color=COLORS[m], alpha=0.88, label=m) for m in METHOD_ORDER
    ]
    fig.legend(handles=legend_patches, loc='lower center', ncol=3,
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.07))

    plt.tight_layout()
    out = str(FIGS_DIR / 'fig9_soldatos_vertical.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig9_vertical()


def make_fig9_horizontal():
    methods = fig9_data['Method']
    colors  = fig9_data['color']
    n = len(methods)
    y = np.arange(n)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    metrics = [
        ('Hits@10',
         fig9_data['Hits10_mean'], fig9_data['Hits10_lo'], fig9_data['Hits10_hi']),
        ('TPR',
         fig9_data['TPR_mean'],    fig9_data['TPR_lo'],    fig9_data['TPR_hi']),
    ]

    SHARED_XLIM = 0.5  # shared scale for both panels

    for ax, (metric, means, lo, hi) in zip(axes, metrics):
        bars = ax.barh(y, means, height=0.55, color=colors, alpha=0.88, zorder=3,
                       xerr=[lo, hi], capsize=6,
                       ecolor='#333333', error_kw={'linewidth': 1.5, 'zorder': 4})
        ax.set_xlim(0, SHARED_XLIM)
        fmt = "{:.2f}"
        for bar, h in zip(bars, hi):
            annotate_bar_horizontal(ax, bar, h, fmt=fmt, offset_frac=0.03)
        ax.set_yticks(y)
        ax.set_yticklabels(methods, fontsize=10)
        ax.set_xlabel(metric, fontsize=11)
        ax.set_title(metric, fontsize=12, fontweight='bold')
        style_ax(ax, 'horizontal')

    fig.suptitle(
        'External Validation — Soldatos et al. Test Set (n = 7,804)\n'
        'Hits@10 collapses for all methods; INF TPR remains resilient',
        fontsize=11, fontweight='bold', y=1.03
    )

    plt.tight_layout()
    out = str(FIGS_DIR / 'fig9_soldatos_horizontal.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

make_fig9_horizontal()

print("\nAll 6 figures saved to:", str(FIGS_DIR))
