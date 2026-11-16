"""
helpers.py — analysis functions extracted from INF_notebook.ipynb.

These are verbatim copies of the notebook implementations, included here
so the test suite can run them independently of the notebook kernel.
All function signatures and logic are unchanged.
"""

import numpy as np
import pandas as pd
from scipy.stats import t
from sklearn.metrics import recall_score, precision_score, f1_score
from sklearn.inspection import permutation_importance
from sklearn.metrics import make_scorer
from typing import Callable, Dict, List


# ---------------------------------------------------------------------------
# Figure 6 helpers  (notebook cells 44–45)
# ---------------------------------------------------------------------------

def ci_margin_of_error(data: pd.Series, confidence: float = 0.95) -> float:
    if len(data) < 2:
        return 0.0
    std_dev = np.std(data, ddof=1)
    n = len(data)
    t_score = t.ppf((1 + confidence) / 2, df=n - 1)
    return t_score * (std_dev / np.sqrt(n))


def analyse_parameter_slices_quantile(
    fold_results: dict,
    num_slices: int = 10,
    balance_slices: bool = False,
    confidence_level: float = 0.95,
):
    slice_data = []
    RANDOM_SEED = 69

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
            pos_df['slice_num'] = pd.qcut(
                pos_df['param_values'], q=num_slices,
                labels=False, duplicates='drop',
            )
        except ValueError:
            continue

        for slice_id in sorted(pos_df['slice_num'].unique()):
            slice_pos_df = pos_df[pos_df['slice_num'] == slice_id]
            min_param = slice_pos_df['param_values'].min()
            max_param = slice_pos_df['param_values'].max()
            n_pos = len(slice_pos_df)

            slice_neg_conformant = fold_df[
                (fold_df['y_val'] == 0) &
                (fold_df['param_values'].between(min_param, max_param))
            ]
            n_val_neg = len(slice_neg_conformant)

            if balance_slices and n_val_neg > n_pos:
                final_slice_negs = slice_neg_conformant.sample(
                    n=n_pos, random_state=RANDOM_SEED
                )
            else:
                final_slice_negs = slice_neg_conformant

            final_slice_df = pd.concat([slice_pos_df, final_slice_negs])
            y_true_slice = final_slice_df['y_val'].values
            y_pred_slice = final_slice_df['y_pred'].values

            recall    = recall_score(y_true_slice, y_pred_slice, zero_division=0)
            precision = precision_score(y_true_slice, y_pred_slice, zero_division=0)
            f_1       = f1_score(y_true_slice, y_pred_slice, zero_division=0)
            n_true_pos = np.sum((y_pred_slice == 1) & (y_true_slice == 1))
            n_pred_pos = np.sum(y_pred_slice)

            slice_data.append({
                'slice_num': slice_id + 1,
                'recall':    recall,
                'precision': precision,
                'f_1':       f_1,
                'min_param': min_param,
                'max_param': max_param,
                'true_pos':  n_true_pos,
                'pred_pos':  n_pred_pos,
                'val_neg':   len(final_slice_negs),
                'pos':       n_pos,
            })

    if not slice_data:
        return pd.DataFrame()

    df = pd.DataFrame(slice_data)
    ci_func = lambda x: ci_margin_of_error(x, confidence=confidence_level)

    result = df.groupby('slice_num').agg(
        recall_mean        =('recall',    'mean'),
        recall_ci          =('recall',    ci_func),
        precision_mean     =('precision', 'mean'),
        precision_ci       =('precision', ci_func),
        f1_mean            =('f_1',       'mean'),
        f1_ci              =('f_1',       ci_func),
        mean_true_pos      =('true_pos',  'mean'),
        mean_pred_pos      =('pred_pos',  'mean'),
        mean_val_neg       =('val_neg',   'mean'),
        mean_pos_per_slice =('pos',       'mean'),
        min_param          =('min_param', 'mean'),
        max_param          =('max_param', 'mean'),
    ).reset_index()

    result[['min_param', 'max_param']] = (
        result[['min_param', 'max_param']].round(0).astype(int)
    )
    return result.round(3)


# ---------------------------------------------------------------------------
# Feature-importance helpers (Figures 10–11; notebook cells 496–497)
# ---------------------------------------------------------------------------

def permimp_ranking(model, X, y_true, feat_names, scorer):
    permimp = permutation_importance(
        model, X, y_true,
        scoring=scorer,
        n_repeats=5,
        random_state=69,
    )['importances_mean']

    permimp_dict = {i: permimp[i] for i in range(len(permimp))}

    return {
        feat_names[feat_index]: round(permimp_dict[feat_index], 3)
        for feat_index in sorted(
            permimp_dict, key=lambda k: permimp_dict[k], reverse=True
        )[:10]
    }


def get_sliced_permutation_importance(
    model: Callable,
    eval_df: pd.DataFrame,
    neg_to_test: pd.DataFrame,
    boundaries_df: pd.DataFrame,
    feature_cols: List[str],
    param_col: str,
    head_col: str,
    tail_col: str,
    k: int = 10,
    optimization_metric: str = 'recall',
    top_n_features: int = 5,
) -> Dict[int, pd.DataFrame]:

    best_configs = boundaries_df.copy()
    for i in range(1, len(best_configs)):
        best_configs.loc[i, 'min_param'] = best_configs.loc[i - 1, 'max_param'] + 1

    full_eval_data = pd.merge(
        eval_df, neg_to_test,
        on=[head_col, tail_col, param_col],
        how='inner',
    )
    if full_eval_data.empty:
        return {}

    def recall_scorer(y_true, y_pred, **kwargs):
        return recall_score(y_true, (y_pred > 0.5).astype(int))

    if optimization_metric == 'recall':
        scorer = make_scorer(recall_scorer, needs_proba=True)
    else:
        raise ValueError("optimization_metric must be 'recall'")

    results = {}
    for slice_idx, row in best_configs.iterrows():
        min_b, max_b = row['min_param'], row['max_param']
        mask = (
            (full_eval_data[param_col] >= min_b) &
            (full_eval_data[param_col] <= max_b)
        )
        data_slice = full_eval_data.loc[mask].copy()
        if data_slice.empty:
            continue

        X_slice = data_slice[feature_cols]
        y_slice = data_slice['true_label']

        importance_ranks = permimp_ranking(
            model=model,
            X=X_slice,
            y_true=y_slice,
            feat_names=feature_cols,
            scorer=scorer,
        )

        importance_df = (
            pd.DataFrame(importance_ranks.items(), columns=['feature', 'importance'])
            .sort_values('importance', ascending=False)
            .head(top_n_features)
            .reset_index(drop=True)
        )
        results[slice_idx] = importance_df

    return results
