"""
test_si_inf_grid.py — Reproduction tests for the §3.2 SI artefacts
produced by ``generate_si_inf_grid.py``.

All three claims below are checked under the same fold-mean estimator
(``analyse_parameter_slices_quantile``) that drives the headline §3.2
numbers and Figure 7 of the main text, so the SI sits on a single
internally consistent footing.

1. Decile-1 best INF config is config 12 with mean CV TPR = 0.424
   (matches ``best_configs_l3_recall.pkl``, the §3.2 headline number).
2. The top-3 INF configurations at decile 1 under the fold-mean
   estimator all use deflator exponent = 0 — supports the worked
   example proposed for the §3.2 rewrite.
3. The best single global INF config (by mean-across-deciles CV TPR,
   excluding raw baseline) is config 20 — i.e. deflator_exp=0.3,
   agg=product, pool=mean — matching the description in §3.2.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config
from metapaths_src.metapaths.starterpack import load_pickle
from tests.helpers import analyse_parameter_slices_quantile


# ── Fixtures ───────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def best_configs():
    return load_pickle(config.BEST_CONFIGS_L3_RECALL_PKL).sort_values(
        "slice_num").reset_index(drop=True)


@pytest.fixture(scope="module")
def cv_l3():
    return load_pickle(config.CV_GS_RESULTS_L3_PKL)


@pytest.fixture(scope="module")
def param_combos():
    return load_pickle(config.PARAM_COMBOS_PKL)


@pytest.fixture(scope="module")
def per_config_per_decile(cv_l3):
    """(35, 10) matrix of per-decile mean CV TPR via the fold-mean estimator."""
    import numpy as np
    matrix = np.full((35, 10), np.nan)
    for cfg_idx in range(35):
        df = analyse_parameter_slices_quantile(cv_l3[cfg_idx][1])
        for _, row in df.iterrows():
            matrix[cfg_idx, int(row["slice_num"]) - 1] = row["recall_mean"]
    return matrix


# ── Claim 1: decile-1 headline number ──────────────────────────────────────
def test_decile_1_best_config_is_12(best_configs):
    row = best_configs.iloc[0]
    assert int(row["slice_num"]) == 1
    assert int(row["param_combo"]) == 12
    assert row["recall_mean"] == pytest.approx(0.424, abs=1e-3)


def test_decile_1_topscore_matches_estimator(per_config_per_decile):
    """Config 12 D1 from the live estimator must match the cached 0.424."""
    assert per_config_per_decile[12, 0] == pytest.approx(0.424, abs=1e-3)


# ── Claim 2: decile-1 top-3 all use deflator_exp = 0 ──────────────────────
def test_decile_1_top3_all_deflator_zero(per_config_per_decile, param_combos):
    import numpy as np
    order = np.argsort(-per_config_per_decile[:, 0])
    top3 = [int(i) for i in order[:3]]

    for cfg_idx in top3:
        assert param_combos[cfg_idx]["path_deflator_exp"] == 0, (
            f"config {cfg_idx} in decile-1 top-3 does not have "
            f"deflator_exp=0: {param_combos[cfg_idx]}"
        )


# ── Claim 3: best single global INF config is config 20 ───────────────────
def test_global_best_config_is_20(per_config_per_decile, param_combos):
    import numpy as np
    overall = np.nanmean(per_config_per_decile, axis=1)
    # Exclude raw baseline (config 0)
    best_cfg = 1 + int(np.argmax(overall[1:]))
    best_score = float(overall[best_cfg])

    assert best_cfg == 20
    assert param_combos[best_cfg] == {
        "path_deflator_exp": 0.3,
        "inf_inflator":      "product",
        "inf_pooling":       "mean",
    }
    # Mean across deciles for config 20 under fold-mean estimator
    assert best_score == pytest.approx(0.872, abs=2e-3)
