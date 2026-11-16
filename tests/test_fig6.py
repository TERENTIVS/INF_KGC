"""
test_fig6.py — Reproduction test for manuscript Figure 6.

Figure 6 is a bar chart comparing CV TPR on the 1st, 2nd, and 10th deciles
for L3 and L4 original (untransformed) metapath count features.

The six values that feed the plot call (notebook cell 218) are the ground
truth for this test.  They are derived deterministically from the CV results
pickles via analyse_parameter_slices_quantile() — no randomness is involved.

Pass criterion: all six values match to 3 decimal places.
"""

import pytest
from tests.helpers import analyse_parameter_slices_quantile

# Ground truth: hardcoded literals from notebook cell 218
# plot_tpr_comparison(
#     [0.202, 0.56,  0.959],   <- L3: deciles 1, 2, 10
#     [0.223, 0.575, 0.966],   <- L4: deciles 1, 2, 10
#     ...
# )
EXPECTED = {
    "l3": {0: 0.202, 1: 0.560, 9: 0.959},
    "l4": {0: 0.223, 1: 0.575, 9: 0.966},
}

TOLERANCE = 0.001   # 3 decimal places


# ---------------------------------------------------------------------------
# L3
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def l3_slices(cv_results_l3):
    # config 0 = raw (untransformed) counts; index 1 = per-fold CV results
    return analyse_parameter_slices_quantile(cv_results_l3[0][1])


def test_l3_decile1_tpr(l3_slices):
    """L3 raw counts: 1st-decile (0–5 AEs) mean CV TPR ≈ 0.202."""
    actual = l3_slices.iloc[0]["recall_mean"]
    assert abs(actual - EXPECTED["l3"][0]) <= TOLERANCE, (
        f"L3 decile 1: expected {EXPECTED['l3'][0]}, got {actual}"
    )


def test_l3_decile2_tpr(l3_slices):
    """L3 raw counts: 2nd-decile mean CV TPR ≈ 0.560."""
    actual = l3_slices.iloc[1]["recall_mean"]
    assert abs(actual - EXPECTED["l3"][1]) <= TOLERANCE, (
        f"L3 decile 2: expected {EXPECTED['l3'][1]}, got {actual}"
    )


def test_l3_decile10_tpr(l3_slices):
    """L3 raw counts: 10th-decile (163–256 AEs) mean CV TPR ≈ 0.959."""
    actual = l3_slices.iloc[9]["recall_mean"]
    assert abs(actual - EXPECTED["l3"][9]) <= TOLERANCE, (
        f"L3 decile 10: expected {EXPECTED['l3'][9]}, got {actual}"
    )


# ---------------------------------------------------------------------------
# L4
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def l4_slices(cv_results_l4):
    return analyse_parameter_slices_quantile(cv_results_l4[0][1])


def test_l4_decile1_tpr(l4_slices):
    """L4 raw counts: 1st-decile mean CV TPR ≈ 0.223."""
    actual = l4_slices.iloc[0]["recall_mean"]
    assert abs(actual - EXPECTED["l4"][0]) <= TOLERANCE, (
        f"L4 decile 1: expected {EXPECTED['l4'][0]}, got {actual}"
    )


def test_l4_decile2_tpr(l4_slices):
    """L4 raw counts: 2nd-decile mean CV TPR ≈ 0.575."""
    actual = l4_slices.iloc[1]["recall_mean"]
    assert abs(actual - EXPECTED["l4"][1]) <= TOLERANCE, (
        f"L4 decile 2: expected {EXPECTED['l4'][1]}, got {actual}"
    )


def test_l4_decile10_tpr(l4_slices):
    """L4 raw counts: 10th-decile mean CV TPR ≈ 0.966."""
    actual = l4_slices.iloc[9]["recall_mean"]
    assert abs(actual - EXPECTED["l4"][9]) <= TOLERANCE, (
        f"L4 decile 10: expected {EXPECTED['l4'][9]}, got {actual}"
    )
