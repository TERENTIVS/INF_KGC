"""
test_figure_data.py — Checks the bundled aggregate values for figures whose
inputs are restricted (Figure 8, Figure 11, Supplementary Figure S4).

Values are those reported in the manuscript (two decimals) and in
data/figure_data/; tolerance ±0.002 allows for the unseeded bootstrap /
ensemble retraining described in data/figure_data/README.md.
"""

import pandas as pd
import pytest

import config

TOL = 0.002
GCG_ADR = ("(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)"
           "<-[r2:Gene_Compound]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)")


@pytest.fixture(scope="module")
def fig8():
    return pd.read_csv(config.FIGURE_DATA_DIR / "fig8_holdout_d1.csv").set_index(["method", "metric"])


@pytest.fixture(scope="module")
def figs4():
    return pd.read_csv(config.FIGURE_DATA_DIR / "figS4_holdout_vs_external.csv").set_index(
        ["test_set", "method", "metric"])


@pytest.fixture(scope="module")
def fig11():
    return pd.read_csv(config.FIGURE_DATA_DIR / "fig11_external_feat_importance_d1.csv")


@pytest.mark.parametrize("method,metric,expected", [
    ("Best INF Config", "TPR", 0.476), ("DWPC", "TPR", 0.256), ("Original", "TPR", 0.201),
    ("Best INF Config", "Hits@10", 0.376), ("DWPC", "Hits@10", 0.376), ("Original", "Hits@10", 0.321),
])
def test_fig8_values(fig8, method, metric, expected):
    """Figure 8: holdout 1st-decile TPR (0.48 / 0.26 / 0.20) and Hits@10."""
    row = fig8.loc[(method, metric)]
    assert abs(row["value"] - expected) <= TOL
    assert row["ci_low"] <= row["value"] <= row["ci_high"]


@pytest.mark.parametrize("test_set,method,metric,expected", [
    ("holdout", "INF Ensemble", "TPR", 0.90), ("holdout", "DWPC", "TPR", 0.868),
    ("holdout", "Original", "TPR", 0.847),
    ("external", "INF Ensemble", "TPR", 0.402), ("external", "DWPC", "TPR", 0.215),
    ("external", "Original", "TPR", 0.217),
    ("holdout", "INF Ensemble", "Hits@10", 0.353), ("external", "INF Ensemble", "Hits@10", 0.064),
])
def test_figS4_values(figs4, test_set, method, metric, expected):
    """Figure S4 / §3.3: overall TPR 0.90 vs 0.87 / 0.85 (holdout), 0.40 vs 0.22 (external)."""
    tol = 0.005 if (test_set, method) == ("holdout", "INF Ensemble") else TOL
    assert abs(figs4.loc[(test_set, method, metric), "value"] - expected) <= tol


def test_fig11_inf_rank1(fig11):
    """Figure 11: rank-1 INF feature on external data is Gene_Compound_Gene_ADR (0.104)."""
    top = fig11[(fig11["model"] == "INF (config 12)") & (fig11["rank"] == 1)].iloc[0]
    assert top["feature"] == GCG_ADR
    assert abs(top["importance"] - 0.104) <= 0.0005


def test_fig11_shape(fig11):
    """Figure 11: five ranked features for each of the three models."""
    assert fig11.groupby("model")["rank"].apply(list).apply(lambda r: r == [1, 2, 3, 4, 5]).all()
    assert fig11["model"].nunique() == 3
