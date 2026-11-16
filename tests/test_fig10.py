"""
test_fig10.py — Reproduction test for manuscript Figure 10.

Figure 10 reports the top-5 metapath features by permutation importance for
TPR (recall) on the holdout test data, 1st decile (targets with 0–5 known AEs),
for the best INF configuration ("INF_CV_Best", config 12) and, for comparison,
the raw-count (config 0) and DWPC benchmarks.

Ground truth is the bundled pickle config_top_feats_by_slice.pkl
(config.FIG10_FEAT_IMP_PKL), produced by the "Top features" section of
INF_notebook_clean.ipynb: dict[config][metric][slice] -> DataFrame; slice 0
is the 1st decile.

Pass criteria:
  1. Config 12, rank-1 feature is Gene_Compound_Gene_ADR — the feature
     discussed in §3.5 as most informative for low-degree targets.
  2. Top-5 rankings for configs 12, 0 and DWPC match the values shown in
     Figure 10 (rank order + importance to 3 d.p.).
"""

import pytest

BEST_INF_RANK1_FEATURE = (
    "(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)"
    "<-[r2:Gene_Compound]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)"
)

EXPECTED = {
    12: [
        ('(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)<-[r2:Gene_Compound]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.13),
        ('(n_source:Gene)-[r1:Gene_Pathway]->(n_1:Pathway)<-[r2:Gene_Pathway]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.081),
        ('(n_source:Gene)-[r1:Gene_MolecularFunction]->(n_1:MolecularFunction)<-[r2:Gene_MolecularFunction]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.079),
        ('(n_source:Gene)-[r1:Gene_BiologicalProcess]->(n_1:BiologicalProcess)<-[r2:Gene_BiologicalProcess]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.079),
        ('(n_source:Gene)-[r1:Gene_Gene]->(n_1:Gene)-[r2:Gene_ADR]->(n_target:ADR)', 0.059),
    ],
    0: [
        ('(n_source:Gene)-[r1:Gene_MolecularFunction]->(n_1:MolecularFunction)<-[r2:Gene_MolecularFunction]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.068),
        ('(n_source:Gene)-[r1:Gene_Gene]->(n_1:Gene)-[r2:Gene_ADR]->(n_target:ADR)', 0.055),
        ('(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)<-[r2:Gene_Compound]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.048),
        ('(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)-[r2:Compound_SideEffect]->(n_2:SideEffect)<-[r3:ADR_SideEffect]-(n_target:ADR)', 0.026),
        ('(n_source:Gene)-[r1:Gene_BiologicalProcess]->(n_1:BiologicalProcess)<-[r2:Gene_BiologicalProcess]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.02),
    ],
    'dwpc': [
        ('(n_source:Gene)-[r1:Gene_Pathway]->(n_1:Pathway)<-[r2:Gene_Pathway]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.057),
        ('(n_source:Gene)-[r1:Gene_MolecularFunction]->(n_1:MolecularFunction)<-[r2:Gene_MolecularFunction]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.037),
        ('(n_source:Gene)-[r1:Gene_Gene]->(n_1:Gene)-[r2:Gene_ADR]->(n_target:ADR)', 0.031),
        ('(n_source:Gene)-[r1:Gene_Compound]->(n_1:Compound)<-[r2:Gene_Compound]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.028),
        ('(n_source:Gene)-[r1:Gene_ADR]->(n_1:ADR)<-[r2:Gene_ADR]-(n_2:Gene)-[r3:Gene_ADR]->(n_target:ADR)', 0.026),
    ],
}

TOLERANCE = 0.001


@pytest.fixture(scope="module")
def feat_imp_fig10():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.FIG10_FEAT_IMP_PKL)


def _slice0(feat_imp, config_key):
    return feat_imp[config_key]["recall"][0]


def test_best_inf_rank1_feature(feat_imp_fig10):
    """Config 12: rank-1 feature is Gene_Compound_Gene_ADR (§3.5)."""
    actual = _slice0(feat_imp_fig10, 12).iloc[0]["feature"]
    assert actual == BEST_INF_RANK1_FEATURE, f"Config 12 rank-1 feature: got {actual}"


@pytest.mark.parametrize("config_key", [12, 0, "dwpc"])
def test_top5_features(feat_imp_fig10, config_key):
    """Top-5 ranking matches Figure 10 for INF (12), raw (0) and DWPC."""
    df = _slice0(feat_imp_fig10, config_key)
    for rank, (expected_feat, expected_imp) in enumerate(EXPECTED[config_key]):
        actual_feat = df.iloc[rank]["feature"]
        actual_imp = df.iloc[rank]["importance"]
        assert actual_feat == expected_feat, (
            f"Config {config_key} rank {rank+1}: feature mismatch\n"
            f"  expected: {expected_feat}\n  got:      {actual_feat}"
        )
        assert abs(actual_imp - expected_imp) <= TOLERANCE, (
            f"Config {config_key} rank {rank+1}: importance {actual_imp} != {expected_imp}"
        )
