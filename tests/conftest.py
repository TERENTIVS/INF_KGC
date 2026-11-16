"""
conftest.py — shared fixtures for reproduction tests.

Fixtures are session-scoped so each pickle is loaded once per test run.
"""

import sys
from pathlib import Path
import pytest

# Allow imports of config and helpers from the project root
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def cv_results_l3():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.CV_GS_RESULTS_L3_BOTTOM_PKL)


@pytest.fixture(scope="session")
def cv_results_l4():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.CV_GS_RESULTS_L4_PKL)


@pytest.fixture(scope="session")
def hitk_dfs():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.HITK_DFS_PKL)


@pytest.fixture(scope="session")
def hitk_partitions():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.HITK_PARTITIONS_PKL)


@pytest.fixture(scope="session")
def dwpc_negs_to_test():
    from metapaths.starterpack import load_pickle
    import config
    return load_pickle(config.DWPC_NEGS_TO_TEST_PKL)
