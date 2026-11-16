"""
generate_fig4_data.py
=====================
Queries the TSA graph for the top 20 Gene nodes by number of Gene_ADR edges,
maps Entrez Gene IDs to common names, and writes the result to:

    data/fig4_top20_targets.csv

This CSV is included in the data deposit so that Figure 4 can be reproduced
without a live Neo4j connection.  Run once, with the read database running:

    python generate_fig4_data.py

Requires the full TSA graph (all Gene_ADR edges present, i.e. test edges
restored).  Expected: 6,657 Gene_ADR edges total.
"""

import pandas as pd
from py2neo import Graph
from pathlib import Path
import config

# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------

QUERY = """
MATCH (n:Gene)-[r:Gene_ADR]-()
RETURN n.name AS Gene, count(r) AS AEs
ORDER BY AEs DESC
LIMIT 20
"""

print(f"Connecting to {config.NEO4J_URI_READ} ...")
graph = Graph(config.NEO4J_URI_READ)

# Verify this is the full graph (all edges present)
total_gene_adr = graph.run(
    "MATCH ()-[r:Gene_ADR]->() RETURN COUNT(r) AS n"
).evaluate()
assert total_gene_adr == 6657, (
    f"Expected 6,657 Gene_ADR edges (full graph), found {total_gene_adr}. "
    "Ensure all test edges are restored before running."
)
print(f"  Gene_ADR edge count verified: {total_gene_adr}")

print("Running top-20 query ...")
top_targets = pd.DataFrame(graph.run(QUERY).data())
print(f"  {len(top_targets)} rows returned.")

# ---------------------------------------------------------------------------
# Map Entrez Gene IDs to common names via the hetionet node mapping file
# ---------------------------------------------------------------------------

print(f"Loading name mapping from {config.MAPPING_FILE} ...")
map_df = pd.read_csv(config.MAPPING_FILE, sep="\t")

top_targets = (
    top_targets
    .merge(map_df[["id", "name"]], left_on="Gene", right_on="id", how="left")
    .drop(columns=["id"])
    .rename(columns={"name": "Gene_name"})
    [["Gene", "Gene_name", "AEs"]]
)

# ---------------------------------------------------------------------------
# Sanity checks (values cited in the paper)
# ---------------------------------------------------------------------------

TOTAL_GT = 6657
top20_sum = int(top_targets["AEs"].sum())
top20_pct = top20_sum / TOTAL_GT

print(f"\nSanity checks:")
print(f"  Top-20 AE sum : {top20_sum}  (paper states 2,997)")
print(f"  Top-20 share  : {top20_pct:.1%}  (paper states ~45%)")

assert top20_sum == 2997, f"Expected 2,997, got {top20_sum}"
assert 0.44 <= top20_pct <= 0.46, f"Expected ~45%, got {top20_pct:.1%}"
print("  All checks passed.")

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------

out_path = Path(__file__).parent / "data" / "fig4_top20_targets.csv"
out_path.parent.mkdir(exist_ok=True)
top_targets.to_csv(out_path, index=False)
print(f"\nSaved to {out_path}")
print(top_targets.to_string(index=False))
