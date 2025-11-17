import os
import time
import glob
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from scipy.sparse import csr_matrix
import shutil
import sys
import miRWfunc1
import miRWfunc2

# ============================================
#            Read input data
# ============================================

# network links
links = pd.read_csv("data_example/Links.csv", index_col=0)

# expression matrix
expression_data = pd.read_csv("data_example/exprSet.csv", index_col=0)

# CellMarkers.csv is optional
cellmarker_path = "data_example/CellMarkers.csv"
has_seed_markers = os.path.exists(cellmarker_path)

if has_seed_markers:
    print(f"\nFound {cellmarker_path}, running both S1 (seed-based) and S2 (seed-free)...")
    seed_nodes_df = pd.read_csv(cellmarker_path, index_col=0)
    print("Seed nodes (CellMarkers):")
    print(seed_nodes_df.head())
else:
    print(f"\n{cellmarker_path} not found. Skipping S1, only running S2 (seed-free).")

# ============================================
#                 S1: with seeds
# ============================================

if has_seed_markers:
    G1, real_original_weights1, expression_data1, seed_nodes_df1 = miRWfunc1.prepare_data_RWR(
        seed_nodes_df=seed_nodes_df,
        links=links,
        expression_data=expression_data
    )
    print(f"[S1] Number of nodes: {G1.number_of_nodes()}, Number of edges: {G1.number_of_edges()}")

    all_edge_S1 = []
    for sample_id in tqdm(expression_data.columns, desc="Processing Samples S1", unit="sample"):
        print(f"\n[S1] Processing Sample: {sample_id}")
        edge_weights_df = miRWfunc1.RWR_S1(
            sample_id=sample_id,
            G=G1,
            real_original_weights=real_original_weights1,
            expression_data=expression_data1,
            seed_nodes_df=seed_nodes_df1
        )
        print(edge_weights_df.head())
        print(edge_weights_df.shape)
        all_edge_S1.append(edge_weights_df)

    SeedS1_edge_weights = pd.concat(all_edge_S1, ignore_index=True)
    print("\n[S1] Combined edge weights (head):")
    print(SeedS1_edge_weights.head())
else:
    SeedS1_edge_weights = None  # 占位，方便后面判断

# ============================================
#                 S2: no seeds
# ============================================

G2, real_original_weights2, expression_data2 = miRWfunc2.prepare_data_RWR(
    links=links,
    expression_data=expression_data
)

print(f"\n[S2] Number of nodes: {G2.number_of_nodes()}, Number of edges: {G2.number_of_edges()}")

all_edge_S2 = []
for sample_id in tqdm(expression_data.columns, desc="Processing Samples S2", unit="sample"):
    print(f"\n[S2] Processing Sample: {sample_id}")
    edge_weights_df = miRWfunc2.RWR_S2(
        sample_id=sample_id,
        G=G2,
        real_original_weights=real_original_weights2,
        expression_data=expression_data2
    )
    print(edge_weights_df.head())
    print(edge_weights_df.shape)
    all_edge_S2.append(edge_weights_df)

SeedS2_edge_weights = pd.concat(all_edge_S2, ignore_index=True)
print("\n[S2] Combined edge weights (head):")
print(SeedS2_edge_weights.head())

# ============================================
#               Save results
# ============================================

# S1 outputs (only if CellMarkers.csv exists)
if SeedS1_edge_weights is not None:
    cols_imp = ["Sample", "link", "miRW-Imp"]
    cols_flow = ["Sample", "link", "miRW-Flow"]

    miRWImpP = SeedS1_edge_weights[cols_imp]
    print("\n[S1] miRWImpP (head):")
    print(miRWImpP.head())
    miRWImpP.to_csv('miRWImpP.csv', index=False)

    miRWFlowP = SeedS1_edge_weights[cols_flow]
    print("\n[S1] miRWFlowP (head):")
    print(miRWFlowP.head())
    miRWFlowP.to_csv('miRWFlowP.csv', index=False)
else:
    print("\n[S1] Skipped (no CellMarkers.csv). No miRWImpP.csv / miRWFlowP.csv generated.")

# S2 outputs (always)
cols_imp = ["Sample", "link", "miRW-Imp"]
cols_flow = ["Sample", "link", "miRW-Flow"]

miRWImpE = SeedS2_edge_weights[cols_imp]
print("\n[S2] miRWImpE (head):")
print(miRWImpE.head())
miRWImpE.to_csv('miRWImpE.csv', index=False)

miRWFlowE = SeedS2_edge_weights[cols_flow]
print("\n[S2] miRWFlowE (head):")
print(miRWFlowE.head())
miRWFlowE.to_csv('miRWFlowE.csv', index=False)

print("\nDone.")
