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

seed_nodes_df = pd.read_csv("CellMarkers.csv", index_col=0)
print(seed_nodes_df.head())

links = pd.read_csv("Links.csv", index_col=0)
print(links.head())

expression_data = pd.read_csv("exprSet.csv", index_col=0)
print(expression_data.iloc[1:5, 1:5])

# ============================================

G1, real_original_weights1, expression_data1, seed_nodes_df1 = miRWfunc1.prepare_data_RWR(
    seed_nodes_df=seed_nodes_df,
    links=links,
    expression_data=expression_data
)
print(f"Number of nodes: {G1.number_of_nodes()}, Number of edges: {G1.number_of_edges()}")

all_edge_S1 = []
for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
    print(f"\nProcessing Sample: {sample_id}")
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

# ============================================

G2, real_original_weights2, expression_data2 = miRWfunc2.prepare_data_RWR(
    links=links,
    expression_data=expression_data
)

print(f"Number of nodes: {G2.number_of_nodes()}, Number of edges: {G2.number_of_edges()}")

all_edge_S2 = []
for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
    print(f"\nProcessing Sample: {sample_id}")
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

print(SeedS1_edge_weights.head())
print(SeedS2_edge_weights.head())

cols = ["Sample", "link", "miRW-Imp"]
miRWImpP = SeedS1_edge_weights[cols]
print(miRWImpP.head())
miRWImpP.to_csv('miRWImpP.csv')

cols = ["Sample", "link", "miRW-Flow"]
miRWFlowP = SeedS1_edge_weights[cols]
print(miRWFlowP.head())
miRWFlowP.to_csv('miRWFlowP.csv')

cols = ["Sample", "link", "miRW-Imp"]
miRWImpE = SeedS2_edge_weights[cols]
print(miRWImpE.head())
miRWImpE.to_csv('miRWImpE.csv')

cols = ["Sample", "link", "miRW-Flow"]
miRWFlowE = SeedS2_edge_weights[cols]
print(miRWFlowE.head())
miRWFlowE.to_csv('miRWFlowE.csv')
