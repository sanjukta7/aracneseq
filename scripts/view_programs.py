
import pandas as pd
import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dataset_dir = os.path.join(base_dir, "datasets", "k562_demo")

usage_path = os.path.join(dataset_dir, "k562_demo.usages.k_10.dt_0_05.consensus.txt")
gene_scores_path = os.path.join(dataset_dir, "k562_demo.gene_spectra_score.k_10.dt_0_05.txt")

usage = pd.read_csv(usage_path, sep='\t', index_col=0)
print(f"Usage Matrix (Cells x Programs): {usage.shape}")
# Expected: (2521, 10) -> One row per cell, one col per program

gene_scores = pd.read_csv(gene_scores_path, sep='\t', index_col=0)
print(f"Gene Score Matrix (Programs x Genes): {gene_scores.shape}")
# Expected: (10, 8248) -> One row per program, one col per gene

print(gene_scores.iloc[1].sort_values(ascending=False).head(5))
