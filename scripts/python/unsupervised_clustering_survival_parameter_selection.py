# =============================================================================
# Unsupervised Clustering for Survival Parameter Selection
# =============================================================================
# This script:
#   1. Reads DE results (padj < 0.05, |log2FC| > 0.5, baseMean > 100) for each RNA type
#   2. Reads RPKM expression matrix and sample IDs
#   3. For each RNA type, selects top N variable genes (range: 10 to 500, step 50)
#   4. Performs clustering (K-means, Hierarchical, Spectral) with varying k (2-5)
#   5. Performs log-rank test to compare survival between clusters
#   6. Saves all results (p-values, cluster sizes) to an Excel file for parameter optimization
# =============================================================================

import os
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# User Configuration (relative to project root)
# =============================================================================
DATA_DIR = "data"
DE_DIR = os.path.join(DATA_DIR, "DE_results")          # DE CSV files
RPKM_FILE = os.path.join(DATA_DIR, "df_rpkm_300W.csv")
CLINICAL_FILE = os.path.join(DATA_DIR, "Meta_clinical_20241208.xlsx")

# Output directory for results
OUTPUT_DIR = os.path.join("results", "clustering_parameter_search")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# Helper Functions
# =============================================================================
def read_file(file_path):
    """Read a text file and return lines as a list (strip newlines)."""
    with open(file_path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f.readlines()]

def read_and_filter_de(file_path):
    """Read DE CSV and filter: padj < 0.05, |log2FoldChange| > 0.5, baseMean > 100."""
    data = pd.read_csv(file_path, index_col=0)
    return data.loc[(data['padj'] < 0.05) & 
                    (data['log2FoldChange'].abs() > 0.5) & 
                    (data['baseMean'] > 100), :]

# =============================================================================
# Load Data
# =============================================================================
# DE files (IPAH vs NOR for each RNA type)
de_file_paths = {
    "mRNA": os.path.join(DE_DIR, "DE_mRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "miRNA": os.path.join(DE_DIR, "DE_miRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "lncRNA": os.path.join(DE_DIR, "DE_lncRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "mttRNA": os.path.join(DE_DIR, "DE_mttRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "tRFs": os.path.join(DE_DIR, "DE_tRFs_NOR_not4pnk_IPAH_not4pnk.csv"),
    "rsRNA": os.path.join(DE_DIR, "DE_rsRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "ysRNA": os.path.join(DE_DIR, "DE_ysRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "snRNA": os.path.join(DE_DIR, "DE_snRNA_NOR_not4pnk_IPAH_not4pnk.csv"),
    "snoRNA": os.path.join(DE_DIR, "DE_snoRNA_NOR_not4pnk_IPAH_not4pnk.csv")
}

DE = {}
for name, path in de_file_paths.items():
    if os.path.exists(path):
        DE[name] = read_and_filter_de(path)
        print(f"Loaded {name}: {len(DE[name])} DE genes")
    else:
        print(f"Warning: {path} not found")

# RPKM matrix
rpkm = pd.read_csv(RPKM_FILE, index_col=0)

# Sample IDs (only PAH subtypes needed for clustering)
sampleID_files = {
    "CHD_not4pnk": os.path.join(DATA_DIR, "CHD_not4pnk_ID.txt"),
    "SLE_not4pnk": os.path.join(DATA_DIR, "SLE_not4pnk_ID.txt"),
    "NOR_not4pnk": os.path.join(DATA_DIR, "NOR_not4pnk_ID.txt"),
    "IPAH_not4pnk": os.path.join(DATA_DIR, "IPAH_not4pnk_ID.txt"),
    "IPAH_t4pnk": os.path.join(DATA_DIR, "IPAH_t4pnk_ID.txt"),
    "NOR_t4pnk": os.path.join(DATA_DIR, "NOR_t4pnk_ID.txt")
}

sampleID = {}
for key, path in sampleID_files.items():
    if os.path.exists(path):
        sampleID[key] = read_file(path)
    else:
        sampleID[key] = []

# RNA name lists (for filtering, though not strictly needed here since we use DE genes)
rna_name_files = {
    "rsRNA": os.path.join(DATA_DIR, "rsRNA.txt"),
    "ysRNA": os.path.join(DATA_DIR, "ysRNA.txt"),
    "tRFs": os.path.join(DATA_DIR, "tRFs.txt"),
    "miRNA": os.path.join(DATA_DIR, "miRNA.txt"),
    "mttRNA": os.path.join(DATA_DIR, "mttRNA.txt"),
    "snRNA": os.path.join(DATA_DIR, "snRNA.txt"),
    "snoRNA": os.path.join(DATA_DIR, "snoRNA.txt"),
    "mRNA": os.path.join(DATA_DIR, "mRNA.txt"),
    "lncRNA": os.path.join(DATA_DIR, "lncRNA.txt")
}

RNAname = {}
for key, path in rna_name_files.items():
    if os.path.exists(path):
        RNAname[key] = read_file(path)
    else:
        RNAname[key] = []

# Clinical data for survival
clinical_df = pd.read_excel(CLINICAL_FILE, sheet_name="20241223")
clinical_df = clinical_df[['Seq_ID', '生存时间2021', '随访结果2021（0=存活，1=死亡，2=失访）']].copy()
clinical_df.columns = ['Seq_ID', 'Time', 'Status']
clinical_df = clinical_df.dropna(subset=['Time'])
clinical_df['Status'] = clinical_df['Status'].astype(int)

# =============================================================================
# Parameter Grid
# =============================================================================
# RNA types to analyze (use those with DE results available)
rna_types = list(DE.keys())
# Number of top variable genes to consider (step 50 from 10 to 500)
num_genes_range = list(range(10, 501, 50))
# Number of clusters to try (2 to 5)
n_cluster_range = list(range(2, 6))

# Clustering methods
clustering_methods = {
    "K-means": lambda k: KMeans(n_clusters=k, random_state=42),
    "Hierarchical": lambda k: AgglomerativeClustering(n_clusters=k),
    "Spectral": lambda k: SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=42)
}

# =============================================================================
# Main Loop
# =============================================================================
# Dictionary to store all results
all_results = {}

for rna_type in rna_types:
    print(f"Processing RNA type: {rna_type}")
    all_results[rna_type] = {}
    
    # Subset samples: only PAH patients (IPAH, CHD, SLE) – no NOR
    sample_ids = []
    for group in ['IPAH_not4pnk', 'CHD_not4pnk', 'SLE_not4pnk']:
        sample_ids.extend(sampleID.get(group, []))
    
    # Filter RPKM matrix to these samples
    rpkm_sub = rpkm.loc[:, [col for col in rpkm.columns if col in sample_ids]]
    
    # Filter to DE genes for this RNA type
    de_genes = DE[rna_type].index.tolist()
    rpkm_de = rpkm_sub.loc[[g for g in de_genes if g in rpkm_sub.index], :]
    
    # Remove genes with zero median (optional)
    gene_medians = rpkm_de.median(axis=1)
    rpkm_de = rpkm_de.loc[gene_medians > 0, :]
    
    # Log2 transform (if not already)
    rpkm_log = np.log2(rpkm_de + 1)
    
    for num_genes in num_genes_range:
        print(f"  Top {num_genes} variable genes...")
        all_results[rna_type][num_genes] = {}
        
        # Select top variable genes
        gene_var = rpkm_log.var(axis=1)
        top_genes = gene_var.nlargest(num_genes).index
        rpkm_top = rpkm_log.loc[top_genes, :]
        
        # Transpose: samples as rows, genes as columns
        X = rpkm_top.T
        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        for n_cluster in n_cluster_range:
            print(f"    n_clusters = {n_cluster}")
            all_results[rna_type][num_genes][n_cluster] = {"p_values": {}, "cluster_sizes": {}}
            
            for method_name, model_func in clustering_methods.items():
                # Perform clustering
                model = model_func(n_cluster)
                try:
                    labels = model.fit_predict(X_scaled)
                except Exception as e:
                    print(f"      Error in {method_name}: {e}")
                    continue
                
                # Record cluster sizes
                unique, counts = np.unique(labels, return_counts=True)
                cluster_sizes = {int(label): int(count) for label, count in zip(unique, counts)}
                all_results[rna_type][num_genes][n_cluster]["cluster_sizes"][method_name] = cluster_sizes
                
                # Prepare survival data
                samples = X.index.tolist()
                cluster_df = pd.DataFrame({"Sample": samples, "Cluster": labels})
                survival_data = pd.merge(clinical_df, cluster_df, left_on='Seq_ID', right_on='Sample', how='inner')
                survival_data = survival_data.dropna(subset=['Cluster'])
                
                # Log-rank test between clusters
                clusters = survival_data['Cluster'].unique()
                if len(clusters) < 2:
                    continue
                
                p_values = {}
                for i in range(len(clusters)):
                    for j in range(i+1, len(clusters)):
                        group1 = survival_data[survival_data['Cluster'] == clusters[i]]
                        group2 = survival_data[survival_data['Cluster'] == clusters[j]]
                        if len(group1) > 0 and len(group2) > 0:
                            result = logrank_test(
                                group1['Time'], group2['Time'],
                                event_observed_A=group1['Status'], event_observed_B=group2['Status']
                            )
                            p_values[f"Cluster_{clusters[i]}_vs_Cluster_{clusters[j]}"] = result.p_value
                        else:
                            p_values[f"Cluster_{clusters[i]}_vs_Cluster_{clusters[j]}"] = np.nan
                
                all_results[rna_type][num_genes][n_cluster]["p_values"][method_name] = p_values

# =============================================================================
# Consolidate Results into a Single DataFrame
# =============================================================================
rows = []
for rna_type, dict1 in all_results.items():
    for num_genes, dict2 in dict1.items():
        for n_cluster, dict3 in dict2.items():
            for method, p_dict in dict3["p_values"].items():
                cluster_sizes = dict3["cluster_sizes"].get(method, {})
                for comp, p_val in p_dict.items():
                    row = {
                        "RNA_type": rna_type,
                        "num_genes": num_genes,
                        "n_cluster": n_cluster,
                        "Method": method,
                        "Comparison": comp,
                        "p_value": p_val
                    }
                    # Add cluster size columns (up to 10 clusters)
                    for i in range(10):
                        row[f"Cluster_{i}_size"] = cluster_sizes.get(i, 0)
                    rows.append(row)

final_df = pd.DataFrame(rows)
output_excel = os.path.join(OUTPUT_DIR, "clustering_survival_parameter_results.xlsx")
final_df.to_excel(output_excel, index=False)
print(f"Results saved to {output_excel}")

# =============================================================================
# Optional: Plot silhouette scores for optimal cluster selection
# (Can be added if needed, but not in original script)
# =============================================================================
print("Parameter search completed.")