# =============================================================================
# Unsupervised Clustering, Survival Analysis, and Clinical Correlation
# =============================================================================
# This script:
#   1. Reads DE results (padj < 0.05, |log2FC| > 0.5) for each RNA type
#   2. Reads RPKM expression matrix and sample IDs
#   3. Selects top 50 most variable genes across all RNA types
#   4. Performs clustering (K-means, Hierarchical, Spectral) with optimal k (determined by silhouette score)
#   5. Performs Kaplan-Meier survival analysis and log-rank tests
#   6. Correlates clusters with clinical parameters (continuous variables: boxplots; categorical: stacked barplots)
#   7. Generates PDF reports for each clustering method
# =============================================================================

import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# User Configuration (relative to project root)
# =============================================================================
DATA_DIR = "data"
DE_DIR = os.path.join(DATA_DIR, "DE_results")
RPKM_FILE = os.path.join(DATA_DIR, "df_rpkm_300W.csv")
CLINICAL_FILE = os.path.join(DATA_DIR, "Meta_clinical_20241208.xlsx")
REVEAL_FILE = os.path.join(DATA_DIR, "Meta_REVEAL_241223.csv")

# Output directories
OUTPUT_RESULTS = os.path.join("results", "clustering")
OUTPUT_FIGURES = os.path.join("figures", "clustering")
os.makedirs(OUTPUT_RESULTS, exist_ok=True)
os.makedirs(OUTPUT_FIGURES, exist_ok=True)

# =============================================================================
# Helper Functions
# =============================================================================
def read_file(file_path):
    """Read a text file and return lines as a list (strip newlines)."""
    with open(file_path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f.readlines()]

def read_and_filter_de(file_path):
    """Read DE CSV and filter: padj < 0.05, |log2FoldChange| > 0.5."""
    data = pd.read_csv(file_path, index_col=0)
    return data.loc[(data['padj'] < 0.05) & (data['log2FoldChange'].abs() > 0.5), :]

def find_optimal_clusters(data, k_min=2, k_max=10, random_state=42):
    """
    Automatically select optimal number of clusters (k) for K-means, Hierarchical, and Spectral.
    Returns best_k and corresponding method.
    """
    methods = {
        "K-means": lambda k: KMeans(n_clusters=k, random_state=random_state),
        "Hierarchical": lambda k: AgglomerativeClustering(n_clusters=k),
        "Spectral": lambda k: SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=random_state)
    }
    results = []
    for name, model_func in methods.items():
        for k in range(k_min, k_max+1):
            try:
                model = model_func(k)
                labels = model.fit_predict(data)
                score = silhouette_score(data, labels)
                results.append({"Method": name, "k": k, "Silhouette": score})
            except Exception as e:
                print(f"{name}, k={k} failed: {e}")
    df = pd.DataFrame(results)
    best = df.loc[df["Silhouette"].idxmax()]
    return best["Method"], int(best["k"]), df

def plot_silhouette_curve(df, output_pdf):
    """Plot silhouette scores for different methods and k."""
    plt.figure(figsize=(4.5/2.54, 4.5/2.54*0.75))
    for method in df["Method"].unique():
        subset = df[df["Method"] == method]
        plt.plot(subset["k"], subset["Silhouette"], marker='o', label=method)
    plt.xlabel("Number of clusters (k)")
    plt.ylabel("Silhouette Score")
    plt.title("Optimal number of clusters")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300)
    plt.close()

def perform_survival_analysis(cluster_df, ax):
    """Perform Kaplan-Meier and log-rank test for clusters."""
    survival_data = pd.merge(clinical_surv, cluster_df, left_on='Seq_ID', right_on='Sample', how='left')
    survival_data = survival_data.dropna(subset=['Cluster'])
    clusters = survival_data['Cluster'].unique()
    kmfs = []
    for cl in sorted(clusters):
        sub = survival_data[survival_data['Cluster'] == cl]
        kmf = KaplanMeierFitter()
        kmf.fit(sub['Time'], event_observed=sub['Status'], label=f"Cluster {cl} (N={len(sub)})")
        kmf.plot(ax=ax)
        kmfs.append(kmf)
    # Log-rank test
    if len(clusters) >= 2:
        p_values = {}
        for i in range(len(clusters)):
            for j in range(i+1, len(clusters)):
                cl1 = clusters[i]; cl2 = clusters[j]
                group1 = survival_data[survival_data['Cluster'] == cl1]
                group2 = survival_data[survival_data['Cluster'] == cl2]
                result = logrank_test(group1['Time'], group2['Time'],
                                      event_observed_A=group1['Status'], event_observed_B=group2['Status'])
                p_values[f"{cl1} vs {cl2}"] = result.p_value
        # Add p-value annotation
        p_text = "\n".join([f"{k}: p={v:.4f}" for k, v in p_values.items()])
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    # Add risk table
    add_at_risk_counts(*kmfs, ax=ax, rows_to_show=['At risk'], xticks=[0,12,24,36,48,60])
    ax.set_title("Kaplan-Meier Survival Curves")
    ax.set_xlabel("Time (Months)")
    ax.set_ylabel("Survival Probability")

def plot_box_with_pval(data, column, cluster_col='Cluster', ax=None):
    """Boxplot for continuous clinical variable with t-test p-value."""
    data_filtered = data.dropna(subset=[column, cluster_col])
    if len(data_filtered[cluster_col].unique()) != 2:
        return
    clusters = sorted(data_filtered[cluster_col].unique())
    group1 = data_filtered[data_filtered[cluster_col] == clusters[0]][column]
    group2 = data_filtered[data_filtered[cluster_col] == clusters[1]][column]
    t_stat, p_val = ttest_ind(group1, group2, equal_var=False)
    if ax is None:
        fig, ax = plt.subplots(figsize=(4*0.3937, 8*0.3937))
    sns.boxplot(x=cluster_col, y=column, data=data_filtered, palette="Set2", ax=ax)
    ax.text(0.5, 0.95, f"p = {p_val:.6f}", transform=ax.transAxes, ha='center', va='top')
    ax.set_title(column)
    ax.set_xlabel("Cluster")
    ax.set_ylabel(column)
    return ax.figure

def plot_stacked_percentage(data, column, cluster_col='Cluster', ax=None):
    """Stacked barplot for categorical clinical variable."""
    data_filtered = data.dropna(subset=[column, cluster_col])
    if len(data_filtered) == 0:
        return
    cross_tab = pd.crosstab(data_filtered[cluster_col], data_filtered[column], normalize='index')
    if ax is None:
        fig, ax = plt.subplots(figsize=(4*0.3937, 7.2*0.3937))
    cross_tab.plot(kind='bar', stacked=True, ax=ax, edgecolor='black', linewidth=0.4)
    ax.set_title(column)
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Cluster")
    ax.set_ylim(0, 1.2)
    # Add sample counts
    for i, cluster in enumerate(sorted(data_filtered[cluster_col].unique())):
        n = len(data_filtered[data_filtered[cluster_col] == cluster])
        ax.text(i, -0.1, f"N={n}", ha='center')
    ax.legend(fontsize=7)
    return ax.figure

# =============================================================================
# Load Data
# =============================================================================
# DE files (IPAH vs NOR for each RNA type)
de_files = {
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
for name, path in de_files.items():
    if os.path.exists(path):
        DE[name] = read_and_filter_de(path)
        print(f"Loaded {name}: {len(DE[name])} DE genes")
    else:
        print(f"Warning: {path} not found")

# RPKM matrix
rpkm = pd.read_csv(RPKM_FILE, index_col=0)

# Sample IDs (only PAH patients: IPAH, CHD, SLE)
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

# RNA name lists (not strictly needed for DE-based selection but used for filtering)
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

# =============================================================================
# Prepare Expression Matrix for Clustering (all PAH samples, top 50 variable genes across all RNA types)
# =============================================================================
# Subset PAH samples only (IPAH, CHD, SLE)
pah_samples = sampleID["IPAH_not4pnk"] + sampleID["CHD_not4pnk"] + sampleID["SLE_not4pnk"]
rpkm_pah = rpkm.loc[:, [col for col in rpkm.columns if col in pah_samples]]

# Collect all DE genes from all RNA types
all_de_genes = set()
for df in DE.values():
    all_de_genes.update(df.index)
print(f"Total DE genes across all RNA types: {len(all_de_genes)}")

# Filter to genes present in RPKM
rpkm_de = rpkm_pah.loc[[g for g in all_de_genes if g in rpkm_pah.index], :]
# Remove genes with median zero
gene_medians = rpkm_de.median(axis=1)
rpkm_de = rpkm_de.loc[gene_medians > 0, :]

# Log2 transform
X = np.log2(rpkm_de + 1).T   # samples as rows, genes as columns
# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Select top 50 most variable genes (based on variance across samples)
gene_var = X.var(axis=0)
top_50_genes = gene_var.nlargest(50).index
X_top = X_scaled[:, [list(X.columns).index(g) for g in top_50_genes]]

# =============================================================================
# Determine optimal number of clusters
# =============================================================================
best_method, best_k, silhouette_df = find_optimal_clusters(X_top, k_min=2, k_max=10)
print(f"Optimal clustering: {best_method} with k={best_k} (silhouette={silhouette_df.loc[silhouette_df['Silhouette'].idxmax(), 'Silhouette']:.4f})")
plot_silhouette_curve(silhouette_df, os.path.join(OUTPUT_FIGURES, "silhouette_plot.pdf"))

# =============================================================================
# Perform clustering with the three methods (using best_k for each? Use same k for all for comparison)
# For consistency, we use the same k for all methods.
n_clusters = best_k
clustering_methods = {
    "K-means": KMeans(n_clusters=n_clusters, random_state=42),
    "Hierarchical": AgglomerativeClustering(n_clusters=n_clusters),
    "Spectral": SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42)
}

cluster_labels = {}
for name, model in clustering_methods.items():
    labels = model.fit_predict(X_top)
    cluster_labels[name] = labels
    print(f"{name} silhouette score: {silhouette_score(X_top, labels):.4f}")

# Save cluster assignments
for name, labels in cluster_labels.items():
    out_df = pd.DataFrame({"Sample": X.index, "Cluster": labels})
    out_df.to_csv(os.path.join(OUTPUT_RESULTS, f"{name}_clustering_results.csv"), index=False)

# =============================================================================
# PCA and t-SNE Visualization
# =============================================================================
# PCA 2D
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_top)
pca_df = pd.DataFrame(X_pca, columns=["PC1", "PC2"])
for name, labels in cluster_labels.items():
    pca_df[name] = labels
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for i, name in enumerate(cluster_labels.keys()):
    sns.scatterplot(x="PC1", y="PC2", hue=name, data=pca_df, palette="viridis", ax=axes[i])
    axes[i].set_title(f"{name} Clustering (PCA)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "PCA_clustering_results.pdf"), dpi=300)
plt.close()

# t-SNE 2D
tsne = TSNE(n_components=2, random_state=42)
X_tsne = tsne.fit_transform(X_top)
tsne_df = pd.DataFrame(X_tsne, columns=["tSNE1", "tSNE2"])
for name, labels in cluster_labels.items():
    tsne_df[name] = labels
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for i, name in enumerate(cluster_labels.keys()):
    sns.scatterplot(x="tSNE1", y="tSNE2", hue=name, data=tsne_df, palette="viridis", ax=axes[i])
    axes[i].set_title(f"{name} Clustering (t-SNE)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "tsne_clustering_results.pdf"), dpi=300)
plt.close()

# =============================================================================
# Survival Analysis
# =============================================================================
# Load clinical survival data
clinical_surv = pd.read_excel(CLINICAL_FILE, sheet_name="20241223")
clinical_surv = clinical_surv[['Seq_ID', '生存时间2021', '随访结果2021（0=存活，1=死亡，2=失访）']].copy()
clinical_surv.columns = ['Seq_ID', 'Time', 'Status']
clinical_surv = clinical_surv.dropna(subset=['Time'])
clinical_surv['Status'] = clinical_surv['Status'].astype(int)

# Survival plots for each clustering method
fig, axes = plt.subplots(3, 1, figsize=(10, 18))
for i, (name, labels) in enumerate(cluster_labels.items()):
    cluster_df = pd.DataFrame({"Sample": X.index, "Cluster": labels})
    perform_survival_analysis(cluster_df, axes[i])
    axes[i].set_title(f"{name} Clustering")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "KM_cluster_plot.pdf"), dpi=300)
plt.close()

# =============================================================================
# Clinical Correlation (Boxplots and Stacked Barplots)
# =============================================================================
# Load full clinical metadata (Sheet1)
clinical_meta = pd.read_excel(CLINICAL_FILE, sheet_name="Sheet1")
# Convert date column if needed
clinical_meta["采血日期"] = clinical_meta["采血日期"].dropna().astype(float)
clinical_meta["采血日期"] = pd.to_datetime(clinical_meta["采血日期"], unit="D", origin="1899-12-30")
clinical_meta = clinical_meta.set_index("Seq_ID")
# Select relevant columns
clinical_meta = clinical_meta.loc[:, ["分类","采血日期","性别","年龄","六mwt","心功能分级","体重rhc","身高rhc","体表面积rhc",
                                      "右心导管sbp","右心导管dbp","右心导管hr","右心房压","肺动脉平均压","肺毛细血管楔压",
                                      "心输出量","心脏指数","肺血管阻力","混合静脉血氧饱和度","ntprobnp","肌酐","尿素氮",
                                      "起始用药（2=双药联合，3=三药联合，4=单药，0=未使用靶向药物）"]]
# Rename columns to English
col_map = {
    "分类": "Category", "采血日期": "Sampling_Date", "性别": "Gender", "年龄": "Age",
    "六mwt": "SixMWT", "心功能分级": "Cardiac_Function", "体重rhc": "Weight_RHC",
    "身高rhc": "Height_RHC", "体表面积rhc": "BSA_RHC", "右心导管sbp": "RHC_SBP",
    "右心导管dbp": "RHC_DBP", "右心导管hr": "RHC_HR", "右心房压": "RAP",
    "肺动脉平均压": "mPAP", "肺毛细血管楔压": "PCWP", "心输出量": "CO",
    "心脏指数": "CI", "肺血管阻力": "PVR", "混合静脉血氧饱和度": "SvO2",
    "ntprobnp": "NT_ProBNP", "肌酐": "Creatinine", "尿素氮": "BUN",
    "起始用药（2=双药联合，3=三药联合，4=单药，0=未使用靶向药物）": "Initial_Treatment"
}
clinical_meta = clinical_meta.rename(columns=col_map)

# Merge with survival to get Time and Status (for filtering)
clinical_meta = pd.merge(clinical_meta, clinical_surv, left_index=True, right_on='Seq_ID', how='left')
clinical_meta = clinical_meta.dropna(subset=['Seq_ID']).set_index('Seq_ID')

# Load REVEAL and COMPERA scores
reveal = pd.read_csv(REVEAL_FILE, encoding='GBK')
reveal = reveal.set_index('Seq_ID')
reveal = reveal[["REVEAL_sub_groups", "COMPERA_2.0"]]
clinical_meta = pd.merge(clinical_meta, reveal, left_index=True, right_index=True, how='left')

# List of continuous variables for boxplots
continuous_vars = ["SixMWT", "Weight_RHC", "Height_RHC", "BSA_RHC", "RHC_SBP", "RHC_DBP", "RHC_HR",
                   "RAP", "mPAP", "PCWP", "CO", "CI", "PVR", "SvO2", "NT_ProBNP", "Creatinine", "BUN",
                   "Sampling_Date", "Age"]
# Categorical variables for stacked barplots
categorical_vars = ["Category", "Gender", "Cardiac_Function", "Initial_Treatment", "REVEAL_sub_groups", "COMPERA_2.0"]

# For each clustering method, generate a PDF with all clinical correlation plots
cluster_files = [
    ("K-means", "K-means_clustering_results.csv"),
    ("Hierarchical", "Hierarchical_clustering_results.csv"),
    ("Spectral", "Spectral_clustering_results.csv")
]

for method_name, csv_file in cluster_files:
    cluster_df = pd.read_csv(os.path.join(OUTPUT_RESULTS, csv_file))
    # Merge clinical data with clusters
    merged = pd.merge(clinical_meta, cluster_df, left_index=True, right_on='Sample', how='left')
    merged = merged.dropna(subset=['Cluster'])
    
    # Create PDF with all plots (boxplots + stacked barplots)
    output_pdf = os.path.join(OUTPUT_FIGURES, f"{method_name}_clinical_correlation.pdf")
    all_columns = continuous_vars + categorical_vars
    plots_per_page = 16
    total_plots = len(all_columns)
    total_pages = (total_plots + plots_per_page - 1) // plots_per_page
    
    with PdfPages(output_pdf) as pdf:
        for page in range(total_pages):
            fig, axes = plt.subplots(4, 4, figsize=(6.3, 12.6))
            axes = axes.flatten()
            start = page * plots_per_page
            end = min((page+1)*plots_per_page, total_plots)
            for idx, col in enumerate(all_columns[start:end]):
                ax = axes[idx]
                if col in continuous_vars:
                    plot_box_with_pval(merged, col, ax=ax)
                else:
                    plot_stacked_percentage(merged, col, ax=ax)
            # Hide unused subplots
            for j in range(end-start, plots_per_page):
                axes[j].axis('off')
            plt.tight_layout()
            pdf.savefig(fig, dpi=300)
            plt.close(fig)
    print(f"Saved clinical correlation plots for {method_name} to {output_pdf}")

# Also generate individual plots (as in original) for each variable separately (optional, but original code did that)
# We'll keep the PDF summary as above, but also create per-variable PDFs if needed.
# Original script had separate loops for boxplots and stacked barplots. We'll also produce them.
for method_name, csv_file in cluster_files:
    cluster_df = pd.read_csv(os.path.join(OUTPUT_RESULTS, csv_file))
    merged = pd.merge(clinical_meta, cluster_df, left_index=True, right_on='Sample', how='left')
    merged = merged.dropna(subset=['Cluster'])
    # Boxplots per continuous variable
    box_dir = os.path.join(OUTPUT_FIGURES, f"{method_name}_boxplots")
    os.makedirs(box_dir, exist_ok=True)
    for col in continuous_vars:
        fig = plot_box_with_pval(merged, col)
        if fig:
            fig.savefig(os.path.join(box_dir, f"{col}_boxplot.pdf"), dpi=300)
            plt.close(fig)
    # Stacked barplots per categorical variable
    stack_dir = os.path.join(OUTPUT_FIGURES, f"{method_name}_stacked")
    os.makedirs(stack_dir, exist_ok=True)
    for col in categorical_vars:
        fig = plot_stacked_percentage(merged, col)
        if fig:
            fig.savefig(os.path.join(stack_dir, f"{col}_stacked.pdf"), dpi=300)
            plt.close(fig)

print("All unsupervised clustering, survival, and clinical correlation analyses completed.")