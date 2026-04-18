# =============================================================================
# Multi-class Classification: PAH Subtypes (I/HPAH vs CHD-PAH vs SLE-PAH)
# =============================================================================
# This script:
#   1. Reads sample IDs, RNA name lists, RPKM expression matrix, and raw counts
#   2. Performs DESeq2 to identify differentially expressed cfRNAs (one-vs-others) using training set only
#   3. Selects stable features using LASSO with bootstrap resampling (frequency >= 75%)
#   4. Trains 14 machine learning classifiers (multi-class: IPAH, CHD, SLE)
#   5. Evaluates models using 30 iterations of stratified train-test split (70/30)
#   6. Saves feature frequency, per-class AUC/ACC results, risk scores, and LASSO coefficients
#   7. Generates confusion matrices and ROC curves for each iteration
# =============================================================================

import os
import numpy as np
import pandas as pd
from collections import Counter
from datetime import datetime
from scipy.stats import sem, t
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.linear_model import LassoCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_curve, auc, accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import label_binarize

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# =============================================================================
# User Configuration (relative to project root)
# =============================================================================
DATA_DIR = "data"
OUTPUT_RESULTS_DIR = "results/ML/3class"
OUTPUT_FIGURES_DIR = "figures/ML/3class"

# Input files
SAMPLE_ID_FILES = {
    "CHD": os.path.join(DATA_DIR, "CHD_not4pnk_ID.txt"),
    "SLE": os.path.join(DATA_DIR, "SLE_not4pnk_ID.txt"),
    "NOR": os.path.join(DATA_DIR, "NOR_not4pnk_ID.txt"),
    "IPAH": os.path.join(DATA_DIR, "IPAH_not4pnk_ID.txt"),
    "IPAH_t4pnk": os.path.join(DATA_DIR, "IPAH_t4pnk_ID.txt"),
    "NOR_t4pnk": os.path.join(DATA_DIR, "NOR_t4pnk_ID.txt")
}
RNA_NAME_FILES = {
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
RPKM_FILE = os.path.join(DATA_DIR, "df_rpkm_300W.csv")
COUNTS_FILE = os.path.join(DATA_DIR, "df_count_20241117.csv")

# Output subdirectories
FEATURE_AUC_ACC_DIR = os.path.join(OUTPUT_RESULTS_DIR, "Tables/Feature_AUC_ACC")
RISK_SCORE_DIR = os.path.join(OUTPUT_RESULTS_DIR, "Tables/Risk_Score")
CONFUSION_MATRIX_DIR = os.path.join(OUTPUT_FIGURES_DIR, "confusion_matrix")
AUC_PLOT_DIR = os.path.join(OUTPUT_FIGURES_DIR, "AUC_plot")

# Create output directories
for d in [FEATURE_AUC_ACC_DIR, RISK_SCORE_DIR, CONFUSION_MATRIX_DIR, AUC_PLOT_DIR]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# Helper Functions
# =============================================================================
def read_file(file_path):
    """Read a text file and return lines as a list (strip newlines)."""
    with open(file_path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f.readlines()]

def ci95(data):
    """Calculate 95% confidence interval (mean, lower, upper)."""
    m = np.mean(data)
    s = sem(data)
    ci = t.interval(0.95, len(data)-1, loc=m, scale=s)
    return m, ci[0], ci[1]

def safe_ci95(x):
    """Safe 95% CI calculation, returns (nan, nan, nan) if all NaN."""
    arr = np.asarray(x, dtype=float)
    if arr.size == 0 or np.all(np.isnan(arr)):
        return (np.nan, np.nan, np.nan)
    vals = arr[~np.isnan(arr)]
    n = vals.size
    mean = np.mean(vals)
    if n < 2:
        return (mean, np.nan, np.nan)
    sem_val = np.std(vals, ddof=1) / np.sqrt(n)
    z = 1.96
    return (mean, mean - z * sem_val, mean + z * sem_val)

def _make_empty_results_analyze(iteration):
    """Return empty placeholder results for failed analyses."""
    all_models = list(models.keys())
    results_placeholder = {m: {f'Class_{i+1}': np.nan for i in range(3)} for m in all_models}
    acc_placeholder = {m: np.nan for m in all_models}
    train_auc_placeholder = {m: {f'Class_{i+1}': np.nan for i in range(3)} for m in all_models}
    train_acc_placeholder = {m: np.nan for m in all_models}
    selected_feats_placeholder = []
    risk_df_placeholder = pd.DataFrame(columns=["Sample", "Model", "RNA_Type", "Risk_Score", "True_Label", "Predicted_Label", "Iteration"])
    lasso_coefs_placeholder = pd.DataFrame()
    return (results_placeholder, acc_placeholder, train_auc_placeholder, train_acc_placeholder,
            selected_feats_placeholder, risk_df_placeholder, lasso_coefs_placeholder)

def process_group13(target_group, reference_groups, sampleID, X_train):
    """
    Perform DESeq2 analysis for one-vs-others (multi-class).
    """
    idx = sampleID[target_group] + [s for g in reference_groups for s in sampleID[g]]
    df_counts = X_train[X_train.index.isin(idx)]
    df_meta = pd.DataFrame(index=df_counts.index)
    df_meta.loc[df_meta.index.isin(sampleID[target_group]), 'Sample_Group'] = target_group
    for g in reference_groups:
        df_meta.loc[df_meta.index.isin(sampleID[g]), 'Sample_Group'] = 'others'
    # Filter low-expressed genes
    genes_to_keep = df_counts.columns[df_counts.median(axis=0) >= 5]
    df_counts = df_counts[genes_to_keep]
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=df_counts, metadata=df_meta, design="~Sample_Group",
                       refit_cooks=True, inference=inference)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["Sample_Group", target_group, "others"], inference=inference)
    ds.summary()
    return ds.results_df

# =============================================================================
# Define 14 Machine Learning Models (multi-class compatible)
# =============================================================================
models = {
    "GLMNETRIDGE": LogisticRegression(penalty='l2', solver='liblinear', C=1.0,
                                      class_weight="balanced", random_state=42, max_iter=500),
    "GLMNETLASSO": LogisticRegression(penalty='l1', solver='liblinear', C=1.0,
                                      class_weight="balanced", random_state=42, max_iter=500),
    "SVMLIN": SVC(kernel="linear", probability=True, class_weight="balanced", random_state=42),
    "SVMRAD": SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=42),
    "RF": RandomForestClassifier(n_estimators=100, max_depth=5, class_weight="balanced", random_state=42),
    "EXTRATREES": ExtraTreesClassifier(n_estimators=100, max_depth=5, class_weight="balanced", random_state=42),
    "NNET": MLPClassifier(hidden_layer_sizes=(100,), max_iter=500, random_state=42),
    "LDA": LinearDiscriminantAnalysis(),
    "C5": DecisionTreeClassifier(max_depth=5, class_weight="balanced", random_state=42),
    "KNN": KNeighborsClassifier(n_neighbors=5),
    "NB": GaussianNB(),
    "RPART": DecisionTreeClassifier(class_weight="balanced", random_state=42),
    "GLM": LogisticRegression(penalty=None, solver='lbfgs', class_weight="balanced", random_state=42, max_iter=500),
    "GBM": GradientBoostingClassifier(n_estimators=100, random_state=42)
}

# =============================================================================
# Main Analysis Function (Multi-class)
# =============================================================================
def analyze_single_DE(DE_data, RNA_type, compare_type, rpkm, sampleID, X_train, X_test,
                      y_train, y_test, n_classes=3, stability_threshold=0.75, n_bootstraps=100,
                      return_metrics=False, return_features=False, iteration=None):
    """
    Feature selection (LASSO stability) and model training for multi-class classification.
    """
    # Collect DE genes from each one-vs-others comparison (top 100 per comparison)
    all_genes = set()
    for key, res in DE_data.items():
        if isinstance(res, pd.DataFrame) and 'padj' in res.columns and 'log2FoldChange' in res.columns:
            sig = res.loc[res['padj'] < 0.05]
            n_take = min(100, sig.shape[0])
            if n_take > 0:
                idxs = sig['log2FoldChange'].abs().nlargest(n_take).index
                sig_top = sig.loc[idxs, :]
                all_genes.update(sig_top.index)
    filtered_genes = list(all_genes)
    if len(filtered_genes) == 0:
        print("Warning: No DE genes found. Skipping.")
        return _make_empty_results_analyze(iteration)

    filtered_genes_in_rpkm = [g for g in filtered_genes if g in rpkm.index]
    if len(filtered_genes_in_rpkm) == 0:
        print("Warning: No DE genes found in rpkm matrix. Skipping.")
        return _make_empty_results_analyze(iteration)

    # Prepare training expression matrix
    X_train = rpkm.loc[filtered_genes_in_rpkm, X_train.index].T
    y_train_df = pd.DataFrame(index=X_train.index)
    group_mapping = {"IPAH_not4pnk": 1, "CHD_not4pnk": 2, "SLE_not4pnk": 3}
    for group, label in group_mapping.items():
        y_train_df.loc[y_train_df.index.isin(sampleID[group]), 'label'] = label
    y_train = y_train_df['label'].astype(int)

    # Log2 transform
    X_train_scaled = np.log2(X_train + 1)
    X_train_scaled_df = pd.DataFrame(X_train_scaled, columns=X_train.columns, index=X_train.index)

    # LASSO stability selection
    print(f"Performing {n_bootstraps} LASSO bootstraps for stability selection...")
    feature_counts = pd.Series(0, index=X_train_scaled_df.columns)
    sss = StratifiedShuffleSplit(n_splits=n_bootstraps, test_size=0.2, random_state=42)
    any_successful = False
    for i, (train_idx, val_idx) in enumerate(sss.split(X_train_scaled_df, y_train)):
        X_boot = X_train_scaled_df.iloc[train_idx]
        y_boot = y_train.iloc[train_idx]
        try:
            lasso = LassoCV(cv=5, random_state=i).fit(X_boot, y_boot)
        except Exception:
            continue
        any_successful = True
        lasso_coef = lasso.coef_
        selected = X_boot.columns[lasso_coef != 0]
        feature_counts[selected] += 1

    if not any_successful:
        print("Warning: All LASSO bootstraps failed. Skipping.")
        return _make_empty_results_analyze(iteration)

    selection_freq = feature_counts / n_bootstraps
    stable_features = selection_freq[selection_freq >= stability_threshold].index.tolist()
    print(f"Stable features (freq >= {stability_threshold}): {len(stable_features)}")
    if len(stable_features) == 0:
        print("Warning: No stable features selected. Skipping.")
        return _make_empty_results_analyze(iteration)

    X_train_selected = X_train_scaled_df[stable_features]

    # Prepare test data
    X_test = rpkm.loc[filtered_genes_in_rpkm, X_test.index].T
    y_test_df = pd.DataFrame(index=X_test.index)
    for group, label in group_mapping.items():
        y_test_df.loc[y_test_df.index.isin(sampleID[group]), 'label'] = label
    y_test = y_test_df['label'].astype(int)

    X_test_scaled = np.log2(X_test + 1)
    X_test_scaled_df = pd.DataFrame(X_test_scaled, columns=X_test.columns, index=X_test.index)
    X_test_selected = X_test_scaled_df[stable_features]

    # Train and evaluate models
    acc_dict = {}
    auc_dict = {model: {} for model in models}
    train_auc_dict = {model: {} for model in models}
    train_acc_dict = {}
    all_risk_scores = []
    lasso_all_coefs = []

    current_date = datetime.now().strftime('%Y-%m-%d')
    if iteration is not None:
        pdf_file_name = os.path.join(CONFUSION_MATRIX_DIR, f"Confusion_Matrix_AllModels_{RNA_type}_{compare_type}_iteration_{iteration}_{current_date}.pdf")
    else:
        pdf_file_name = os.path.join(CONFUSION_MATRIX_DIR, f"Confusion_Matrix_AllModels_{RNA_type}_{compare_type}_{current_date}.pdf")

    n_models = len(models)
    ncols = 3
    nrows = int(np.ceil(max(1, n_models) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 6 * nrows))
    axes = axes.flatten()

    with PdfPages(pdf_file_name) as pdf:
        for idx, (name, model) in enumerate(models.items()):
            # Train model
            try:
                model.fit(X_train_selected, y_train)
            except Exception:
                print(f"Warning: Model {name} failed to fit. Skipping.")
                acc_dict[name] = np.nan
                train_acc_dict[name] = np.nan
                for cls in range(1, 4):
                    auc_dict[name][f'Class_{cls}'] = np.nan
                    train_auc_dict[name][f'Class_{cls}'] = np.nan
                continue

            # Training set performance
            try:
                y_train_prob = model.predict_proba(X_train_selected)
                y_train_bin = label_binarize(y_train, classes=[1, 2, 3])
                for cls in range(3):
                    fpr, tpr, _ = roc_curve(y_train_bin[:, cls], y_train_prob[:, cls])
                    train_auc_dict[name][f'Class_{cls+1}'] = auc(fpr, tpr)
                y_train_pred = model.predict(X_train_selected)
                train_acc_dict[name] = accuracy_score(y_train, y_train_pred)
            except Exception:
                train_acc_dict[name] = np.nan
                for cls in range(1, 4):
                    train_auc_dict[name][f'Class_{cls}'] = np.nan

            # Test set performance (multi-class: argmax)
            try:
                y_prob = model.predict_proba(X_test_selected)
                y_test_bin = label_binarize(y_test, classes=[1, 2, 3])
                for cls in range(3):
                    fpr, tpr, _ = roc_curve(y_test_bin[:, cls], y_prob[:, cls])
                    auc_dict[name][f'Class_{cls+1}'] = auc(fpr, tpr)
                y_pred = np.argmax(y_prob, axis=1) + 1
                acc_dict[name] = accuracy_score(y_test, y_pred)
            except Exception:
                acc_dict[name] = np.nan
                for cls in range(1, 4):
                    auc_dict[name][f'Class_{cls}'] = np.nan
                y_pred = np.full(len(y_test), np.nan)

            # Risk score (simplified: use maximum predicted probability as risk)
            try:
                risk_score = np.max(y_prob, axis=1)
                risk_df = pd.DataFrame({
                    "Sample": X_test_selected.index,
                    "Model": name,
                    "RNA_Type": RNA_type,
                    "Risk_Score": risk_score,
                    "True_Label": y_test.values,
                    "Predicted_Label": y_pred,
                    "Iteration": iteration
                })
                all_risk_scores.append(risk_df)
            except Exception:
                pass

            # Draw confusion matrix
            try:
                cm = confusion_matrix(y_test, y_pred)
                ax = axes[idx]
                im = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
                ax.figure.colorbar(im, ax=ax)
                ax.set(xticks=np.arange(3), yticks=np.arange(3),
                       xticklabels=["IPAH", "CHD", "SLE"], yticklabels=["IPAH", "CHD", "SLE"],
                       title=f'Confusion Matrix for {name}', ylabel='True label', xlabel='Predicted label')
                thresh = cm.max() / 2. if cm.size > 0 else 0
                for i0 in range(cm.shape[0]):
                    for j0 in range(cm.shape[1]):
                        ax.text(j0, i0, format(cm[i0, j0], 'd'),
                                ha="center", va="center",
                                color="white" if cm[i0, j0] > thresh else "black")
            except Exception:
                pass

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    # Prepare return values
    try:
        risk_concat = pd.concat(all_risk_scores, axis=0, ignore_index=True) if all_risk_scores else pd.DataFrame()
    except Exception:
        risk_concat = pd.DataFrame()

    if return_metrics and return_features:
        return auc_dict, acc_dict, train_auc_dict, train_acc_dict, stable_features, risk_concat, lasso_all_coefs
    elif return_metrics:
        return auc_dict, acc_dict
    elif return_features:
        return stable_features
    else:
        return None

# =============================================================================
# Main Loop: Repeated Classification (r=30 iterations)
# =============================================================================
def run_repeated_classification(df_merged, RNA_type, rpkm, sampleID, r=30):
    """
    Run multi-class classification r times with stratified train-test splits.
    """
    all_models = list(models.keys())
    auc_DE13 = {model: {f'Class_{i+1}': [] for i in range(3)} for model in all_models}
    acc_DE13 = {model: [] for model in all_models}
    train_auc_DE13 = {model: {f'Class_{i+1}': [] for i in range(3)} for model in all_models}
    train_acc_DE13 = {model: [] for model in all_models}
    selected_feature_list = []
    all_risk_score_dfs = []
    all_lasso_coef_dfs = []

    for i in range(r):
        print(f"===== Iteration {i+1}/{r} =====")
        # Prepare counts matrix for this RNA type
        df_counts = df_merged[RNA_type].T.astype(int)
        meta = pd.DataFrame(index=df_counts.index)
        target_groups = ["IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk"]
        for group in target_groups:
            meta.loc[meta.index.isin(sampleID[group]), 'Sample_Group'] = group
        df_counts = df_counts[meta['Sample_Group'].notna()]
        df_counts['Sample_Group'] = meta['Sample_Group']

        # Train-test split (stratified)
        train_df, test_df = train_test_split(df_counts, test_size=0.3, stratify=df_counts['Sample_Group'], random_state=i)
        X_train = train_df.drop(columns=['Sample_Group'])
        y_train = train_df['Sample_Group']
        X_test = test_df.drop(columns=['Sample_Group'])
        y_test = test_df['Sample_Group']

        # DESeq2 one-vs-others for each subtype (using training set only)
        DE13 = {}
        all_groups = ["IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk"]
        for group in all_groups:
            rest_groups = [g for g in all_groups if g != group]
            DE13[group.split("_not4pnk")[0]] = process_group13(group, rest_groups, sampleID, X_train)

        compare_type = f"1v3_multi_class_train_test_iter{i+1}"
        auc_dict, acc_dict, train_auc_dict, train_acc_dict, selected, risk_df, coef_df = analyze_single_DE(
            DE13, RNA_type, compare_type, rpkm, sampleID,
            X_train, X_test, y_train, y_test,
            return_metrics=True, return_features=True, iteration=i+1
        )
        all_risk_score_dfs.append(risk_df)
        # all_lasso_coef_dfs.append(coef_df)  # coef_df not collected here; can be added

        for model in all_models:
            for cls in auc_dict[model]:
                auc_DE13[model][cls].append(auc_dict[model][cls])
                train_auc_DE13[model][cls].append(train_auc_dict[model][cls])
            acc_DE13[model].append(acc_dict[model])
            train_acc_DE13[model].append(train_acc_dict[model])
        selected_feature_list.extend(selected)

    # Print mean performance with 95% CI
    print("\n=== DE13 Test Set Performance (mean ± 95% CI) ===")
    for model in all_models:
        acc_m, acc_l, acc_u = safe_ci95(acc_DE13[model])
        print(f"{model} ACC = {acc_m:.3f} [{acc_l:.3f}, {acc_u:.3f}]")
        for cls in auc_DE13[model]:
            auc_m, auc_l, auc_u = safe_ci95(auc_DE13[model][cls])
            print(f"  {cls}: AUC = {auc_m:.3f} [{auc_l:.3f}, {auc_u:.3f}]")

    print("\n=== DE13 Training Set Performance (mean ± 95% CI) ===")
    for model in all_models:
        acc_m, acc_l, acc_u = safe_ci95(train_acc_DE13[model])
        print(f"{model} ACC = {acc_m:.3f} [{acc_l:.3f}, {acc_u:.3f}]")
        for cls in train_auc_DE13[model]:
            auc_m, auc_l, auc_u = safe_ci95(train_auc_DE13[model][cls])
            print(f"  {cls}: AUC = {auc_m:.3f} [{auc_l:.3f}, {auc_u:.3f}]")

    # Save results
    current_date = datetime.now().strftime('%Y-%m-%d')
    # Feature frequency
    feature_counter = Counter(selected_feature_list)
    feature_df = pd.DataFrame(feature_counter.items(), columns=["Feature", "Frequency"]).sort_values(by="Frequency", ascending=False)
    feature_df.to_csv(os.path.join(FEATURE_AUC_ACC_DIR, f"features_frequency_DE13_{RNA_type}_{current_date}.csv"), index=False)

    # Save AUC (test and train) and ACC
    def save_auc_csv(auc_dict, set_type):
        rows = []
        for it in range(r):
            row = {'Iteration': it + 1}
            for model in all_models:
                for cls in auc_dict[model]:
                    row[f'{model}_{cls}'] = auc_dict[model][cls][it]
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(FEATURE_AUC_ACC_DIR, f"AUC_{set_type}_DE13_{RNA_type}_{current_date}.csv"), index=False)

    def save_acc_csv(acc_dict, set_type):
        df = pd.DataFrame(acc_dict)
        df['Iteration'] = range(1, r + 1)
        df = df[['Iteration'] + all_models]
        df.to_csv(os.path.join(FEATURE_AUC_ACC_DIR, f"ACC_{set_type}_DE13_{RNA_type}_{current_date}.csv"), index=False)

    save_auc_csv(auc_DE13, "test")
    save_auc_csv(train_auc_DE13, "train")
    save_acc_csv(acc_DE13, "test")
    save_acc_csv(train_acc_DE13, "train")

    # Risk scores
    all_risk_df = pd.concat(all_risk_score_dfs, axis=0, ignore_index=True) if all_risk_score_dfs else pd.DataFrame()
    all_risk_df.to_csv(os.path.join(RISK_SCORE_DIR, f"RiskScore_DE13_All_{RNA_type}_{current_date}.csv"), index=False)

    return auc_DE13, train_auc_DE13, acc_DE13, train_acc_DE13, feature_df

# =============================================================================
# Data Loading and Preparation
# =============================================================================
# Read sample IDs
sampleID = {}
for group, path in SAMPLE_ID_FILES.items():
    sampleID[group] = read_file(path)

# Read RNA name lists
RNAname = {}
for rna, path in RNA_NAME_FILES.items():
    RNAname[rna] = read_file(path)

# Read RPKM matrix
rpkm = pd.read_csv(RPKM_FILE, index_col=0)
# Keep only samples from all groups of interest (including NOR for completeness, but only three disease groups will be used later)
all_sample_groups = ['CHD_not4pnk', 'SLE_not4pnk', 'NOR_not4pnk', 'IPAH_not4pnk', 'IPAH_t4pnk', 'NOR_t4pnk']
columns_to_select = []
for group in all_sample_groups:
    columns_to_select.extend(sampleID.get(group, []))
rpkm = rpkm[columns_to_select]

# Read counts matrix
counts = pd.read_csv(COUNTS_FILE, index_col=0)
counts = counts.loc[:, rpkm.columns]

# Build merged count data per RNA type (for all samples, but will subset later)
rna_types = ['rsRNA', 'ysRNA', 'tRFs', 'miRNA', 'mttRNA', 'snRNA', 'snoRNA', 'mRNA', 'lncRNA']
groups_for_counts = ['CHD_not4pnk', 'SLE_not4pnk', 'NOR_not4pnk', 'IPAH_not4pnk', 'IPAH_t4pnk', 'NOR_t4pnk']
counts_data = {}
for group in groups_for_counts:
    for rna_type in rna_types:
        df_name = f'counts_{group}_{rna_type}'
        RNAname_idx = RNAname[rna_type]
        sampleID_idx = sampleID[group]
        counts_data[df_name] = counts.loc[RNAname_idx, sampleID_idx]

# Merge counts for each RNA type across groups (only the three disease groups will be used in classification)
# But we keep all for consistency; run_repeated_classification will filter to IPAH/CHD/SLE only.
df_merged = {}
for rna_type in rna_types:
    dfs = [counts_data[f'counts_{group}_{rna_type}'] for group in ['CHD_not4pnk', 'SLE_not4pnk', 'NOR_not4pnk', 'IPAH_not4pnk']]
    merged_df = pd.concat(dfs, axis=1)
    df_merged[rna_type] = merged_df

# =============================================================================
# Run for each RNA type
# =============================================================================
run_repeated_classification(df_merged, RNA_type="miRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="rsRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="ysRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="snRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="snoRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="tRFs", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="mRNA", rpkm=rpkm, sampleID=sampleID, r=30)
run_repeated_classification(df_merged, RNA_type="lncRNA", rpkm=rpkm, sampleID=sampleID, r=30)