# =============================================================================
# Multi-modal Integration: Stacking and Weighted Fusion for PAH vs NOR
# =============================================================================
# This script:
#   1. Loads pre-computed risk scores from three modalities (cfRNA, mutation, microbe)
#   2. Selects the best base model per modality based on test AUC (from individual classification results)
#   3. Aggregates risk scores across iterations (median) and aligns samples
#   4. Performs intermediate stacking (logistic regression with cross-validation) to compute OOF predictions
#   5. Compares fusion performance with single modalities
#   6. Generates ROC curves, coefficient bar plots, and permutation importance
#   7. Explores weighted fusion by grid search over weight combinations
#   8. Visualizes weight-AUC heatmap
# =============================================================================

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import seaborn as sns
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.inspection import permutation_importance
from collections import defaultdict

# =============================================================================
# User Configuration (relative to project root)
# =============================================================================
DATA_DIR = "data"
RESULTS_DIR = "results/ML"
FIGURES_DIR = "figures/ML"

# Paths to risk score and AUC result files (relative to project root)
paths = {
    "microbe": {
        "risk": os.path.join(RESULTS_DIR, "microbe/Tables/Risk_Score/RiskScore_All_*.csv"),
        "auc":  os.path.join(RESULTS_DIR, "microbe/Tables/Feature_AUC_ACC/AUC_test_results_*.csv"),
    },
    "mutation": {
        "risk": os.path.join(RESULTS_DIR, "mutation/Tables/Risk_Score/RiskScore_All_*.csv"),
        "auc":  os.path.join(RESULTS_DIR, "mutation/Tables/Feature_AUC_ACC/AUC_test_results_*.csv"),
    },
    "cfRNA": {
        "risk": os.path.join(RESULTS_DIR, "Tables/Risk_Score/NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100/RiskScore_All_all_lasso_selected_feature_*.csv"),
        "auc":  os.path.join(RESULTS_DIR, "Tables/Feature_AUC_ACC/NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100/AUC_test_results_all_lasso_selected_feature_*.csv"),
    }
}

# Output directories for fusion results
OUTPUT_TABLES = os.path.join(RESULTS_DIR, "Tables/Combine_cfRNA_mutation_microbe")
OUTPUT_FIGURES = os.path.join(FIGURES_DIR, "Combine_cfRNA_mutation_microbe")

# Create output directories
for d in [OUTPUT_TABLES, OUTPUT_FIGURES]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# Helper Functions
# =============================================================================
def latest_file(pattern):
    """Return the most recently modified file matching the pattern."""
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No file matches pattern: {pattern}")
    return max(files, key=os.path.getmtime)

# =============================================================================
# Step 1: Identify best base model per modality (based on test AUC)
# =============================================================================
best_model_by_modality = {}
for modality, pp in paths.items():
    auc_file = latest_file(pp["auc"])
    auc_df = pd.read_csv(auc_file)
    # Columns are model names; Iteration column present
    model_cols = [c for c in auc_df.columns if c != "Iteration"]
    means = auc_df[model_cols].mean(axis=0).sort_values(ascending=False)
    best_model_by_modality[modality] = means.index[0]
    print(f"Best model for {modality}: {best_model_by_modality[modality]} (AUC={means.iloc[0]:.3f})")

# =============================================================================
# Step 2: Load risk scores for the best model of each modality
# =============================================================================
risk_tables = {}
label_by_sample = None

for modality, pp in paths.items():
    risk_file = latest_file(pp["risk"])
    df = pd.read_csv(risk_file)
    # Keep only the best model for this modality
    df = df[df["Model"] == best_model_by_modality[modality]].copy()

    # Map True_Label: 1=NOR, 2=PAH -> convert to 0/1 (0=NOR, 1=PAH)
    lbl = df.groupby("Sample")["True_Label"].median().map({1: 0, 2: 1})
    if label_by_sample is None:
        label_by_sample = lbl
    else:
        # Combine labels (should be consistent)
        label_by_sample = label_by_sample.combine_first(lbl)

    # Aggregate risk scores across iterations (use median for robustness)
    risk = df.groupby("Sample")["Risk_Score"].median().rename(f"risk_{modality}")
    risk_tables[modality] = risk

# Combine all modalities (inner join: keep samples present in all three)
X = pd.concat(risk_tables.values(), axis=1, join="inner")
y = label_by_sample.reindex(X.index).astype(int)

print(f"Fusion input shape: {X.shape}")
print(f"Positive class rate: {y.mean():.3f}")

# =============================================================================
# Step 3: Intermediate stacking (logistic regression with GroupKFold)
# =============================================================================
gkf = GroupKFold(n_splits=5)
groups = X.index  # each sample is its own group (to avoid leakage)

oof_pred = pd.Series(index=X.index, dtype=float)

for fold, (tr, te) in enumerate(gkf.split(X, y, groups), 1):
    X_tr, y_tr = X.iloc[tr], y.iloc[tr]
    X_te, y_te = X.iloc[te], y.iloc[te]

    base_lr = LogisticRegression(solver="liblinear", penalty="l2", class_weight="balanced", max_iter=1000)
    param_grid = {"C": [0.01, 0.03, 0.1, 0.3, 1, 3, 10]}
    clf = GridSearchCV(base_lr, param_grid, scoring="roc_auc", cv=5, n_jobs=-1)
    clf.fit(X_tr, y_tr)

    proba = clf.predict_proba(X_te)[:, 1]
    oof_pred.iloc[te] = proba
    fold_auc = roc_auc_score(y_te, proba)
    print(f"Fold {fold} AUC = {fold_auc:.3f}, best C = {clf.best_params_['C']}")

overall_auc = roc_auc_score(y, oof_pred)
print(f"\n=== Fusion (stacking) OOF AUC = {overall_auc:.3f} ===")

# =============================================================================
# Step 4: Single-modality baseline AUC (using the same aggregated risk scores)
# =============================================================================
baseline_auc = {}
for col in X.columns:
    baseline_auc[col.replace("risk_", "")] = roc_auc_score(y, X[col])
print("Single-modality AUCs:", baseline_auc)

# =============================================================================
# Step 5: Plot ROC curves (fusion vs single modalities)
# =============================================================================
plt.figure(figsize=(6, 6))
# Fusion curve
fpr, tpr, _ = roc_curve(y, oof_pred)
plt.plot(fpr, tpr, lw=2, label=f"Fusion (AUC={auc(fpr, tpr):.3f})")

# Single-modality curves
for col in X.columns:
    fpr, tpr, _ = roc_curve(y, X[col])
    label = col.replace("risk_", "")
    plt.plot(fpr, tpr, lw=1, linestyle="--", label=f"{label} (AUC={auc(fpr, tpr):.3f})")

plt.plot([0, 1], [0, 1], lw=1, color="gray", linestyle="-")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("PAH vs NOR – Intermediate Integration")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "fusion_ROC.pdf"), dpi=300)
print(f"Saved: {os.path.join(OUTPUT_FIGURES, 'fusion_ROC.pdf')}")
plt.close()

# =============================================================================
# Step 6: Save fusion OOF predictions for downstream use
# =============================================================================
out_csv = os.path.join(OUTPUT_TABLES, "fusion_oof_predictions.csv")
pd.DataFrame({
    "Sample": X.index,
    "TrueLabel": y.values,
    "FusionScore": oof_pred.values,
    "risk_microbe": X["risk_microbe"].values,
    "risk_mutation": X["risk_mutation"].values,
    "risk_cfRNA": X["risk_cfRNA"].values,
}).to_csv(out_csv, index=False)
print(f"Saved OOF predictions to {out_csv}")

# =============================================================================
# Step 7: Weighted fusion grid search (exhaustive weights)
# =============================================================================
# Load the saved predictions (or use current data)
df = pd.read_csv(out_csv)
y = df["TrueLabel"].values
risk_cfRNA = df["risk_cfRNA"].values
risk_mutation = df["risk_mutation"].values
risk_microbe = df["risk_microbe"].values

weights = np.arange(0, 1.01, 0.1)
results = []
for w_cfRNA in weights:
    for w_microbe in weights:
        w_mutation = 1.0 - w_cfRNA - w_microbe
        if 0 <= w_mutation <= 1.0 and abs(w_cfRNA + w_mutation + w_microbe - 1.0) < 1e-6:
            fusion_score = w_cfRNA * risk_cfRNA + w_mutation * risk_mutation + w_microbe * risk_microbe
            auc_val = roc_auc_score(y, fusion_score)
            results.append({
                "w_cfRNA": round(w_cfRNA, 1),
                "w_mutation": round(w_mutation, 1),
                "w_microbe": round(w_microbe, 1),
                "AUC": round(auc_val, 4)
            })

result_df = pd.DataFrame(results)
result_df = result_df.sort_values(by="AUC", ascending=False).reset_index(drop=True)
print("Best weight combination (highest AUC):")
print(result_df.head(1))

# Save weight combination results
result_df.to_csv(os.path.join(OUTPUT_TABLES, "weight_combination_auc_results.csv"), index=False)

# =============================================================================
# Step 8: Heatmap of weighted fusion AUC (w_cfRNA vs w_microbe)
# =============================================================================
pivot_df = result_df.pivot(index="w_cfRNA", columns="w_microbe", values="AUC")

plt.figure(figsize=(10, 8))
ax = sns.heatmap(pivot_df, annot=True, cmap="YlGnBu", fmt=".3f",
                 cbar_kws={"label": "AUC"},
                 annot_kws={"fontsize": 10})
plt.title("Weighted Fusion AUC (w_mutation = 1 - w_cfRNA - w_microbe)", fontsize=12)
plt.xlabel("w_microbe (0~1, step 0.1)")
plt.ylabel("w_cfRNA (0~1, step 0.1)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "weight_auc_heatmap.pdf"), dpi=300)
print(f"Saved: {os.path.join(OUTPUT_FIGURES, 'weight_auc_heatmap.pdf')}")
plt.close()

# =============================================================================
# Step 9: Logistic regression coefficients and permutation importance
# =============================================================================
# Use the best weight from grid search (or use stacking predictions? Here we use weighted fusion)
best_weights = result_df.iloc[0]
fusion_score_best = (best_weights["w_cfRNA"] * risk_cfRNA +
                     best_weights["w_mutation"] * risk_mutation +
                     best_weights["w_microbe"] * risk_microbe)

# Logistic regression on the three risk scores (to compare with weighted fusion)
X_log = np.vstack([risk_cfRNA, risk_mutation, risk_microbe]).T
log_reg = LogisticRegression(solver="liblinear", class_weight="balanced", max_iter=1000)
log_reg.fit(X_log, y)
coefs = log_reg.coef_[0]
features = ["cfRNA", "Mutation", "Microbe"]

# Bar plot of coefficients
plt.figure(figsize=(5, 4))
plt.bar(features, coefs, color=["#1f77b4", "#ff7f0e", "#2ca02c"])
plt.ylabel("Coefficient")
plt.title("Logistic Regression Coefficients (Direct Integration)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "fusion_logreg_coeffs.pdf"), dpi=300)
print(f"Saved: {os.path.join(OUTPUT_FIGURES, 'fusion_logreg_coeffs.pdf')}")
plt.close()

# Permutation importance
result_perm = permutation_importance(log_reg, X_log, y, n_repeats=30, random_state=42)
importances = result_perm.importances_mean
stds = result_perm.importances_std

plt.figure(figsize=(5, 4))
plt.bar(features, importances, yerr=stds, color=["#1f77b4", "#ff7f0e", "#2ca02c"],
        alpha=0.8, edgecolor="black", linewidth=0.5)
plt.ylabel("Permutation Importance (mean decrease in score)")
plt.title("Feature Importance (Permutation)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "fusion_permutation_importance.pdf"), dpi=300)
print(f"Saved: {os.path.join(OUTPUT_FIGURES, 'fusion_permutation_importance.pdf')}")
plt.close()

# =============================================================================
# Optional: Plot weighted fusion ROC curve (best weights)
# =============================================================================
plt.figure(figsize=(6, 6))
fpr, tpr, _ = roc_curve(y, fusion_score_best)
plt.plot(fpr, tpr, lw=2, label=f"Weighted Fusion (AUC={auc(fpr, tpr):.3f})")
for name, score in [("cfRNA", risk_cfRNA), ("Mutation", risk_mutation), ("Microbe", risk_microbe)]:
    fpr, tpr, _ = roc_curve(y, score)
    plt.plot(fpr, tpr, lw=1.5, linestyle="--", label=f"{name} (AUC={roc_auc_score(y, score):.3f})")
plt.plot([0, 1], [0, 1], lw=1, color="gray", linestyle="-")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("PAH vs NOR – Weighted Fusion (Optimal Weights)")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_FIGURES, "fusion_weighted_ROC.pdf"), dpi=300)
print(f"Saved: {os.path.join(OUTPUT_FIGURES, 'fusion_weighted_ROC.pdf')}")
plt.close()

print("All fusion analyses completed.")