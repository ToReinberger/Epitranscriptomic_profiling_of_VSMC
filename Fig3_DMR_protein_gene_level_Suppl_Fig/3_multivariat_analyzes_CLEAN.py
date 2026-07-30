import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

file = "direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER_NEW_GENENAMES.tsv"

cond = "C4"
# ---------------------------------------------------------
# filter for non-zero values
# ---------------------------------------------------------

df = pd.read_csv(file, sep="\t")
df = df[
    (df[f"gene_level_ont_{cond}"] > 0)
    & (df[f"protein_level_{cond}"] > 0)
    ]

# --------------------------------------------------------
# Construct biologically meaningful variables
# ---------------------------------------------------------

df["has_dmr"] = (df["dmr_count"] > 0).astype(int)

df["is_deg"] = (df["ont_gene_qval"] <= 0.05).astype(int)
# or for protein:
df["is_dap"] = (df["protein_qval"] <= 0.05).astype(int)

df["mean_gene_abundance"] = (
    df["gene_level_ont_C3"] + df["gene_level_ont_C4"]
) / 2

df["mean_protein_abundance"] = (
    df["protein_level_C3"] + df["protein_level_C4"]
) / 2

df["abs_gene_log2fc"] = df["ont_gene_log2FC_C4vsC3"].abs()
df["abs_protein_log2fc"] = df["protein_log2FC_C4vsC3"].abs()

variables = [
    "transcript_count",
    "median_transcript_len",
    "mean_gene_abundance",
    "mean_protein_abundance",
]

predictors = variables

model_df = df[
    ["dmr_count", "has_dmr", "is_deg", "is_dap"] + predictors
].replace([np.inf, -np.inf], np.nan).dropna()
scaler = StandardScaler()
standardized = scaler.fit_transform(model_df[predictors])
z_predictors = ["z_" + col for col in predictors]
model_df[z_predictors] = standardized


# --------------------------------------------------------
# Run binomial generalized linear model / logistic regression model
# ---------------------------------------------------------

forumla_nb = "is_dap ~ has_dmr + " + " + ".join(z_predictors)
print(forumla_nb)

deg_model = smf.glm(
    formula=forumla_nb,
    data=model_df,
    family=sm.families.Binomial()
).fit(cov_type="HC3")

"""
HC3 is a conservative robust covariance estimator. It is often preferred over HC0/HC1 in smaller 
or moderately sized datasets because it penalizes high-leverage observations more strongly.
"""


deg_results = pd.DataFrame({
    "variable": deg_model.params.index,
    "coefficient": deg_model.params.values,
    "odds_ratio": np.exp(deg_model.params.values),
    "p_value": deg_model.pvalues.values,
    "ci_low": np.exp(deg_model.conf_int()[0].values),
    "ci_high": np.exp(deg_model.conf_int()[1].values),
})


"""
Genes/proteins with DMRs have approximately 1.38-fold higher odds of being differentially abundant 
proteins compared with genes/proteins without DMRs, after adjustment for transcript count, transcript length, 
gene abundance and protein abundance.
"""

print(deg_model.summary())
print(deg_results)
deg_results.to_excel("DAPs_covariates_Logistic regression.xlsx", index=False)
# Because all predictors are z-standardized, each coefficient corresponds to a one standard deviation increase in that predictor.