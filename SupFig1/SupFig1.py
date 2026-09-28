# -*- coding: utf-8 -*-
"""
FigSup1: The initial selection coefficient of fitness functions correlates with the length of the gene of interest.
"""
#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy import stats
import seaborn as sns
#%% Importing data

def load_dfs_to_dict(folder_path):
    dfs_dict = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            key = os.path.splitext(filename)[0]  # Remove .csv extension
            filepath = os.path.join(folder_path, filename)
            dfs_dict[key] = pd.read_csv(filepath)
    return dfs_dict

processed_data = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Processed_data/Processed_data_SC/"
updated_dfs = load_dfs_to_dict(processed_data)

# Add the expression values from deboer
DeBoer_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/DeBoer_Expression.csv"
Expression_file = pd.read_csv(DeBoer_file, sep=';')

expression_df = Expression_file[['Seq', 'Expression']]

# Merge 'Expression' into each DataFrame in the dictionary
for key in updated_dfs:
    updated_dfs[key] = updated_dfs[key].merge(expression_df, on='Seq', how='left')
    
# Remove some extreme values
updated_dfs["CDC25"] = updated_dfs["CDC25"].drop([407, 402, 463]) #Really low s_T0_T7, really different from the two other replicates
updated_dfs["PDC2"] = updated_dfs["PDC2"].drop([383, 191, 0]) # Promoter #1 is very low and breaks the relationship. Not biologically possible points, the expression is probably not the same.

#%% Figsup1
def compute_y_at_one(dfs_dict, frac=0.4):
    results = []

    for name, df in dfs_dict.items():
        df = df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["Expression", "s_T0_T7", "Replicate"]
        )
        df_comp = df[df["Type"] == "Comp"]

        x_data = df_comp["Expression"].values
        y_data = df_comp["s_T0_T7"].values

        if len(x_data) == 0:
            results.append({"dataset": name, "initial_s": np.nan})
            continue

        # LOESS fit
        loess_fit = lowess(y_data, x_data, frac=frac, return_sorted=True)
        x_fit = loess_fit[:, 0]
        y_fit = loess_fit[:, 1]

        # Interpolate at pSyn1
        y1 = np.interp(1.489829144232, x_fit, y_fit)

        results.append({"Gene": name, "initial_s": y1})

    return pd.DataFrame(results)

y1_df = compute_y_at_one(updated_dfs)

length_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig1/Gene_length.csv"
Gene_length = pd.read_csv(length_file, sep=';')

data_bias = pd.merge(y1_df, Gene_length, on='Gene')

#Spearman correlation
X = data_bias["Length_nuc"]
y = data_bias["initial_s"]

res = stats.spearmanr(X, y)

print(res.statistic)
print(res.pvalue)
# r(9) = -0.691, p = 0.019

# Calculate linear regression parameters
slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)

# Create the regression line function
def regression_line(x_val):
    return slope * x_val + intercept

# Generate y-values for the regression line
regression_y = list(map(regression_line, X))

# Plotting
ax = sns.regplot(x="Length_nuc", y="initial_s", data=data_bias, scatter_kws={'alpha':0.7}, line_kws={'color':'black'}, ci=None)
ax.text(0.1, 0.2,
    r"$r$(9) = -0.691, p = 0.019",
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment='top'
)
ax.margins(x=0.1, y=0.1)

# Sort by gene length
df_sorted = data_bias.sort_values("Length_nuc").reset_index(drop=True)

for i, row in df_sorted.iterrows():

    if i < 5:
        dx = 5
        dy = 7
        ha = "left"
    else:
        dx = 5
        dy = -7
        ha = "left"

    ax.annotate(
        row["Gene"],
        (row["Length_nuc"], row["initial_s"]),
        xytext=(dx, dy),
        textcoords="offset points",
        ha=ha,
        va="center",
        fontsize=9,
        fontstyle="italic"
    )

plt.xlabel('Gene length (nucleotides)')
plt.ylabel('Initial selection coefficient')
plt.tight_layout()
plt.show()
