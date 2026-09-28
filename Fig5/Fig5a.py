# -*- coding: utf-8 -*-
"""
Fig5a: Fitness function of ESP1 in SC glucose
"""
#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
#%% Importing data

def load_dfs_to_dict(folder_path):
    dfs_dict = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            key = os.path.splitext(filename)[0]  # Remove .csv extension
            filepath = os.path.join(folder_path, filename)
            dfs_dict[key] = pd.read_csv(filepath)
    return dfs_dict

processed_data = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Processed_data/Processed_data_SC"
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

#%% Fig5a

gene = "ESP1"
df = updated_dfs[gene]

# Clean data
df = (
    df.replace([np.inf, -np.inf], np.nan)
      .dropna(subset=["Expression", "s_T0_T7", "Replicate"])
)

df_comp = df[df["Type"] == "Comp"]

x_data = df_comp["Expression"].values
y_data = df_comp["s_T0_T7"].values

# Figure
fig, ax = plt.subplots(figsize=(5,4))

# Raw points
ax.scatter(
    x_data,
    y_data,
    alpha=0.25,
    s=10
)

# LOWESS fit
frac = 0.4
n_bootstrap = 200
ci = 0.95

loess_fit = lowess(
    y_data,
    x_data,
    frac=frac,
    return_sorted=True
)

x_fit = loess_fit[:, 0]
y_fit = loess_fit[:, 1]

# Bootstrap CI
y_boots = np.zeros((n_bootstrap, len(x_fit)))

rng = np.random.default_rng()

for b in range(n_bootstrap):

    idx = rng.choice(
        len(x_data),
        size=len(x_data),
        replace=True
    )

    x_b = x_data[idx]
    y_b = y_data[idx]

    loess_b = lowess(
        y_b,
        x_b,
        frac=frac,
        return_sorted=True
    )

    y_boots[b, :] = np.interp(
        x_fit,
        loess_b[:, 0],
        loess_b[:, 1]
    )

# Confidence intervals
lower = np.percentile(
    y_boots,
    100 * (1 - ci) / 2,
    axis=0
)

upper = np.percentile(
    y_boots,
    100 * (1 + ci) / 2,
    axis=0
)

# Plot LOWESS + CI
ax.plot(
    x_fit,
    y_fit,
    color='tab:orange',
    linewidth=2
)

ax.fill_between(
    x_fit,
    lower,
    upper,
    color='tab:orange',
    alpha=0.2
)

# Formatting
ax.set_title(gene, fontsize=16, style="italic")
ax.set_xlabel("Promoter strength", fontsize=12)
ax.set_ylabel("Selection coefficient", fontsize=12)

ax.set_ylim(-0.08, 0.05)

ax.yaxis.grid(False)
ax.xaxis.grid(False)

ax.axvline(
    x= 11.3211665856487, #pSyn65 expression
    color='black',
    linestyle='--',
    linewidth=1.25
)

plt.tight_layout()
plt.show()