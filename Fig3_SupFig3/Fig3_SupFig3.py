"""
Grids of fitness functions (Fig3). The code can be modified to look at the control promoters (SupFig3)
"""

#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import (ttest_ind, shapiro, levene)
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

#%% t-test to verify if the selection coefficient varies between lowest and strongest promoters for each gene.
test_dfs = {}

for name, df in updated_dfs.items():

    # Keep only competition data
    comp = df[df["Type"] == "Comp"].copy()

    # Promoters that are actually present
    promoters = sorted(comp["Prom_Nb"].unique())

    # Five lowest and five highest promoter numbers
    lower_prom = promoters[:5]
    higher_prom = promoters[-5:]

    # Keep only those promoters
    filtered = comp[
        comp["Prom_Nb"].isin(lower_prom + higher_prom)
    ].copy()

    # Label groups
    filtered["Limit"] = np.where(
        filtered["Prom_Nb"].isin(lower_prom),
        "Lower",
        "Higher"
    )

    test_dfs[name] = filtered
    combined_df = pd.concat(test_dfs.values(), ignore_index=True)

    # Average biological replicates for each promoter
    combined_df = (
        combined_df
        .groupby(["Gene", "Prom_Nb", "Limit"], as_index=False)
        .agg(
        s_T0_T7=("s_T0_T7", "mean")
        )
)

# Run one t-test per gene

results = []

# Run one analysis per gene
for gene, group in combined_df.groupby("Gene"):
    
    lower = group[group["Limit"] == "Lower"]["s_T0_T7"].dropna()
    higher = group[group["Limit"] == "Higher"]["s_T0_T7"].dropna()
    
    # Skip incomplete groups
    if len(lower) < 2 or len(higher) < 2:
        continue

    # Shapiro normality tests
    shapiro_lower_p = shapiro(lower).pvalue if len(lower) >= 3 else np.nan
    shapiro_higher_p = shapiro(higher).pvalue if len(higher) >= 3 else np.nan
    
    # Equal variance test
    levene_p = levene(lower, higher).pvalue
    
    # Conditions
    normal_lower = shapiro_lower_p > 0.05 if not np.isnan(shapiro_lower_p) else False
    normal_higher = shapiro_higher_p > 0.05 if not np.isnan(shapiro_higher_p) else False
    equal_var = levene_p > 0.05
    
    ttest_ok = normal_lower and normal_higher

    # t-test
    t_stat, t_pvalue = ttest_ind(
        lower,
        higher,
        equal_var=equal_var,
        nan_policy="omit"
    )

    # Store results
    results.append({
        "gene": gene,
        
        "n_lower": len(lower),
        "n_higher": len(higher),
        
        "mean_lower": lower.mean(),
        "mean_higher": higher.mean(),
        
        "shapiro_lower_p": shapiro_lower_p,
        "shapiro_higher_p": shapiro_higher_p,
        "levene_p": levene_p,
        
        "normal_lower": normal_lower,
        "normal_higher": normal_higher,
        "equal_variance": equal_var,
        "ttest_conditions_ok": ttest_ok,
        
        "t_stat": t_stat,
        "ttest_pvalue": t_pvalue,
    })

# Results df
stats_results = pd.DataFrame(results)

print(stats_results)

#%% #%% Fig 3 (if using df["Type"] == "Comp") or Fig Sup 3 ( if using: df["Type"] == "Ctrl")

def plot_grid(dfs_dict, n_bootstrap=200, frac=0.4,
              ci=0.95, ymin=-0.15, ymax=0.05):

    n_plots = len(dfs_dict)
    n_cols = 3
    n_rows = int(np.ceil(n_plots / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4*n_cols, 3*n_rows),
        sharex=True,
        sharey=True
    )

    axes = axes.flatten()

    for ax in axes:
        ax.tick_params(axis='x', labelbottom=True)

    for i, (name, df) in enumerate(dfs_dict.items()):

        ax = axes[i]

        # ---- Clean data ----
        df = (
            df.replace([np.inf, -np.inf], np.nan)
              .dropna(subset=["Expression", "s_T0_T7", "Replicate"])
        )

        df_comp = df[df["Type"] == "Comp"] # Change here to see control "Ctrl" or competition "Comp" data

        x_data = df_comp["Expression"].values
        y_data = df_comp["s_T0_T7"].values

        # ---- RAW POINTS (faint) ----
        ax.scatter(
            x_data,
            y_data,
            alpha=0.25,
            s=10
        )

        # ---- LOESS FIT ----
        if len(x_data) > 0:

            loess_fit = lowess(
                y_data,
                x_data,
                frac=frac,
                return_sorted=True
            )

            x_fit = loess_fit[:, 0]
            y_fit = loess_fit[:, 1]

            # ---- Bootstrap for CI ----
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

                # Bootstrap LOWESS
                loess_b = lowess(
                    y_b,
                    x_b,
                    frac=frac,
                    return_sorted=True
                )

                # Interpolate onto original grid
                y_boots[b, :] = np.interp(
                    x_fit,
                    loess_b[:, 0],
                    loess_b[:, 1]
                )

            # ---- Confidence intervals ----
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

            # ---- Plot LOWESS + CI ----
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

        ax.set_title(name, size=16, style="italic")
        ax.set_ylim(ymin, ymax)
        ax.grid(True)

    # ---- Remove empty panels ----
    for j in range(n_plots, len(axes)):
        fig.delaxes(axes[j])

    fig.supxlabel("Promoter strength", fontsize=20)
    fig.supylabel("Selection coefficient", fontsize=20)

    plt.tight_layout()
    plt.show()


plot_grid(updated_dfs, ymin=-0.12, ymax=0.05) #0.05 for Comp, 0.06 for Ctrl
