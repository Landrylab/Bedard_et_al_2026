"""
Fig sup 4: Baseline correction across different conditions for ESP1 and GCR1
"""
#%% Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
import seaborn as sns

#%% Prepare dataframes
def load_dfs_to_dict(folder_path):
    dfs_dict = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            key = os.path.splitext(filename)[0]
            filepath = os.path.join(folder_path, filename)
            dfs_dict[key] = pd.read_csv(filepath)
    return dfs_dict

processed_data = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Processed_data/Processed_data_Conditions"
updated_dfs = load_dfs_to_dict(processed_data)

# Add the expression values from deboer
DeBoer_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/DeBoer_Expression.csv"
Expression_file = pd.read_csv(DeBoer_file, sep=';')

expression_df = Expression_file[['Seq', 'Expression']]

# Merge 'Expression' into each DataFrame in the dictionary
for key in updated_dfs:
    updated_dfs[key] = updated_dfs[key].merge(expression_df, on='Seq', how='left')
    
order = ["ESP1_Glucose50%", "ESP1_Galactose", "ESP1_Raffinose",
         "ESP1_DMSO", "ESP1_Benomyl", "ESP1_MMS", "ESP1_GlucoseAgit",
         "GCR1_Glucose50%", "GCR1_Galactose", "GCR1_Raffinose",
         "GCR1_DMSO", "GCR1_Benomyl", "GCR1_MMS", "GCR1_GlucoseAgit"]

updated_dfs = {k: updated_dfs[k] for k in order if k in updated_dfs}

#%%
esp1_dfs = {k: v for k, v in updated_dfs.items() if k.startswith("ESP1")}
gcr1_dfs = {k: v for k, v in updated_dfs.items() if k.startswith("GCR1")}

#%% Compute LOESS + CI for each dataset
def compute_loess(dfs_dict, frac=0.4, n_bootstrap=200, ci=0.95):

    loess_results = {}
    rng = np.random.default_rng()

    for name, df in dfs_dict.items():

        # Clean data
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=["Expression", "s_T0_T7", "Replicate"])
        df_comp = df[df["Type"] == "Comp"]

        x_data = df_comp["Expression"].values
        y_data = df_comp["s_T0_T7"].values

        if len(x_data) == 0:
            continue

        # Sort
        order = np.argsort(x_data)
        x_sorted = x_data[order]
        y_sorted = y_data[order]

        # Loess fit
        loess_fit = lowess(y_sorted, x_sorted, frac=frac, return_sorted=True)
        x_fit = loess_fit[:, 0]
        y_fit = loess_fit[:, 1]

        # Bootstrap CI
        y_boots = np.zeros((n_bootstrap, len(x_fit)))

        for b in range(n_bootstrap):
            idx = rng.choice(len(x_data), size=len(x_data), replace=True)
            x_b = x_data[idx]
            y_b = y_data[idx]
            loess_b = lowess(y_b, x_b, frac=frac, xvals=x_fit) if hasattr(lowess, 'xvals') else lowess(y_b, x_b, frac=frac)[:,1]
            y_boots[b, :] = loess_b

        lower = np.percentile(y_boots, 100*(1-ci)/2, axis=0)
        upper = np.percentile(y_boots, 100*(1+ci)/2, axis=0)

        # Find max
        idx_max = np.argmax(y_fit)
        x_max = x_fit[idx_max]
        y_max = y_fit[idx_max]

        # Store results
        loess_results[name] = {
            "x_fit": x_fit,
            "y_fit": y_fit,
            "lower": lower,
            "upper": upper,
            "x_max": x_max,
            "y_max": y_max
        }

    return loess_results

esp1_loess = compute_loess(esp1_dfs, frac=0.4)
gcr1_loess = compute_loess(gcr1_dfs, frac=0.4)

#%% LOESS + CI per replicate with datapoints in grid
def compute_loess_per_replicate(dfs_dict, frac=0.4, n_bootstrap=200, ci=0.95):

    rng = np.random.default_rng()
    results = {}

    for name, df in dfs_dict.items():
        results[name] = {}
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=["Expression", "s_T0_T7", "Replicate"])
        df_comp = df[df["Type"] == "Comp"]

        replicate_ids = df_comp["Replicate"].unique()
        for rep in replicate_ids:
            df_rep = df_comp[df_comp["Replicate"] == rep]

            x_data = df_rep["Expression"].values
            y_data = df_rep["s_T0_T7"].values

            if len(x_data) == 0:
                continue

            # Sort
            order = np.argsort(x_data)
            x_sorted = x_data[order]
            y_sorted = y_data[order]

            # LOESS
            loess_fit = lowess(y_sorted, x_sorted, frac=frac, return_sorted=True)
            x_fit = loess_fit[:,0]
            y_fit = loess_fit[:,1]

            # Bootstrap CI
            y_boots = np.zeros((n_bootstrap, len(x_fit)))
            for b in range(n_bootstrap):
                idx = rng.choice(len(x_data), size=len(x_data), replace=True)
                x_b = x_data[idx]
                y_b = y_data[idx]
                loess_b = lowess(y_b, x_b, frac=frac, xvals=x_fit) if hasattr(lowess,'xvals') else lowess(y_b, x_b, frac=frac)[:,1]
                y_boots[b,:] = loess_b

            lower = np.percentile(y_boots, 100*(1-ci)/2, axis=0)
            upper = np.percentile(y_boots, 100*(1+ci)/2, axis=0)

            # Max
            idx_max = np.argmax(y_fit)
            x_max = x_fit[idx_max]
            y_max = y_fit[idx_max]


            results[name][rep] = {
                "x_fit": x_fit,
                "y_fit": y_fit,
                "lower": lower,
                "upper": upper,
                "x_max": x_max,
                "y_max": y_max,
                "x_data": x_data,
                "y_data": y_data
            }

    return results

esp1_loess_reps = compute_loess_per_replicate(esp1_dfs, frac=0.4, n_bootstrap=500)
gcr1_loess_reps = compute_loess_per_replicate(gcr1_dfs, frac=0.4, n_bootstrap=500)

#%% Extract baseline corrections
def extract_baselines(loess, loess_reps):

    rows = []

    # Baselines with all data
    for name, res in loess.items():

        baseline = res["y_fit"][0]

        # split dataset name
        if "_" in name:
            gene, condition = name.split("_", 1)
        else:
            gene = name
            condition = "NA"

        rows.append({
            "gene": gene,
            "condition": condition,
            "baseline": baseline,
            "source": "all",
            "replicate": "all"
        })

    # Baselines per replicate
    for name, reps_dict in loess_reps.items():

        if "_" in name:
            gene, condition = name.split("_", 1)
        else:
            gene = name
            condition = "NA"

        for rep_id, res in reps_dict.items():

            baseline = res["y_fit"][0]

            rows.append({
                "gene": gene,
                "condition": condition,
                "baseline": baseline,
                "source": "replicate",
                "replicate": rep_id
            })

    return pd.DataFrame(rows)

df_esp1 = extract_baselines(esp1_loess, esp1_loess_reps)

df_gcr1 = extract_baselines(gcr1_loess, gcr1_loess_reps)

baseline_df = pd.concat([df_esp1, df_gcr1], ignore_index=True)

#%% SupFig 4

plt.figure(figsize=(10,6))

# per replicate
rep_df = baseline_df[baseline_df["source"] == "replicate"]

sns.stripplot(
    data=rep_df,
    x="condition",
    y="baseline",
    hue="gene",
    dodge=False,
    jitter=True,
    alpha=0.25,
    size=8
)

# all data
pooled_df = baseline_df[baseline_df["source"] == "all"]

sns.scatterplot(
    data=pooled_df,
    x="condition",
    y="baseline",
    hue="gene",
    s=100,
    edgecolor="black",
    linewidth=1.2,
    legend=True
)

# clean legend
handles, labels = plt.gca().get_legend_handles_labels()

# keep only one entry per gene
unique = dict(zip(labels, handles))

plt.legend(
    unique.values(),
    unique.keys(),
    title="",
    fontsize=16
)

plt.xticks(rotation=0, fontsize=14)
plt.yticks(fontsize=14)

plt.ylabel("Initial bias", fontsize=16)
plt.xlabel("Condition", fontsize=16)

plt.tight_layout()
plt.show()
