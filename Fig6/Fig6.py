"""
Figure 6: Fitness functions across different condiditions for ESP1 and GCR1
Fig sup 4: Baseline correction across different conditions
"""
#%% Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
import scikit_posthocs as sp

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
#%% Fig 6a and 6b
def plot_loess_results(loess_results, title, ymin=-0.2, ymax=0.1, show_max=True):

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.tab10.colors

    baseline_dict = {}

    for i, (gene_name, res) in enumerate(loess_results.items()):
        color = colors[i % len(colors)]

        x_fit = res["x_fit"]
        y_fit = res["y_fit"]
        lower = res["lower"]
        upper = res["upper"]

        # Baseline correction
        baseline = y_fit[0]
        baseline_dict[gene_name] = baseline

        y_shifted = y_fit - baseline
        lower_shifted = lower - baseline
        upper_shifted = upper - baseline
        
        # Show max
        if show_max:
            idx_max = np.argmax(y_shifted)
            x_max = x_fit[idx_max]
            y_max = y_shifted[idx_max]

        # CI
        ax.fill_between(
            x_fit,
            lower_shifted,
            upper_shifted,
            color=color,
            alpha=0.2
        )
        
        name_legend = gene_name.split("_", 1)[-1]

        # Fitness function
        ax.plot(
            x_fit,
            y_shifted,
            color=color,
            label=name_legend
        )

        if show_max:
            ax.scatter(
                x_max,
                y_max,
                color=color,
                zorder=5
            )

    ax.set_xlabel("Promoter strength", fontsize=16)
    ax.set_ylabel("Corrected selection coefficient", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_title(title, fontsize=22, style="italic")

    ax.set_ylim(ymin, ymax)

    #ax.legend(loc='lower left', fontsize=13)
    ax.grid(True, alpha=0.5)

    plt.tight_layout()
    plt.show()

    return baseline_dict

baselines_esp1 = plot_loess_results(
    esp1_loess,
    "ESP1",
    ymin=-0.1,
    ymax=0.1
)

baselines_gcr1 = plot_loess_results(
    gcr1_loess,
    "GCR1",
    ymin=-0.15,
    ymax=0.025,
    show_max=False
)

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

#%% Extract replicate-level maxima for stats

def extract_maxima(loess_results):

    rows = []

    for dataset, reps_dict in loess_results.items():

        for rep_id, res in reps_dict.items():

            x_fit = res["x_fit"]
            y_fit = res["y_fit"]

            # Baseline correction
            y_shifted = y_fit - y_fit[0]

            idx_max = np.argmax(y_shifted)
            x_max = x_fit[idx_max]
            y_max = y_shifted[idx_max]

            rows.append({
                "gene": dataset,
                "replicate": rep_id,
                "x_max": x_max,
                "y_max": y_max
            })

    return pd.DataFrame(rows)

esp1_stats_df = extract_maxima(esp1_loess_reps)


#%% Dunn's test (x_max for optimal expression, y_max for highest s)

groups = [group["x_max"].values for _, group in esp1_stats_df.groupby("gene")]

dunn = sp.posthoc_dunn(
    esp1_stats_df,
    val_col="x_max",
    group_col="gene",
    p_adjust="fdr_bh"
)

mask = np.tril(np.ones(dunn.shape), k=-1).astype(bool)

dunn_filtered = dunn.where(mask)

# Long format
dunn_long = (
    dunn_filtered.stack()
    .reset_index()
    .rename(columns={"level_0": "group1", "level_1": "group2", 0: "p_value"})
)

# Keep only significant comparisons
significant = dunn_long[dunn_long["p_value"] < 0.05]

# Sort by p-value
significant = significant.sort_values("p_value")

print(significant)

#%% Fig 6c
def plot_maxima(loess_results):
    
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.tab10.colors

    summary = {}

    for i, (dataset, reps_dict) in enumerate(loess_results.items()):

        color = colors[i % len(colors)]

        label = dataset.split("_", 1)[1] if "_" in dataset else dataset

        xs = []
        ys = []

        for rep_id, res in reps_dict.items():

            x_fit = res["x_fit"]
            y_fit = res["y_fit"]

            # Baseline correction
            y_shifted = y_fit - y_fit[0]

            idx_max = np.argmax(y_shifted)
            x_max = x_fit[idx_max]
            y_max = y_shifted[idx_max]

            xs.append(x_max)
            ys.append(y_max)

            # Plot replicate
            ax.scatter(
                x_max,
                y_max,
                color=color,
                alpha=0.4,
                s=60
            )

        xs = np.array(xs)
        ys = np.array(ys)

        # Mean of replicates
        x_mean = np.mean(xs)
        y_mean = np.mean(ys)

        # 95% CI
        x_lower, x_upper = np.percentile(xs, [2.5, 97.5])
        y_lower, y_upper = np.percentile(ys, [2.5, 97.5])

        # Convert to asymmetric errorbars
        x_err = [[x_mean - x_lower], [x_upper - x_mean]]
        y_err = [[y_mean - y_lower], [y_upper - y_mean]]

        # Store
        summary[dataset] = {
            "x_mean": x_mean,
            "y_mean": y_mean,
            "x_ci": (x_lower, x_upper),
            "y_ci": (y_lower, y_upper)
        }

        # Plot mean + CI
        ax.errorbar(
            x_mean,
            y_mean,
            xerr=x_err,
            yerr=y_err,
            fmt='o',
            color=color,
            ecolor=color,
            elinewidth=2,
            capsize=4,
            markersize=9,
            markeredgecolor='black',
            label = label
        )

    # Formatting
    ax.set_xlabel("Optimal promoter strength", fontsize=16)
    ax.set_xlim(8, 12.5)
    ax.set_ylabel("Corrected maximal selection coefficient", fontsize=16)
    ax.set_title("ESP1", fontsize=22, style="italic")
    ax.set_ylim(0.0, 0.08)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, alpha=0.5)

    # Clean legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #ax.legend(by_label.values(), by_label.keys(), fontsize=12, loc='center left')

    plt.tight_layout()

    return summary, fig, ax

summary, fig, ax = plot_maxima(esp1_loess_reps)

y_min, y_max = ax.get_ylim()
x_min, x_max = ax.get_xlim()

dy = 0.05 * (y_max - y_min)   # vertical spacing
dx = 0.05 * (x_max - x_min)   # horizontal spacing

# Add comparisons

x_comparisons = [
    ("ESP1_Galactose", "ESP1_MMS", "*", 0.073),
    ("ESP1_Galactose", "ESP1_Raffinose", "*", 0.007),
]

for d1, d2, text, y in x_comparisons:

    if d1 not in summary or d2 not in summary:
        continue

    x1 = summary[d1]["x_mean"]
    x2 = summary[d2]["x_mean"]

    ax.plot([x1, x2], [y, y], color='black', lw=1.5, zorder=10)

    ax.text((x1+x2)/2, y + 0.001, text,
            ha='center', va='bottom', fontsize=12, fontweight="bold")
    
y_comparisons = [
    ("ESP1_MMS", "ESP1_Galactose", "**", 11.65),
    ("ESP1_MMS", "ESP1_Benomyl", "*", 11.9),
    ("ESP1_DMSO", "ESP1_Galactose", "*", 12.15),
]

for d1, d2, text, x in y_comparisons:

    if d1 not in summary or d2 not in summary:
        continue

    y1 = summary[d1]["y_mean"]
    y2 = summary[d2]["y_mean"]

    ax.plot([x, x], [y1, y2], color='black', lw=1.5, zorder=10)

    ax.text(
        x + 0.05,
        (y1 + y2) / 2,
        text,
        ha='left',
        va='center',
        fontsize=12,
        fontweight="bold"
    )

plt.show()