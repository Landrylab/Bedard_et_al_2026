"""
Fig sup 5: Uncorrected fitness functions of ESP1 and GCR1 measured in different experimental conditions for each replicate
"""
#%% Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
from matplotlib.gridspec import GridSpec

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

#%% SupFig 5

def plot_loess_replicates(
    loess_results,
    title,
    n_cols=3,
    ymin=None,
    ymax=None
):


    datasets = list(loess_results.keys())
    n_panels = len(datasets)
    n_rows = int(np.ceil(n_panels / n_cols))

    # Figure + GridSpec
    fig = plt.figure(figsize=(4*n_cols, 3*n_rows))
    gs = GridSpec(
    n_rows,
    n_cols,
    figure=fig,
    wspace=0.05,
    hspace=0.3    
)

    axes = []
    colors = plt.cm.tab10.colors

    # Create axes
    for i, ds_name in enumerate(datasets):

        row = i // n_cols
        col = i % n_cols

        # Center final panel
        if (
            n_panels % n_cols == 1 and
            i == n_panels - 1
        ):
            col = n_cols // 2

        # Shared axes
        if len(axes) == 0:
            ax = fig.add_subplot(gs[row, col])
        else:
            ax = fig.add_subplot(
                gs[row, col],
                sharex=axes[0],
                sharey=axes[0]
            )

        axes.append(ax)

        replicates = loess_results[ds_name]

        # Plot replicates
        for j, (rep_id, res) in enumerate(replicates.items()):

            color = colors[j % len(colors)]

            # CI shaded area
            ax.fill_between(
                res["x_fit"],
                res["lower"],
                res["upper"],
                color=color,
                alpha=0.08
            )

            # LOESS line
            ax.plot(
                res["x_fit"],
                res["y_fit"],
                color=color,
                linewidth=2
            )

            # Raw points
            ax.scatter(
                res["x_data"],
                res["y_data"],
                color=color,
                alpha=0.35,
                s=12,
                edgecolor='none'
            )

        # Panel titles
        short_name = (
            ds_name.split("_", 1)[1]
            if "_" in ds_name else ds_name
        )

        ax.set_title(short_name, fontsize=18, pad=4)

        # Show x tick labels
        ax.tick_params(axis='x', labelbottom=True)

        # Y-axis label
        last_row = n_rows - 1

        show_y = (
            col == 0 or
            row == last_row
        )

        if not show_y:
            ax.tick_params(
                axis='y',
                labelleft=False
            )

        ax.grid(True, alpha=0.3)

        if ymin is not None and ymax is not None:
            ax.set_ylim(ymin, ymax)

    # Labels
    fig.supxlabel("Promoter strength", fontsize=20, y=0.05)

    fig.supylabel("Selection coefficient", fontsize=20, x=0.05)
    
    fig.suptitle(title, fontsize=24, style="italic", y=0.95)

    #Layout
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.show()

plot_loess_replicates(esp1_loess_reps, "ESP1", ymin=-0.15, ymax=0.05)

plot_loess_replicates(gcr1_loess_reps, "GCR1", ymin=-0.2, ymax=0.025)
