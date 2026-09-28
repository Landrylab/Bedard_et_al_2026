"""
SupFig 6: Correlation between replicates for selection coefficients measured in different condititions
(ESP1 and GCR1 only)
"""

#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter

#%%
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

order = ["ESP1_Glucose50%", "ESP1_Galactose", "ESP1_Raffinose",
         "ESP1_DMSO", "ESP1_Benomyl", "ESP1_MMS", "ESP1_GlucoseAgit",
         "GCR1_Glucose50%", "GCR1_Galactose", "GCR1_Raffinose",
         "GCR1_DMSO", "GCR1_Benomyl", "GCR1_MMS", "GCR1_GlucoseAgit"]

updated_dfs = {k: updated_dfs[k] for k in order if k in updated_dfs}

ESP1_dfs = {k: v for k, v in updated_dfs.items() if k.startswith("ESP1")}
GCR1_dfs = {k: v for k, v in updated_dfs.items() if k.startswith("GCR1")}

#%% Replicates correlability

# Function: replicate matrix
def plot_replicate_matrix_into_gs(
    df,
    fig,
    outer_gs,
    gene_name=None,
    min_n=100,
    palette=None,
    hue_order=None
):
    replicates = [1, 2, 3]

    # Clean data
    df = df.copy()
    df["Replicate"] = df["Replicate"].astype(int)
    df["Type"] = df["Type"].astype(str)

    # Pivot
    df_pivot = (
        df
        .pivot_table(
            index=["Prom_Nb", "Type"],
            columns="Replicate",
            values="s_T0_T7",
            aggfunc="mean"
        )
        .reset_index()
        .replace([np.inf, -np.inf], np.nan)
    )

    # Grid (3×3)
    inner_gs = GridSpecFromSubplotSpec(
        3, 3,
        subplot_spec=outer_gs,
        wspace=0.2,
        hspace=0.2
    )

    for i, rep_y in enumerate(replicates):
        for j, rep_x in enumerate(replicates):
            ax = fig.add_subplot(inner_gs[i, j])

            # Diagonal
            if i == j:
                ax.text(
                    0.5, 0.5, f"Rep. {rep_x}",
                    ha="center", va="center",
                    fontsize=15, fontweight="bold",
                    transform=ax.transAxes
                )

            # Scatter
            elif i > j:
                plot_df = df_pivot[[rep_x, rep_y, "Type"]].dropna()
                n = len(plot_df)

                if n < min_n:
                    ax.set_facecolor("lightgray")
                else:
                    sns.scatterplot(
                        data=plot_df,
                        x=rep_x,
                        y=rep_y,
                        hue="Type",
                        palette=palette,
                        hue_order=hue_order,
                        s=8,
                        alpha=0.6,
                        ax=ax,
                        legend=False
                    )

                    min_val = plot_df[[rep_x, rep_y]].min().min()
                    max_val = plot_df[[rep_x, rep_y]].max().max()

                    ax.set_xlim(min_val, max_val)
                    ax.set_ylim(min_val, max_val)
                    ax.axline((min_val, min_val), slope=1,
                              linestyle="--", linewidth=0.8)

            # Correlation
            else:
                corr_df = df_pivot[[rep_x, rep_y]].dropna()
                n = len(corr_df)

                if n < min_n:
                    ax.set_facecolor("lightgray")
                else:
                    r, p = pearsonr(corr_df[rep_x], corr_df[rep_y])
                    ax.text(
                        0.5, 0.5,
                        f"r = {r:.2f}\n(n = {n})",
                        ha="center", va="center",
                        fontsize=14,
                        transform=ax.transAxes
                    )

            # Clean axes
            for spine in ax.spines.values():
                if i > j:
                    ax.spines[['top', 'right']].set_visible(False)
                else:
                    spine.set_visible(False)
                
            ax.tick_params(
                left=False, bottom=False,
                labelleft=False, labelbottom=False
            )
            ax.set_xlabel("")
            ax.set_ylabel("")
            
            # Show ticks only for scatter plots
            if i > j:

                # y-axis only for first column
                if j == 0:
                    ax.tick_params(left=True, labelleft=True)

                # x-axis only for bottom row
                if i == 2:
                    ax.tick_params(bottom=True, labelbottom=True)
                    
            # Reduce tick label size and decimal precision
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=10
                )

            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            ax.set_xlabel("")
            ax.set_ylabel("")

    # Title
    if gene_name:

        # keep only text after first "_"
        short_name = gene_name.split("_", 1)[1] if "_" in gene_name else gene_name

        title_ax = fig.add_subplot(outer_gs)

        title_ax.set_title(
            short_name,
            fontsize=24,
            y=0.96
        )

        title_ax.axis("off")


# Main figure
grouped_dfs = {
    "ESP1": {k: v for k, v in updated_dfs.items() if k.startswith("ESP1")},
    "GCR1": {k: v for k, v in updated_dfs.items() if k.startswith("GCR1")}
}

ncols = 4
nrows_per_group = max(
    int(np.ceil(len(dfs) / ncols))
    for dfs in grouped_dfs.values()
)

fig = plt.figure(figsize=(5*ncols, 5*nrows_per_group*2))

# One row for ESP1, one for GCR1
outer_groups = GridSpec(
    2, 1,
    figure=fig,
    hspace=0.18
)

for group_idx, (group_name, dfs_dict) in enumerate(grouped_dfs.items()):

    # Consistent Type mapping
    all_types = sorted(
        set(
            str(t)
            for df in dfs_dict.values()
            for t in df["Type"].dropna().unique()
        )
    )

    palette_list = sns.color_palette("tab10", n_colors=len(all_types))
    type_to_color = dict(zip(all_types, palette_list))

    genes = list(dfs_dict.keys())

    n = len(genes)
    nrows = int(np.ceil(n / ncols))

    group_gs = GridSpecFromSubplotSpec(
        nrows,
        ncols,
        subplot_spec=outer_groups[group_idx],
        wspace=0.2,
        hspace=0.2
    )

    # Plot each gene
    for idx, gene in enumerate(genes):
        row = idx // ncols
        col = idx % ncols

        plot_replicate_matrix_into_gs(
            dfs_dict[gene],
            fig,
            group_gs[row, col],
            gene_name=gene,
            min_n=100,
            palette=type_to_color,
            hue_order=all_types
        )

    # Legend
    total_slots = nrows * ncols

    if n < total_slots:
        legend_idx = n
        row = legend_idx // ncols
        col = legend_idx % ncols

        legend_ax = fig.add_subplot(group_gs[row, col])
        legend_ax.axis("off")

        handles = [
            mpatches.Patch(color=type_to_color[t], label=t)
            for t in all_types
        ]

        legend_ax.legend(
            handles=handles,
            title="Type",
            loc="center",
            frameon=False,
            fontsize=20,
            title_fontsize=24
        )

    # Group title
    group_ax = fig.add_subplot(outer_groups[group_idx])
    group_ax.set_title(group_name, fontsize=28, style="italic", y=1.04)
    group_ax.axis("off")

plt.tight_layout()
plt.show()

#fig.savefig("/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig6/SupFig6.png", dpi=300, bbox_inches="tight", format="png")
#fig.savefig(""/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig6/SupFig6.svg", bbox_inches="tight", format="svg")
