"""
Supplementary figure 2: Correlation between replicates for selection coefficients measured in SC glucose
"""
#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

#%%
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

#%% SupFig 2

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

    # Clean df
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
            if gene_name == "WBP1":
                # Only the bottom-left subplot keeps spines and ticks
                if i == 2 and j == 0:
                    ax.spines[['top', 'right']].set_visible(False)
                    ax.tick_params(
                        left=True,
                        bottom=True,
                        labelleft=True,
                        labelbottom=True
                        )
                else:
                    for spine in ax.spines.values():
                        spine.set_visible(False)
                            
                    ax.tick_params(
                        left=False,
                        bottom=False,
                        labelleft=False,
                        labelbottom=False
                        )

            else:
                # Normal behavior for all other genes
                for spine in ax.spines.values():
                    if i > j:
                        ax.spines[['top', 'right']].set_visible(False)
                    else:
                        spine.set_visible(False)

                ax.tick_params(
                    left=False,
                    bottom=False,
                    labelleft=False,
                    labelbottom=False
                    )

                if i > j and j == 0:
                    ax.tick_params(left=True, labelleft=True)

                if i > j and i == 2:
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
            style="italic",
            y=1
        )

        title_ax.axis("off")


genes = list(updated_dfs.keys())

fig = plt.figure(figsize=(20, 15))
outer_gs = GridSpec(3, 4, figure=fig, wspace=0.2, hspace=0.2)

for idx, gene in enumerate(genes):
    row = idx // 4
    col = idx % 4

    plot_replicate_matrix_into_gs(
        updated_dfs[gene],
        fig,
        outer_gs[row, col],
        gene_name=gene,
        min_n=100
    )

# Legend in the empty panel
legend_ax = fig.add_subplot(outer_gs[2, 3])
legend_ax.axis("off")

# Collect all unique Type values across all dataframes
all_types = sorted(
    set(
        t
        for df in updated_dfs.values()
        for t in df["Type"].dropna().unique()
    )
)

palette = sns.color_palette(n_colors=len(all_types))
handles = [mpatches.Patch(color=palette[i], label=type_name)
           for i, type_name in enumerate(all_types)]

legend_ax.legend(
    handles=handles,
    title="Type",
    loc="center",
    frameon=False,
    fontsize=20,
    title_fontsize=24
)

plt.show()

#fig.savefig("/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig2/SupFig2.png", dpi=300, bbox_inches="tight", format="png")
#fig.savefig("/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig2/SupFig2.svg", bbox_inches="tight", format="svg")
