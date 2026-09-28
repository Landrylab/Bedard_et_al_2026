"""
Process NGS data from the competition assays of ESP1 and GCR1 in different conditions.
"""

#%% 
import pandas as pd
from Bio import SeqIO
import numpy as np
import os

#%%
def parse_fasta_sizes(fasta_file):
    """Parse the FASTA file to extract sequence sizes."""
    size_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        seq_size = int(header.split("size=")[-1])
        size_dict[str(record.seq)] = seq_size
    return size_dict

def add_size_to_dataframe(df, size_dict, seq_col):
    """Add size information to the dataframe based on promoter sequence matches."""
    df['size'] = df[seq_col].map(size_dict).fillna(0)  # Fill missing with 0 or NaN
    return df

#%% Loading files

sample_plan_file_Screen_Conditions = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Plan_Screen_Conditions_agg.csv"
sample_plan_df_Screen_Conditions = pd.read_csv(sample_plan_file_Screen_Conditions, sep=';')

fasta_base_path_Screen_Conditions = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Data/Screen_2025_11/Aggregated_reads_2025_11/" 

promoter_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Prom_list.csv"
df_promoters = pd.read_csv(promoter_file, sep=';')

#%%
merged_dfs_Screen_Conditions = {}

# Loop through each sample in the plan
for _, row in sample_plan_df_Screen_Conditions.iterrows():
    sample_name = row["Sample"]
    timepoint = row["Timepoint"]
    gene = row["Gene"]
    replicate = row["Replicate"]
    media = row["Media"]
    fasta_file = fasta_base_path_Screen_Conditions + row["Fasta"]

    size_dict = parse_fasta_sizes(fasta_file)
    df_temp = add_size_to_dataframe(df_promoters.copy(), size_dict, "Seq")
    df_temp.rename(columns={"size": f"size_{timepoint}"}, inplace=True)

    # Calculate the total size and percentage
    total_size = df_temp[f"size_{timepoint}"].sum()
    df_temp[f"{timepoint}_percentage"] = df_temp[f"size_{timepoint}"] / total_size

   # Add metadata columns
    df_temp['Replicate'] = int(replicate)
    df_temp['Gene'] = gene
    df_temp['Media'] = media

    # Use a single string key "Gene_Media_Replicate" for merging
    merge_key = f"{gene}_{media}_{replicate}"

    if merge_key not in merged_dfs_Screen_Conditions:
        merged_dfs_Screen_Conditions[merge_key] = df_temp
    else:
        merged_dfs_Screen_Conditions[merge_key] = pd.merge(
            merged_dfs_Screen_Conditions[merge_key],
            df_temp,
            on=["Type", "Prom_Nb", "Seq"],
            how="outer"
        )

# Merge all replicates for the same Gene x Media combination
final_merged_dfs_Screen_Conditions = {}

for merge_key, df in merged_dfs_Screen_Conditions.items():
    
    gene, media, _ = merge_key.split('_')
    key = f"{gene}_{media}"

    if key not in final_merged_dfs_Screen_Conditions:
        final_merged_dfs_Screen_Conditions[key] = df
    else:
        final_merged_dfs_Screen_Conditions[key] = pd.concat([final_merged_dfs_Screen_Conditions[key], df], ignore_index=True)

# Filter out rows with less than 100 reads at T0
for key, df in final_merged_dfs_Screen_Conditions.items():
    if 'size_T0' in df.columns:
        final_merged_dfs_Screen_Conditions[key] = df[df['size_T0'] >= 100]
        
#%% Calculate selection coefficient between T0 and each timepoint

Generation_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Gentime_Screen_Conditions.csv"
Gen_file = pd.read_csv(Generation_file, sep=';')

Gen_file['Gene'] = Gen_file['Gene'].astype(str)
Gen_file['Replicate'] = Gen_file['Replicate'].astype(int)
Gen_file['Timepoint'] = Gen_file['Timepoint'].str.replace("T", "").astype(int)
Gen_file['Generation'] = pd.to_numeric(Gen_file['Generation'], errors='coerce')

updated_dfs = {}

timepoints = [0, 5, 7]

for df_name, df in final_merged_dfs_Screen_Conditions.items():

    Gene = str(df["Gene"].iloc[0])
    Media = df["Media"].iloc[0]

    for Replicate in df["Replicate"].unique():

        for i in range(len(timepoints)):
            for j in range(i + 1, len(timepoints)):

                t1 = timepoints[i]
                t2 = timepoints[j]

                # Get generations
                Gen_1_query = Gen_file[
                    (Gen_file['Gene'] == Gene) &
                    (Gen_file['Replicate'] == Replicate) &
                    (Gen_file['Timepoint'] == t1)
                ]

                Gen_2_query = Gen_file[
                    (Gen_file['Gene'] == Gene) &
                    (Gen_file['Replicate'] == Replicate) &
                    (Gen_file['Timepoint'] == t2)
                ]

                if Gen_1_query.empty or Gen_2_query.empty:
                    continue

                Gen_1 = Gen_1_query['Generation'].iloc[0]
                Gen_2 = Gen_2_query['Generation'].iloc[0]

                Generation = Gen_2 - Gen_1

                if Generation == 0:
                    continue

                # WT normalization
                WT_subset = df[
                    (df["Type"] == "Ctrl") &
                    (df["Replicate"] == Replicate)
                ]

                WT_FC_median = (
                    WT_subset[f"T{t2}_percentage"] /
                    WT_subset[f"T{t1}_percentage"]
                ).median()

                # log2FC
                df[f"log2FC_T{t1}_T{t2}"] = np.log2(
                    df[f"T{t2}_percentage"] / df[f"T{t1}_percentage"]
                )

                # Selection coefficient
                df[f"s_T{t1}_T{t2}"] = np.log2(
                    (df[f"T{t2}_percentage"] / df[f"T{t1}_percentage"]) /
                    WT_FC_median
                ) / Generation

    updated_dfs[df_name] = df

for df_name, df in final_merged_dfs_Screen_Conditions.items():

    # Extract Gene/Media from dataframe
    Gene = str(df["Gene"].iloc[0])
    Media = df["Media"].iloc[0]

    for Replicate in df["Replicate"].unique():

        for Timepoint in range(1, 8):

            Gen_i_query = Gen_file[
                (Gen_file['Gene'] == Gene) &
                (Gen_file['Replicate'] == Replicate) &
                (Gen_file['Timepoint'] == Timepoint)
            ]

            Gen_0_query = Gen_file[
                (Gen_file['Gene'] == Gene) &
                (Gen_file['Replicate'] == Replicate) &
                (Gen_file['Timepoint'] == 0)
            ]

            if Gen_i_query.empty or Gen_0_query.empty:
                continue

            Gen_i = Gen_i_query['Generation'].iloc[0]
            Gen_0 = Gen_0_query['Generation'].iloc[0]

            Generation = Gen_i - Gen_0

            # WT normalization
            WT_i = df[
                (df["Type"] == "Ctrl") &
                (df["Replicate"] == Replicate)
            ]

            WT_FC_median = (
                WT_i[f"T{Timepoint}_percentage"] /
                WT_i["T0_percentage"]
            ).median()

            # Calculate log2FC
            df[f"log2FC_T{Timepoint}"] = np.log2(
                df[f"T{Timepoint}_percentage"] / df["T0_percentage"]
            )

            # Selection coefficient
            df[f"s_T0_T{Timepoint}"] = np.log2(
                (df[f"T{Timepoint}_percentage"] / df["T0_percentage"]) /
                WT_FC_median
            ) / Generation

    updated_dfs[df_name] = df

#%%
def save_dict_of_dfs(dict_of_dfs, folder_path):
    os.makedirs(folder_path, exist_ok=True)
    for key, df in dict_of_dfs.items():
        filename = f"{key}.csv"
        filepath = os.path.join(folder_path, filename)
        df.to_csv(filepath, index=False)

processed_data = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Processed_data/Processed_data_Conditions/"
save_dict_of_dfs(updated_dfs, processed_data)
