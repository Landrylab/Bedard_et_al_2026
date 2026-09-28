"""
Process NGS data from the competition assays in SC glucose for 11 genes.
"""

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
    df['size'] = df[seq_col].map(size_dict).fillna(0)  # Fill missing with 0
    return df

#%% Loading files

#for first assay
sample_plan_file_Screen_SC_1 = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Plan_Screen_SC_1_agg.csv"
sample_plan_df_Screen_SC_1 = pd.read_csv(sample_plan_file_Screen_SC_1, sep=';')
 
fasta_base_path_Screen_SC_1 = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Data/NGS/Comp/Aggregated_reads/"

#For second assay
sample_plan_file_Screen_SC_2 = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Plan_Screen_SC_2_agg.csv" 
sample_plan_df_Screen_SC_2 = pd.read_csv(sample_plan_file_Screen_SC_2, sep=';')

fasta_base_path_Screen_SC_2 = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Data/NGS/Comp/Aggregated_reads_2/"

promoter_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Prom_list.csv"
df_promoters = pd.read_csv(promoter_file, sep=';')

#%%first assay
merged_dfs_Screen_SC_1 = {}

# Loop through each sample in the plan
for _, row in sample_plan_df_Screen_SC_1.iterrows():
    sample_name = row["Sample"]
    timepoint = row["Timepoint"]
    fasta_file = fasta_base_path_Screen_SC_1 + row["Fasta"]

    size_dict = parse_fasta_sizes(fasta_file)
    df_temp = add_size_to_dataframe(df_promoters.copy(), size_dict, "Seq")
    df_temp.rename(columns={"size": f"size_{timepoint}"}, inplace=True)

    # Percentage of reads
    total_size = df_temp[f"size_{timepoint}"].sum()
    df_temp[f"{timepoint}_percentage"] = df_temp[f"size_{timepoint}"] / total_size

    # Extract gene_replicate
    gene_part = sample_name.split('_')[0] + '_' + sample_name.split('_')[1]
    replicate = sample_name.split('_')[1]
    gene = sample_name.split('_')[0]

    # Store the updated DataFrame by gene_part
    if gene_part not in merged_dfs_Screen_SC_1:
        merged_dfs_Screen_SC_1[gene_part] = df_temp
    else:
        merged_dfs_Screen_SC_1[gene_part] = pd.merge(
            merged_dfs_Screen_SC_1[gene_part], 
            df_temp, 
            on=["Type", "Prom_Nb", "Seq"], 
            how="outer"
        )
        
    df_temp['Replicate'] = replicate
    df_temp['Replicate'] = df_temp['Replicate'].astype(int)
    df_temp['Gene'] = gene
    
final_merged_dfs_Screen_SC_1 = {}

for key, df in merged_dfs_Screen_SC_1.items():
    base_name = key.split('_')[0]  # Extract the gene
    
    if base_name not in final_merged_dfs_Screen_SC_1:
        final_merged_dfs_Screen_SC_1[base_name] = df
    else:
        final_merged_dfs_Screen_SC_1[base_name] = pd.concat([final_merged_dfs_Screen_SC_1[base_name], df], ignore_index=True)

for key, df in final_merged_dfs_Screen_SC_1.items():
    final_merged_dfs_Screen_SC_1[key] = df.drop(df[~(df['size_T0'] >= 100)].index) #Remove the lines where there is less than 100 reads at T0

#%%second assay
merged_dfs_Screen_SC_2 = {}

# Loop through each sample
for _, row in sample_plan_df_Screen_SC_2.iterrows():
    sample_name = row["Sample"]
    timepoint = row["Timepoint"]
    fasta_file = fasta_base_path_Screen_SC_2 + row["Fasta"]

    size_dict = parse_fasta_sizes(fasta_file)
    df_temp = add_size_to_dataframe(df_promoters.copy(), size_dict, "Seq")
    df_temp.rename(columns={"size": f"size_{timepoint}"}, inplace=True)

    #Percentage reads
    total_size = df_temp[f"size_{timepoint}"].sum()
    df_temp[f"{timepoint}_percentage"] = df_temp[f"size_{timepoint}"] / total_size

    # Extract gene and replicate
    gene, replicate = sample_name.split('_')[:2]
    gene_part = f"{gene}_{replicate}"

    # Add Replicate and Gene columns
    df_temp['Replicate'] = int(replicate)
    df_temp['Gene'] = gene

    # Store merged DataFrame
    if gene_part not in merged_dfs_Screen_SC_2:
        merged_dfs_Screen_SC_2[gene_part] = df_temp
    else:
        merged_dfs_Screen_SC_2[gene_part] = pd.merge(
            merged_dfs_Screen_SC_2[gene_part], 
            df_temp.drop(columns=['Replicate', 'Gene']),
            on=["Type", "Prom_Nb", "Seq"], 
            how="outer"
        )

# Merge all replicates per gene
final_merged_dfs_Screen_SC_2 = {}

for key, df in merged_dfs_Screen_SC_2.items():
    base_name = key.split('_')[0]  # Gene name
    if base_name not in final_merged_dfs_Screen_SC_2:
        final_merged_dfs_Screen_SC_2[base_name] = df
    else:
        final_merged_dfs_Screen_SC_2[base_name] = pd.concat([final_merged_dfs_Screen_SC_2[base_name], df], ignore_index=True)

# Final filtering step
for key, df in final_merged_dfs_Screen_SC_2.items():
    final_merged_dfs_Screen_SC_2[key] = df.drop(df[~(df['size_T0'] >= 100)].index)
    
#%%
final_merged_dfs_Screen_SC_2.update(final_merged_dfs_Screen_SC_1)

#%% Calculate selection coefficient between T0 and each timepoint

Generation_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/Gentime_Screens_SC.csv"
Gen_file = pd.read_csv(Generation_file, sep=';')

Gen_file['Gene'] = Gen_file['Gene'].astype(str)
Gen_file['Replicate'] = Gen_file['Replicate'].astype(int)
Gen_file['Timepoint'] = Gen_file['Timepoint'].astype(int)  # Ensure it's an integer
Gen_file['Generation'] = pd.to_numeric(Gen_file['Generation'], errors='coerce')

updated_dfs = {}

for df_name, df in final_merged_dfs_Screen_SC_2.items():
    for _, row in df.iterrows():
        for Timepoint in range(1, 8):  # Loop through Timepoints 1 to 7
            
            # Extract Gene and Replicate
            Gene = str(row["Gene"])
            Replicate = (row["Replicate"])  

            # Query Gen_file for the current and previous timepoint
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

            # Extract Generation values
            Gen_i = Gen_i_query['Generation'].iloc[0]
            Gen_0 = Gen_0_query['Generation'].iloc[0]

            # 
            WT_i = df[(df["Type"] == "Ctrl") & (df["Replicate"] == Replicate) ]
            
            WT_FC_median = (WT_i[f"T{Timepoint}_percentage"] / WT_i["T0_percentage"]).median()
            Generation = Gen_i - Gen_0

            df[f"log2FC_T{Timepoint}"] = np.log2(
                (df[f"T{Timepoint}_percentage"] / df["T0_percentage"]))

            # Calculate selection coefficient
            df[f"s_T0_T{Timepoint}"] = np.log2(
                (df[f"T{Timepoint}_percentage"] / df["T0_percentage"]) / WT_FC_median
            ) / Generation

    # Save the updated dataframe to the dictionary
    updated_dfs[df_name] = df
    
#%% Add the expression values from deboer
DeBoer_file = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Data_preparation/DeBoer_Expression.csv"
Expression_file = pd.read_csv(DeBoer_file, sep=';')

expression_df = Expression_file[['Seq', 'Expression']]

# Merge 'Expression' into each DataFrame in the dictionary
for key in updated_dfs:
    updated_dfs[key] = updated_dfs[key].merge(expression_df, on='Seq', how='left')

#%% Save the files as csv
def save_dict_of_dfs(dict_of_dfs, folder_path):
    os.makedirs(folder_path, exist_ok=True)
    for key, df in dict_of_dfs.items():
        filename = f"{key}.csv"
        filepath = os.path.join(folder_path, filename)
        df.to_csv(filepath, index=False)

processed_data = "/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fitness_Functions/Processed_data/Processed_data_SC"
save_dict_of_dfs(updated_dfs, processed_data)