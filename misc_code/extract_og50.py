import pandas as pd

# Path to the TSV file
file_path = "C:/Users/victo/Documents/1000spider/shiny_app/filtered_with_ids_no_dupes.tsv"

# Load the dataset
data = pd.read_csv(file_path, sep="\t", index_col=0)

# Extract gene IDs for the orthogroup OG0000050 and the species Larinioides_sclopetarius
orthogroup = "OG0000050"
species = "Larinioides_sclopetarius.fasta"  # Column name for the species

if orthogroup in data.index and species in data.columns:
    # Get the cell for the specific orthogroup and species, drop NaN if present
    gene_ids = data.loc[orthogroup, species]
    if pd.notna(gene_ids):
        gene_ids_list = gene_ids.split(", ")  # Split multiple IDs in a cell
        print(f"Gene IDs for {orthogroup} in {species}:")
        print(gene_ids_list)
    else:
        print(f"No gene IDs found for {orthogroup} in {species}.")
else:
    if orthogroup not in data.index:
        print(f"Orthogroup {orthogroup} not found in the dataset.")
    if species not in data.columns:
        print(f"Species {species} not found in the dataset.")
