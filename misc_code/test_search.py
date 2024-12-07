import pandas as pd

# Load the file
file_path = "C:/Users/victo/Documents/1000spider/shiny_app/filtered_with_ids_no_dupes.tsv"
data = pd.read_csv(file_path, sep="\t", index_col=0)  # Adjust index_col based on the file format

# Search for a specific entry
gene_id = "IASV01035240.1_1"  # Example gene ID

# Function to check if the gene_id exists in any cell
def search_gene_in_data(gene_id, data):
    for orthogroup, row in data.iterrows():
        for cell in row:
            if pd.notna(cell) and gene_id in cell.split(", "):  # Split the cell values by ', ' and check
                return orthogroup
    return None

# Search for the matching orthogroup
matching_orthogroup = search_gene_in_data(gene_id, data)

# Print the result
if matching_orthogroup:
    print(f"Gene {gene_id} belongs to orthogroup {matching_orthogroup}")
else:
    print(f"Gene {gene_id} not found in the dataset")

