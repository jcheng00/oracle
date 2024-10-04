import requests
import json
import pandas as pd
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description='Fetch orthologs using DIOPT API.')
parser.add_argument('input_species_id', type=int, help='Input species ID (e.g., 9606 for human)')
parser.add_argument('output_species_id', type=int, help='Output species ID (e.g., 7227 for fruit fly)')
parser.add_argument('entrez_id', type=str, help='Entrez ID (e.g., 51816 for ADA2)')

args = parser.parse_args()

input_species_id = args.input_species_id
output_species_id = args.output_species_id
entrez_id = args.entrez_id

url = f"https://www.flyrnai.org/tools/diopt/web/diopt_api/v8/get_orthologs_from_entrez/{input_species_id}/{entrez_id}/{output_species_id}/none"
req = requests.get(url, verify=False)
data = req.json()
gene_name = data["search_details"]["gene_details"][0]["symbol"]
df = pd.DataFrame(data["results"][entrez_id])

# Transpose the DataFrame to switch rows and columns
df = df.transpose()

# Reset the index to make it the first column
df.reset_index(inplace=True)

# Rename the first column to 'entrez_id'
df.rename(columns={'index': 'entrez_id'}, inplace=True)

# Reorder columns to make 'symbol' the second column
cols = ['entrez_id', 'symbol'] + [col for col in df.columns if col not in ['entrez_id', 'symbol']]
df = df[cols]

folder = "ortholog_output"
os.makedirs(folder, exist_ok=True)

df.to_csv(f"{folder}/{gene_name}_fly_orthologs.csv", index=False)
