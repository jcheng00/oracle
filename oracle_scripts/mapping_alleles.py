import requests
import os
import pandas as pd
import re
import argparse

def extract_numbers(s):
    """
    Extracts all numbers from a given string and returns them as a list of integers.

    Parameters:
    s (str): The input string from which to extract numbers.

    Returns:
    list: A list of integers extracted from the string.
    """
    # Use regular expression to find all sequences of digits in the string
    numbers = re.findall(r'\d+', s)
    
    # Convert the extracted sequences to integers
    return [int(num) for num in numbers]

def create_color_dict(df, key_col, value_col):
    """
    Creates a dictionary where the color is set to red if the string contains 'pathogenic' or 'Pathogenic'.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing the key and value columns.
    key_col (str): The name of the column to be used as keys in the dictionary.
    value_col (str): The name of the column to be checked for 'pathogenic' or 'Pathogenic'.

    Returns:
    dict: A dictionary with keys from key_col and values as 'red' if the value_col contains 'pathogenic' or 'Pathogenic'.
    """
    color_dict = {
    row[key_col]: 'red' if 'pathogenic' in row[value_col].lower() else 
    'green' if 'benign' in row[value_col].lower() else None
    for _, row in df.iterrows()
    }
    return color_dict

def generate_pymol_script(df, position_col, color_col, output_file):
    """
    Generates a PyMOL script to color-code amino acid positions based on a key.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing amino acid positions and color keys.
    position_col (str): The name of the column containing amino acid positions.
    color_col (str): The name of the column containing color keys.
    output_file (str): The path to the output PyMOL script file.
    """
    with open(output_file, 'w') as f:
        f.write("from pymol import cmd\n\n")
        
        for index, row in df.iterrows():
            position = row[position_col]
            color = row[color_col]
            f.write(f"cmd.color('{color}', 'resi {position}')\n")
        
        f.write("\ncmd.show('cartoon')\n")
        f.write("cmd.bg_color('white')\n")
        f.write("cmd.zoom()\n")

def map_known_alleles(gene_id):
    """
    Processes gene data for the given gene symbol.

    Parameters:
    gene_id (str): The gene symbol to process.
    """
    url = "http://v1.marrvel.org/data/clinvar"
    req = requests.get(url, params={"geneSymbol": gene_id})
    df = pd.read_json(req.text)

    # Filter the DataFrame to include only rows where the title contains "(p."
    filtered_df = df[df['title'].str.contains(r'\(p\.', na=False)]

    # Reset the index without keeping the old index column
    filtered_df = filtered_df.reset_index(drop=True)

    # Extract the string between parentheses that contains "p."
    filtered_df['protein_change'] = filtered_df['title'].str.extract(r'\(([^)]*p\.[^)]*)\)')

    # Remove the "p." prefix from the extracted protein change
    filtered_df['protein_change'] = filtered_df['protein_change'].str.replace('p.', '', regex=False)

    significance_description = []
    for i, row in filtered_df.iterrows():
        desc = filtered_df["significance"][i]["description"]
        significance_description.append(desc)
    filtered_df["significance_description"] = significance_description

    amino_acid_position = []
    for i in filtered_df["protein_change"].to_list():
        position = extract_numbers(i)
        if len(position) > 1:
            raise ValueError("Something went wrong. There should not be more than one amino acid position")
        amino_acid_position.append(position[0])
    filtered_df["amino_acid_position"] = amino_acid_position

    # Generate PyMOL script
    color_dict = create_color_dict(filtered_df, 'amino_acid_position', 'significance_description')
    filtered_df["color"] = filtered_df["amino_acid_position"].map(color_dict)
    generate_pymol_script(filtered_df, "amino_acid_position", "color", "color_code_script.pml")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process a gene symbol.")
    parser.add_argument('gene_id', type=str, help='The gene symbol to process')
    args = parser.parse_args()
    
    map_known_alleles(args.gene_id)
