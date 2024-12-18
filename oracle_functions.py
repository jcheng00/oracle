import requests
import json
import pandas as pd
import os
from Bio import Phylo
import matplotlib.pyplot as plt
import re

#########################
#   GENERAL FUNCTIONS   #
#########################

def iterate_files_in_folder(folder_path):
    file_list = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_list.append(os.path.join(root, file))
    return file_list

#################
#   API CALLS   #
#################

def pull_diopt_orthologs(input_species_id, output_species_id, entrez_id, output_folder):
    '''
    Fetches orthologous protein data from the DIOPT API for a given Entrez gene ID and species pair, 
    processes the data into a pandas DataFrame, and saves it as a CSV file.

    Parameters:
    - input_species_id (str): The species ID for the input species.
    - output_species_id (str): The species ID for the output species.
    - entrez_id (str): The Entrez gene ID for which orthologs are to be fetched.

    Returns:
    - pd.DataFrame: A DataFrame containing the orthologous gene data, with 'entrez_id' and 'symbol' as the first two columns, followed by other data columns.
    - file_name (str): The name of the CSV file where the DataFrame is saved.

    The function performs the following steps:
    1. Constructs a URL to query the DIOPT API and sends a GET request to retrieve JSON response.
    2. Reorders dataframe for easier reading
    3. Creates a directory named 'ortholog_output' if it doesn't exist.
    4. Saves the DataFrame to a CSV file named after the gene symbol.
    5. Prints a success message if orthologs are found, or an error message if not.

    Note:
    - The function suppresses SSL verification warnings when making the API request.
    - If an error occurs during the process, a generic error message is printed.
    '''
    url = f"https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez/{input_species_id}/{str(entrez_id)}/{output_species_id}/none"
    req = requests.get(url, verify=False)
    data = req.json()
    gene_name = data["search_details"]["gene_details"][0]["symbol"]
    df = pd.DataFrame(data["results"][str(entrez_id)])

    # Transpose the DataFrame to switch rows and columns
    df = df.transpose()

    # Reset the index to make it the first column
    df.reset_index(inplace=True)

    # Rename the first column to 'entrez_id'
    df.rename(columns={'index': 'entrez_id'}, inplace=True)

    # Reorder columns to make 'symbol' the second column
    cols = ['entrez_id', 'symbol'] + [col for col in df.columns if col not in ['entrez_id', 'symbol']]
    df = df[cols]

    file_name = f"{gene_name}_fly_orthologs.csv"

    df.to_csv(f"{output_folder}/{file_name}", index=False)

    print("Found DIOPT orthologs")

    return df, file_name

#############################
#   ORTHOLOG AND ALIGNMENT  #
#############################
        
def filter_diopt_results(df, file_name, output_folder):
    '''
    Filters a DataFrame to include rows that are likely the best ortholog for a given protein/gene

    Parameters:
    - df (pd.DataFrame): The input DataFrame, typically the output from the `pull_diopt_orthologs` function,
      containing orthologous gene data with columns such as 'best_score', 'best_score_rev', and 'confidence'.

    Returns:
    - pd.DataFrame: A new DataFrame containing only the rows from the input DataFrame where:
      - 'best_score' is "Yes", or
      - 'best_score_rev' is "Yes", or
      - 'confidence' is either "high" or "moderate".

    The function iterates over each row of the input DataFrame and checks the specified conditions.
    If a row meets any of these conditions, it is added to the output DataFrame.
    '''
    output_df = pd.DataFrame()
    for _, row in df.iterrows():
        if row["best_score"] == "Yes" or row["best_score_rev"] == "Yes" or row["confidence"] in ["high", "moderate"]:
            output_df = pd.concat([output_df, pd.DataFrame([row])], ignore_index=True)
    output_file = f"filtered_{file_name}"
    output_df.to_csv(f"{output_folder}/{output_file}", index = False)
    return output_df, output_file

#################
#   EVOLUTION  #
################

def visualize_phylogenetic_tree(file_path, file_format='newick', output_file=None):
    '''
    Visualizes a phylogenetic tree from a file and either displays it or saves it to an output file.

    Parameters:
    - file_path (str): The path to the input file containing the phylogenetic tree. 
      The expected file format is a .dnd file, which is typically in Newick format.
    - file_format (str, optional): The format of the input file. Defaults to 'newick'.
    - output_file (str, optional): The path to save the visualized tree as an image file. 
      If not provided, the tree will be displayed on the screen.

    The function performs the following steps:
    1. Reads the phylogenetic tree from the specified file using the Biopython Phylo module.
    2. Creates a matplotlib figure and draws the tree.
    3. If an output file path is provided, saves the figure to the specified file.
    4. If no output file is specified, displays the tree using matplotlib's show function.

    Note:
    - The function uses the 'newick' format by default, which is common for .dnd files.
    - The figure size is set to 10x8 inches by default, but this can be adjusted as needed.
    '''
    # Read the tree from the file
    tree = Phylo.read(file_path, file_format)
    
    # Create a figure and draw the tree
    fig = plt.figure(figsize=(10, 8))  # You can adjust the size as needed

    if output_file:
        # Save the figure to a file
        plt.savefig(output_file)
        print(f"Tree saved to {output_file}")
    else:
        # Show the plot if no output file is specified
        plt.show()

########################
#   MUTATION EFFECTS   #
########################

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

def generate_pymol_script_alleles(df, position_col, color_col, output_file):
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

def generate_pymol_script_domains(df, output_file='color_sites.pml'):
    """
    Generates a PyMOL script to color active and binding sites.

    Parameters:
    - df (pd.DataFrame): DataFrame containing site information with columns 'type', 
      'location.start.value', and 'location.end.value'.
    - output_file (str): The name of the output PyMOL script file.

    # Example usage:
    generate_pymol_script(df_flattened)
    """
    with open(output_file, 'w') as f:
        f.write("from pymol import cmd\n")
        
        for index, row in df.iterrows():
            start = row['location.start.value']
            end = row['location.end.value']
            site_type = row['type']
            
            if site_type == 'Active site':
                color = 'yellow'
            elif site_type == 'Binding site':
                color = 'blue'
            else:
                continue  # Skip if the type is neither 'Active site' nor 'Binding site'
            
            # Write the PyMOL command to color the specified range
            f.write(f"cmd.select('site_{index}', 'resi {start}-{end}')\n")
            f.write(f"cmd.color('{color}', 'site_{index}')\n")
        
        f.write("cmd.show('sticks', 'site_*')\n")
    
    print(f"PyMOL script saved to {output_file}")