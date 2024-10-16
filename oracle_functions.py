import requests
import json
import pandas as pd
import os

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

def pull_diopt_orthologs(input_species_id, output_species_id, entrez_id):
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
    try: 
        url = f"https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez/{input_species_id}/{entrez_id}/{output_species_id}/none"
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

        folder = "ortholog_output"
        os.makedirs(folder, exist_ok=True)

        file_name = f"{gene_name}_fly_orthologs.csv"

        df.to_csv(file_name, index=False)

        print("Found DIOPT orthologs")

        return df, file_name
    except:
        print(f"No DIOPT Orthologs found for {entrez_id}")

#############################
#   ORTHOLOG AND ALIGNMENT  #
#############################
        
def filter_diopt_results(df, file_name):
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
    output_df.to_csv(output_file, index = False)
    return output_df, output_file