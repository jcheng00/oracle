import re
import argparse

def clean_text_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    cleaned_content = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return cleaned_content

def split_on_transition(s):
    return re.findall(r'[A-Za-z]+|\d+', s)

def replace_letter_at_index(s, index, new_letter, original_letter):
    if index < 0 or index >= len(s):
        raise ValueError("Index is out of bounds")
    if s[index] != original_letter:
        raise ValueError(f"Original letter at position {index+1} is not '{original_letter}'")
    return s[:index] + new_letter + s[index + 1:]

def make_mutant_fasta_header(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    mutant_fasta_header = ''
    for line in lines:
        if line.startswith('>'):
            mutant_fasta_header = line.replace('>', '>MUTANT_', 1)
            mutant_fasta_header = mutant_fasta_header.strip()
    return mutant_fasta_header

def write_strings_to_file(file_path, first_string, second_string):
    with open(file_path, 'w') as file:
        # Write the first string as the first line
        file.write(first_string + '\n')
        
        # Break the second string into lines of 50 characters
        for i in range(0, len(second_string), 50):
            file.write(second_string[i:i+50] + '\n')

def main(input_protein_file, mutations_of_interest):
    """
    Processes a protein file and applies specified mutations.

    Parameters:
    input_protein_file (str): The path to the input protein file.
    mutations_of_interest (list of str): A list of mutations to apply.

    Example command to run the script:
    python oracle/scripts/make_mutation_fasta.py ADA2.txt G47A Y453C
    """
    # Import protein text
    protein = clean_text_file(input_protein_file)
    mutant_protein = clean_text_file(input_protein_file)

    # Parse desired mutations
    mutation_list = []
    for mutation in mutations_of_interest:
        result = split_on_transition(mutation)
        mutation_list.append(result)

    # Create mutations
    for mutation in mutation_list:
        original_residue = mutation[0]
        position = int(mutation[1])
        new_residue = mutation[2]

        mutant_protein = replace_letter_at_index(mutant_protein, position-1, new_residue, original_residue)  # -1 for zero-indexing

    # Make output file
    output_protein_file = input_protein_file.split(".")[0] + "_mutant.txt"
    new_fasta_header = make_mutant_fasta_header(input_protein_file)
    write_strings_to_file(output_protein_file, new_fasta_header, mutant_protein)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a protein file and mutations.")
    parser.add_argument('input_protein_file', type=str, help='The input protein file')
    parser.add_argument('mutations_of_interest', type=str, nargs='+', help='List of mutations of interest')
    args = parser.parse_args()
    
    main(args.input_protein_file, args.mutations_of_interest)
