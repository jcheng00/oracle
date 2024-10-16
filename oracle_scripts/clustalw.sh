#!/bin/bash
#SBATCH --job-name=msa_test1
#SBATCH --output=output_%j.log  # Standard output log (%j is replaced by job ID)
#SBATCH --error=error_%j.log  # Standard error log (%j is replaced by job ID)
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacie.cheng@nih.gov # Send-to address

module load clustalw

# Define the output folder
OUTPUT_FOLDER="ortholog_output"

# Check if the output folder exists, if not, create it
if [ ! -d "$OUTPUT_FOLDER" ]; then
  mkdir -p "$OUTPUT_FOLDER"
fi

# Get the input file from the command-line arguments
while getopts ":i:" opt; do
  case $opt in
    i) infile="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
  esac
done

# Check if infile is set
if [ -z "$infile" ]; then
  echo "Usage: $0 -i <input_file>"
  exit 1
fi

# Run clustalw with the specified input and output file
clustalw -INFILE="$infile" -ALIGN -OUTFILE="$OUTPUT_FOLDER/alignment.aln"
