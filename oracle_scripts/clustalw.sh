#!/bin/bash
#SBATCH --job-name=msa_test1
#SBATCH --output=output_%j.log  # Standard output log (%j is replaced by job ID)
#SBATCH --error=error_%j.log  # Standard error log (%j is replaced by job ID)
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacie.cheng@nih.gov # Send-to address

module load clustalw

# Get the input file from the command-line arguments
while getopts ":INFILE:" opt; do
  case $opt in
    INFILE) infile="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

clustalw -INFILE=$infile -ALIGN