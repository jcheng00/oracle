import sys
import io
import os
import logging
import argparse
from typing import List
from zipfile import ZipFile
from Bio import SeqIO
from collections import defaultdict

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi.api.gene_api import GeneApi as DatasetsGeneApi

# Set up logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Download and extract protein sequences.")
parser.add_argument("gene_ids", type=int, nargs='+', help="List of gene IDs to download")
parser.add_argument("zipfile_name", type=str, help="Name of the zip file to save the dataset")
parser.add_argument("output_file_name", type=str, help="Name of the output file to save protein sequences")
args = parser.parse_args()

# download the data package using the DatasetsGeneApi, and then print out protein sequences for A2M and GNAS
with DatasetsApiClient() as api_client:
    gene_api = DatasetsGeneApi(api_client)
    try:
        gene_dataset_download = gene_api.download_gene_package(
            args.gene_ids,
            include_annotation_type=["FASTA_GENE", "FASTA_PROTEIN"],
        )
        with open(args.zipfile_name, "wb") as f:
            f.write(gene_dataset_download.read())
    except DatasetsApiException as e:
        sys.exit(f"Exception when calling GeneApi: {e}\n")

try:
    with ZipFile(args.zipfile_name) as dataset_zip:
        zinfo = dataset_zip.getinfo("ncbi_dataset/data/protein.faa")
        with io.TextIOWrapper(dataset_zip.open(zinfo), encoding="utf8") as fh:
            # Dictionary to store the longest sequence for each gene ID
            longest_sequences = defaultdict(lambda: None)

            for record in SeqIO.parse(fh, "fasta"):
                # Extract gene ID from the record description or ID
                # This assumes the gene ID is part of the record's description or ID
                # You may need to adjust this based on the actual format
                gene_id = record.id.split('|')[0]  # Example: adjust this line as needed

                # Update the longest sequence for the gene ID
                if longest_sequences[gene_id] is None or len(record.seq) > len(longest_sequences[gene_id].seq):
                    longest_sequences[gene_id] = record

            # Write the longest sequences to the output file
            with open(args.output_file_name, "w") as output_file:
                SeqIO.write(longest_sequences.values(), output_file, "fasta")
except KeyError as e:
    logger.error("File %s not found in zipfile: %s", "protein.faa", e)
except FileNotFoundError as e:
    logger.error("Zipfile %s not found: %s", args.zipfile_name, e)