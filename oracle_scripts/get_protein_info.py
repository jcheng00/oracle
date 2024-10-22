import sys
import io
import os
import logging
import argparse
from typing import List
from zipfile import ZipFile

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
            with open(args.output_file_name, "w") as output_file:
                output_file.write(fh.read())
except KeyError as e:
    logger.error("File %s not found in zipfile: %s", "protein.faa", e)
except FileNotFoundError as e:
    logger.error("Zipfile %s not found: %s", args.zipfile_name, e)