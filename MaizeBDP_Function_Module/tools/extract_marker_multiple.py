import argparse
import pandas as pd
import re
import logging
from base_function_multiple import run_command

def setup_logger(log_file):
    """
    Set up the logger to write logs to a file.

    Args:
        log_file (str): Path to the log file.
    """
    logging.basicConfig(
        level=logging.INFO,  
        format='%(asctime)s - %(levelname)s - %(message)s',  
        handlers=[
            logging.FileHandler(log_file),  # write
            logging.StreamHandler()  # print out
        ]
    )

def extract_marker(input_file, region):
    """
    Extract the markers in the specified region.

    input_file (str):  Path to the input genotype data file.
    region (str): The region for extracting the marker.
    """
    # log
    log_file = "extract_marker.log"
    setup_logger(log_file)
    
    # read file
    command_compress = f"bgzip {input_file} && tabix {input_file}.gz"
    run_command(command_compress)
    command_getmarker = f"bcftools view -r {region} {input_file}.gz -Ov -o submarker.vcf"
    run_command(command_getmarker)
    logging.info(f"Extracted markers in the specified region")

def main():
    parser = argparse.ArgumentParser(description="Extract the markers in the specified region.")
    parser.add_argument("-i","--input_file", type=str, help="Path to the input genotype data file")
    parser.add_argument("-r","--region", type=str, help="The region for extracting the marker (eg. 1:10000-20000)")
    args = parser.parse_args()

    extract_marker(args.input_file, args.region)

if __name__ == "__main__":
    main()