import argparse
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

def extract_sample(input_file, sample_list):
    """
    Extract the genotype of the specified sample.

    input_file (str): Path to the input genotype data file.
    list (str): Path to the sample list.
    """
    # log
    log_file = "extract_sample.log"
    setup_logger(log_file)
    
    # read file
    command_compress = f"bgzip {input_file} && tabix {input_file}.gz"
    run_command(command_compress)
    command_getsample = f"bcftools view -S {sample_list} {input_file}.gz -Ov -o subsample.vcf"
    run_command(command_getsample)
    logging.info(f"Extracted genotype of the specified sample")



def main():
    parser = argparse.ArgumentParser(description="Extract the genotype of the specified sample")
    parser.add_argument("-i","--input_file", type=str, help="Path to the input genotype data file")
    parser.add_argument("-l","--list", type=str, help="Path to the sample list")
    args = parser.parse_args()

    extract_sample(args.input_file, args.list)

if __name__ == "__main__":
    main()