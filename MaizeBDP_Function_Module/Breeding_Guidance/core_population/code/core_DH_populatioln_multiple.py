import argparse
import logging
from base_function_multiple import *

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


def p_distance_redundancy(input_file):
    """
    Calculate the p-distance between different samples.

    Args:
        input_file (str): Path to the input genotype data file.
    """
    # log
    log_file = "core_DH_selection.log"
    setup_logger(log_file)
    logging.info("Start to calculate IBS...")
    
    command_filter = f"plink --vcf {input_file} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter"
    run_command(command_filter)
    logging.info(f"Filter the genotype data and save to filter.vcf")
    
    # calculate the ibs
    ibs_matrix = ibs_distance_matrix('filter.vcf')
    #print(ibs_matrix.head())
    sample_df = pd.read_csv('ibs.mdist.id', sep='\t', header=None)
    sample_names = sample_df.iloc[:, 0].tolist()
    chara_list = list(range(1, ibs_matrix.shape[0] + 1))
    ibs_matrix.insert(0, 'NAME', chara_list)
    ibs_matrix.insert(0, 'ID', sample_names)
    # Write to file
    ibs_matrix.to_csv('output.DH.dist.txt', sep=',', index=False, header=True, float_format='%.3f')

def main():
    parser = argparse.ArgumentParser(description="Calculate p-distance between different samples")
    parser.add_argument("-i", "--input_file", type=str, help="Path to the input genotype data file")
    args = parser.parse_args()

    # Calculate p-distance and write to file
    p_distance_redundancy(args.input_file)


if __name__ == "__main__":
    main()

