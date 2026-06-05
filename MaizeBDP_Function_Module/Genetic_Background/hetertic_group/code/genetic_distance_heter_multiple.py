import argparse
import logging
from base_function_multiple import run_command, ibs_distance_matrix

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

def ibs_distance_heter(input_file, panel_prefix):
    """
    Calculate the p-distance between different samples.

    Args:
        input_file (str): Path to the input genotype data file.
    """
    # log
    log_file = "ibs_dist.log"
    setup_logger(log_file)
    
    # read data
    logging.info("Filtering user vcf...")
    command_filter = f"plink --vcf {input_file} --snps-only --biallelic-only --double-id --keep-allele-order --recode vcf-iid --out filter"
    command_bgzip = f"bgzip filter.vcf && tabix filter.vcf.gz"
    run_command(command_filter)
    run_command(command_bgzip)
    logging.info(f"Filter the genotype data and save to filter.vcf.gz")

    # read panel data
    logging.info("Combine user vcf and panel vcf...")
    command_combine_1 = f"bcftools isec filter.vcf.gz {panel_prefix}.panel.vcf.gz -n =2 -c all > common.site"
    command_combine_2 = f"bcftools merge -R common.site filter.vcf.gz {panel_prefix}.panel.vcf.gz -o merge.vcf"
    run_command(command_combine_1)
    run_command(command_combine_2)
    logging.info(f"Combine the user vcf and panel vcf and save to merge.vcf")
    
    # calculate IBS distance with plink
    logging.info("Calculate IBS distance with plink...")
    ibs_matrix = ibs_distance_matrix('merge.vcf')
    return ibs_matrix


def main():
    parser = argparse.ArgumentParser(description="Calculate p-distance between different samples")
    parser.add_argument("-i", "--input_file", type=str, help="Path to the input genotype data file")
    parser.add_argument("-p", "--panel_prefix", type=str, help="Path to the panel prefix")
    args = parser.parse_args()

    # Calculate ibs and write to file
    dist_file = ibs_distance_heter(args.input_file, args.panel_prefix)
    dist_file.to_csv('ibs.txt', sep='\t', index=True, header=True, float_format='%.6f')

if __name__ == "__main__":
    main()