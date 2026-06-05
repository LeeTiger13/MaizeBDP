import pandas as pd
import numpy as np
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

def parent_prefer(parent_file, offspring_file):
    """
    Determine whether the offspring are more inclined towards the mother or the father.

    parent_file: Path to the input parental line file.
    offspring_file: Path to the input offspring DH line file.
    output_file: Path to the output file of p_dist
    """
    # log
    log_file = "biparent_bias.log"
    setup_logger(log_file)
    logging.info("Start to calculate the genetic distance between parental and offspring lines.")

    # read file
    logging.info("Reading the parental and offspring files...")
    command_parent = f"plink --vcf {parent_file} --keep-allele-order --snps-only --biallelic-only --double-id --recode vcf-iid --out parent"
    command_parent_bgzip = f"bgzip parent.vcf && tabix parent.vcf.gz"
    command_paretn_id = f"bcftools query -l parent.vcf.gz > parent.id"
    command_offspring = f"plink --vcf {offspring_file} --keep-allele-order --snps-only --biallelic-only --double-id --recode vcf-iid --out offspring"
    command_offspring_bgzip = f"bgzip offspring.vcf && tabix offspring.vcf.gz"
    command_offspring_id = f"bcftools query -l offspring.vcf.gz > offspring.id"
    run_command(command_parent)
    run_command(command_parent_bgzip)
    run_command(command_paretn_id)
    run_command(command_offspring)
    run_command(command_offspring_bgzip)
    run_command(command_offspring_id)
    logging.info(f"Filter the parental and offspring files and save to parent.vcf.gz and offspring.vcf.gz")

    # combine
    logging.info("Combining the parental and offspring files...")
    command_combine_1 = f"bcftools isec offspring.vcf.gz parent.vcf.gz -n =2 -c all > common.site"
    command_combine_2 = f"bcftools merge -R common.site offspring.vcf.gz parent.vcf.gz -o merge.vcf"
    run_command(command_combine_1)
    run_command(command_combine_2)
    logging.info(f"Combine the parental and offspring files and save to merge.vcf")

    # extract sample name
    parent_df = pd.read_csv('parent.id', sep='\t', header=None)
    parent_id = parent_df.iloc[:, 0].tolist()
    offspring_df = pd.read_csv('offspring.id', sep='\t', header=None)
    dh_id = offspring_df.iloc[:, 0].tolist()

    ibs_distance_matrix_df = ibs_distance_matrix('merge.vcf')
    ibs_distance_matrix_df_md = ibs_distance_matrix_df.loc[dh_id, parent_id]

    # judge the parent
    ibs_distance_matrix_df_md['parent_prefer'] = ''

    for index, row in ibs_distance_matrix_df_md.iterrows():
        if row[0] < row[1]:
            ibs_distance_matrix_df_md.at[index, 'parent_prefer'] = parent_id[0]
        elif row[0] > row[1]:
            ibs_distance_matrix_df_md.at[index, 'parent_prefer'] = parent_id[1]
        else:
            ibs_distance_matrix_df_md.at[index, 'parent_prefer'] = 'Both'

    ibs_distance_matrix_df_md.to_csv('output.biparental.bias.txt', sep='\t', index=True, header=True, float_format='%.3f')

def main():
    parser = argparse.ArgumentParser(description="Comparing different samples with genome window.")
    parser.add_argument("-p", "--parent_file", help="Input the parent chip file.")
    parser.add_argument("-d", "--offspring_file", help="Input the offspring chip file.")
    args = parser.parse_args()

    parent_prefer(args.parent_file, args.offspring_file)

if __name__ == "__main__":
    main()



