import pandas as pd
import numpy as np
from itertools import combinations
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
    

def simulation_and_compare(pattern_data, hybrid_data):
    """
    Simulate the genotype of the recommended variety pattern offspring and calculate the distance from the hybrid

    pattern_data (str): The variety pattern genotype
    hybrid_data (str): The hybrid genotype
    """
    # log
    log_file = "simulation.log"
    setup_logger(log_file)
    logging.info("Starting reading vcf...")
    
    # read vcf
    if pattern_data is None:
        logging.error("Error: The user genotype is None.")
        return None
    if hybrid_data is None:
        logging.error("Error: The hybrid data is None.")
        return None
    
    # read hybrid genotype
    command_filter_hybrid = f"plink --vcf {hybrid_data} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter_hybrid"
    command_bgzip_hybrid = f"bgzip filter_hybrid.vcf && tabix filter_hybrid.vcf.gz"
    command_hybird_sample = f"bcftools query -l filter_hybrid.vcf.gz > hybrid_sample.txt"
    run_command(command_filter_hybrid)
    run_command(command_bgzip_hybrid)
    run_command(command_hybird_sample)
    logging.info(f"Filter the hybrid genotype")
    hybrid_name = pd.read_csv('hybrid_sample.txt', header=None, sep='\t')
    hybrid_name = hybrid_name[0].tolist()
    
    
    # read parent genotype
    pattern_geno = pretreatment_out(pattern_data, 0, 0, 0)
    pattern_columns = pattern_geno.columns.tolist()
    
    # parent name list
    pattern_name = pattern_columns[4:]
    # create all group info
    combinations_list = list(combinations(pattern_name, 2))
    pattern_df = pd.DataFrame(combinations_list, columns=['parent1', 'parent2'])
    # loop
    result_list = []
    
    for index, row in pattern_df.iterrows():
        parent_id1 = row[0]
        parent_id2 = row[1]
        combine_geno = pd.concat([pattern_geno.iloc[:, 0:4], pattern_geno[[parent_id1, parent_id2]]], axis=1)
        combine_geno.columns = ['ID', 'Chr', 'Pos', 'Ref', 'Parent1', 'Parent2']
        # process the genotype
        # delete NA
        combine_geno = combine_geno.dropna()
        # delete the heterozygous loci
        def is_homozygous(genotype):
            return genotype[0] == genotype[1]
        combine_geno['parent1_homo'] = combine_geno.iloc[:, 4].apply(is_homozygous)
        combine_geno['parent2_homo'] = combine_geno.iloc[:, 5].apply(is_homozygous)
        combine_geno_homo = combine_geno[(combine_geno['parent1_homo']) & (combine_geno['parent2_homo'])]
        combine_geno_homo.drop(columns=['parent1_homo', 'parent2_homo'], inplace=True)
        # generate the offspring genotype by parents
        def generate_offspring_geno(parent1, parent2, ref):
            if parent1 == parent2:
                offspring_geno = parent1[0] + parent2[0]
            else:
                if parent1[0] == ref:
                    offspring_geno = parent1[0] + parent2[0]
                else:
                    offspring_geno = parent2[0] + parent1[0]
            return offspring_geno
        
        combine_geno_homo['Offspring'] = combine_geno_homo.apply(lambda row: generate_offspring_geno(row[4], row[5], row[3]), axis=1)
        combine_geno_homo.drop(columns=['Parent1', 'Parent2'], inplace=True)
        combine_geno_homo.to_csv('temp_combine_geno_homo.txt', sep='\t', header=True, index=False)
        chip2vcf_temp('temp_combine_geno_homo.txt', 'temp_combine_geno_homo.vcf')
        command_bgzip_vcf = f"bgzip temp_combine_geno_homo.vcf && tabix temp_combine_geno_homo.vcf.gz"
        run_command(command_bgzip_vcf)
        
        # combine homo and offspring
        command_combine_1 = f"bcftools isec temp_combine_geno_homo.vcf.gz filter_hybrid.vcf.gz -n =2 -c all > temp_common.site"
        command_conbine_2 = f"bcftools merge -R temp_common.site temp_combine_geno_homo.vcf.gz filter_hybrid.vcf.gz > temp_combine_geno_homo_hybrid.vcf"
        run_command(command_combine_1)
        run_command(command_conbine_2)
        
        temp_ibs_df = ibs_distance_matrix('temp_combine_geno_homo_hybrid.vcf')
        for hybrid in hybrid_name:
            result_list.append({
                'parent1': parent_id1,
                'parent2': parent_id2,
                'control_hybrid': hybrid,
                'distance_simulation_hybrid': temp_ibs_df.loc['Offspring', hybrid]
            })
        # delete the temp file
        command_delete_temp = f"rm -rf temp_combine_geno_homo.txt temp_combine_geno_homo.vcf.gz temp_combine_geno_homo.vcf.gz.tbi temp_combine_geno_homo_hybrid.vcf temp_common.site"
        run_command(command_delete_temp)
    
    result_df = pd.DataFrame(result_list)
    result_df.to_csv('output.simulation.info.txt', sep='\t', header=True, index=False, float_format='%.6f')

def main():
    parser = argparse.ArgumentParser(description="Simulate the genotype and calculate the distance with the hybrid")
    parser.add_argument("-p", "--pattern_geno", type=str, help="Path to the input pattern genotype data file")
    parser.add_argument("-b", "--hybrid_geno", type=str, help="Path to the input hybrid genotype data file")
    args = parser.parse_args()

    simulation_and_compare(args.pattern_geno, args.hybrid_geno)

if __name__ == "__main__":
    main()