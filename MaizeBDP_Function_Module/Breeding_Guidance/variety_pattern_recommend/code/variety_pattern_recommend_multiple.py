import pandas as pd
import numpy as np
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

# variety pattern recommendation
def pattern_recommend(input_data, pattern_data, pattern_info, threshold):
    """
    Recommended variety patterns based on genetic distance

    input_data (str): The input user's inbred genotype
    pattern_data (str): The variety pattern genotype
    pattern_info (str): The information of variety pattern and hybrid
    threshold (float): The threshold to judge the similar variety
    """
    log_file = "pattern.log"
    setup_logger(log_file)
    
    if not isinstance(threshold, float):
        logging.error("Error: The parameter 'threshold' must be of type float.")
        return None

    if not 0 <= threshold <= 1:
        logging.error("Error: The parameter 'threshold' must be between 0 and 1.")
        return None

    logging.info("Starting variety pattern recommendation...")
    if input_data is not None:
        logging.info("Reading user vcf...")
        if pattern_data is not None:
            logging.info("Reading pattern vcf...")
            if pattern_info is not None:
                logging.info("Reading pattern info...")
                # filter vcfs
                command_filter_user = f"plink --vcf {input_data} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter_user"
                command_bgzip_user = f"bgzip filter_user.vcf && tabix filter_user.vcf.gz"
                run_command(command_filter_user)
                run_command(command_bgzip_user)
                logging.info(f"Filter the genotype data and save to filter_user.vcf.gz")

                command_filter_pattern = f"plink --vcf {pattern_data} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter_pattern"
                command_bgzip_pattern = f"bgzip filter_pattern.vcf && tabix filter_pattern.vcf.gz"
                run_command(command_filter_pattern)
                run_command(command_bgzip_pattern)
                logging.info(f"Filter the genotype data and save to filter_pattern.vcf.gz")
                
                # combine user and pattern
                command_combine_1 = f"bcftools isec filter_user.vcf.gz filter_pattern.vcf.gz -n =2 -c all > common.site"
                command_combine_2 = f"bcftools merge -R common.site filter_user.vcf.gz filter_pattern.vcf.gz -o merge.vcf"
                run_command(command_combine_1)
                run_command(command_combine_2)
                logging.info(f"Combine the user vcf and pattern vcf and save to merge.vcf")
                
                # calculate the genetic distance
                logging.info("Calculate the genetic distance...")
                command_get_userid = f"bcftools query -l filter_user.vcf.gz > user_id.txt"
                command_get_patternid = f"bcftools query -l filter_pattern.vcf.gz > pattern_id.txt"
                run_command(command_get_userid)
                run_command(command_get_patternid)
                user_id_df = pd.read_csv('user_id.txt', sep='\t', header=None)
                pattern_id_df = pd.read_csv('pattern_id.txt', sep='\t', header=None)
                uesr_id_list = user_id_df[0].tolist()
                pattern_id_list = pattern_id_df[0].tolist()
                # calculate the ibs matrix
                ibs_dis_df = ibs_distance_matrix('merge.vcf')
                #print(ibs_dis_df.head())
                ibs_dis_df_md = ibs_dis_df.loc[pattern_id_list, uesr_id_list]
                # write out the p_dist matrix
                ibs_dis_df_md.to_csv('genetic_distance.txt', sep='\t', index=True, header=True, float_format='%.6f')

                # filtering with the threshold  
                mask = ibs_dis_df_md < threshold
                positions = mask.stack()[mask.stack()].index.tolist()
                # Consolidate combinations that are smaller than the threshold in the new data box
                df_pattern = pd.DataFrame(positions, columns=['pattern_parent1', 'user_variety'])
                df_pattern = df_pattern.iloc[:, [1, 0]]
                if df_pattern.shape[0] < 1:
                    logging.error("The threshold setting is too small to recommend variety patterns.")
                else:
                    # Based on the information of one of the parents, match in the breed pattern information
                    # add the name of the other parent and Recommendation to the data box after
                    df_pattern_md = pd.DataFrame(columns=['user_variety', 'actual_pattern_parent1', 'Distance',
                                                          'actual_pattern_parent2', 'Recommendation_variety_pattern'])
                    for index, row in df_pattern.iterrows():
                        sample_id = row['user_variety']
                        parent_id = row['pattern_parent1']

                        # search a match in the breed pattern information
                        pattern_info_df = pd.read_csv(pattern_info, sep='\t', header=0)
                        match_parent1 = pattern_info_df[pattern_info_df['parent1'] == parent_id]
                        match_parent2 = pattern_info_df[pattern_info_df['parent2'] == parent_id]
                        # combine two match
                        matches = pd.concat([match_parent1, match_parent2])
                        if not matches.empty:
                            for _, row1 in matches.iterrows():
                                if row1['parent1'] == parent_id:
                                    parent1_id = row1['parent1']
                                    parent2_id = row1['parent2']
                                else:
                                    parent1_id = row1['parent2']
                                    parent2_id = row1['parent1']
                                
                                s_distance = ibs_dis_df_md.loc[parent1_id, sample_id]
                                recommendation = f'{sample_id} x {parent2_id}'
                                
                                # add the info to new dataframe
                                new_row = pd.DataFrame({
                                    'user_variety': [sample_id],
                                    'actual_pattern_parent1': [parent1_id],
                                    'Distance': [s_distance],
                                    'actual_pattern_parent2': [parent2_id],
                                    'Recommendation_variety_pattern': [recommendation] 
                                })
                                
                                df_pattern_md = pd.concat([df_pattern_md, new_row], ignore_index=True)
                    df_pattern_md.to_csv('recommand_info.txt', sep='\t', index=False, header=True)
                    logging.info(f"Variety pattern recommendation finished. The result is saved to 'recommand_info.txt'")
                        
            else:
                logging.error("Error: The pattern info is None.")
        else:
            logging.error("Error: The pattern vcf is None.")
    else:
        logging.error("Error: The user's genotype is None.")



def main():
    parser = argparse.ArgumentParser(description="Variety pattern recommendation")
    parser.add_argument("-i", "--input_data", type=str, help="Path to the input inbred genotype data file")
    parser.add_argument("-p", "--pattern_geno", type=str, help="Path to the input pattern genotype data file")
    parser.add_argument("-f", "--pattern_info", type=str,help="Path to the information of variety pattern and hybrid")
    parser.add_argument("-t", "--threshold", type=float, help="The threshold to judge the similar variety")
    args = parser.parse_args()

    pattern_recommend(args.input_data, args.pattern_geno, args.pattern_info, args.threshold)

if __name__ == "__main__":
    main()