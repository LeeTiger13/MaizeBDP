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


def simulation_and_compare(input_data, panel_prefix, threshold):
    """
    Recommended variety patterns based on genetic distance
    Simulate the genotype of the recommended variety pattern offspring and calculate the distance from the hybrid

    input_data (str): The input user's inbred genotype
    threshold (float): The threshold to judge the similar variety
    """
    # log
    log_file = "pattern_and_simulation.log"
    setup_logger(log_file)
    logging.info("Starting reading vcf...")
    
    if not isinstance(threshold, float):
        logging.error("Error: The parameter 'threshold' must be of type float.")
        return None

    if not 0 <= threshold <= 1:
        logging.error("Error: The parameter 'threshold' must be between 0 and 1.")
        return None

    logging.info("Starting variety pattern recommendation...")
    if input_data is not None:
        logging.info("Reading user vcf...")
        # filter user genotype
        command_filter_user = f"plink --vcf {input_data} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter_user"
        command_bgzip_user = f"bgzip filter_user.vcf && tabix filter_user.vcf.gz"
        run_command(command_filter_user)
        run_command(command_bgzip_user)
        logging.info(f"Filter the genotype data.")
        
        # combine user data and system data
        command_combine_1 = f"bcftools isec filter_user.vcf.gz {panel_prefix}/{panel_prefix}.variety_pattern_genotype.vcf.gz -n =2 -c all > common.site"
        command_combine_2 = f"bcftools merge -R common.site filter_user.vcf.gz {panel_prefix}/{panel_prefix}.variety_pattern_genotype.vcf.gz -o merge.vcf"
        run_command(command_combine_1)
        run_command(command_combine_2)
        logging.info(f"Combine the user genotype and variety pattern genotype.")
        
        # calculate the genetic distance
        logging.info(f"Calculate the genetic distance.")
        command_get_userid = f"bcftools query -l filter_user.vcf.gz > user_id.txt"
        run_command(command_get_userid)
        user_id_df = pd.read_csv('user_id.txt', sep='\t', header=None)
        user_id = user_id_df[0].tolist()
        # read pattern info
        pattern_info = pd.read_csv(f'{panel_prefix}/{panel_prefix}.pattern_hybrid_info.txt', sep='\t', header=0)
        pattern_sample = pattern_info['parent1'].tolist() + pattern_info['parent2'].tolist()
        # calcualte 
        ibs_dis_df = ibs_distance_matrix('merge.vcf')
        ibs_dis_df_md = ibs_dis_df.loc[pattern_sample, user_id]

        # read pattern genotype df
        pattern_df = pretreatment_out('merge.vcf', 0, 0, 0)
        # modify last 4 column name
        #new_columns = ['Z58', 'Chang7-2', 'PH6WC', 'PH4CV']
        #pattern_df.columns = pattern_df.columns[:-4].tolist() + new_columns
        
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
            # add the name of the other parent and offspring to the data box after
            df_pattern_md = pd.DataFrame(columns=['user_variety', 'actual_pattern_parent1', 'Distance',
                                                  'actual_pattern_parent2', 'Recommendation_variety_pattern', 'Hybrid'])
            for index, row in df_pattern.iterrows():
                sample_id = row['user_variety']
                parent_id = row['pattern_parent1']

                # search a match in the breed pattern information
                match_parent1 = pattern_info[pattern_info['parent1'] == parent_id]
                match_parent2 = pattern_info[pattern_info['parent2'] == parent_id]
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
                        hybrid = row1['hybrid']

                        # add the info to new dataframe
                        new_row = pd.DataFrame({
                            'user_variety': [sample_id],
                            'actual_pattern_parent1': [parent1_id],
                            'Distance': [s_distance],
                            'actual_pattern_parent2': [parent2_id],
                            'Recommendation_variety_pattern': [recommendation],
                            'Hybrid': [hybrid]
                        })
                        df_pattern_md = pd.concat([df_pattern_md, new_row], ignore_index=True)
            
            # read the hybrid genotype and the user genotype
            result_list = []
            for index, row in df_pattern_md.iterrows():
                new_p1 = row[0]
                #print(new_p1)
                old_p1 = row[1]
                distance_1 = row[2]
                new_p2 = row[3]
                #print(new_p2)
                recommend = row[4]
                hybrid = row[5]
                # combine new_p1 and new_p2
                combine_geno = pd.concat([pattern_df.iloc[:, 0:4], pattern_df[[new_p1, new_p2]]], axis=1)
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
                command_combine_1 = f"bcftools isec temp_combine_geno_homo.vcf.gz {panel_prefix}/{panel_prefix}.{hybrid}.vcf.gz -n =2 -c all > temp_common.site"
                command_conbine_2 = f"bcftools merge -R temp_common.site temp_combine_geno_homo.vcf.gz {panel_prefix}/{panel_prefix}.{hybrid}.vcf.gz > temp_combine_geno_homo_hybrid.vcf"
                run_command(command_combine_1)
                run_command(command_conbine_2)

                temp_ibs_df = ibs_distance_matrix('temp_combine_geno_homo_hybrid.vcf')
                result_list.append({
                    'user_variety': new_p1,
                    'actual_pattern_parent1': old_p1,
                    'Distance': distance_1, 
                    'actual_pattern_parent2': new_p2,
                    'Recommendation_variety_pattern': recommend,
                    'Hybrid': hybrid,
                    'distance_simulation_hybrid': temp_ibs_df.loc['Offspring', hybrid]
                })
                # delete the temp file
                command_delete_temp = f"rm -rf temp_combine_geno_homo.txt temp_combine_geno_homo.vcf.gz temp_combine_geno_homo.vcf.gz.tbi temp_combine_geno_homo_hybrid.vcf temp_common.site"
                run_command(command_delete_temp)

            result_df = pd.DataFrame(result_list)
            result_df.to_csv('output.pattern_and_simulation.info.txt', sep='\t', header=True, index=False, float_format='%.6f')
            logging.info("The variety pattern recommendation and simulation are completed.")

def main():
    parser = argparse.ArgumentParser(description="Simulate the genotype and calculate the distance with the hybrid")
    parser.add_argument("-i", "--input_data", type=str, help="Path to the input inbred genotype data file")
    parser.add_argument("-p", "--panel_prefix", type=str, help="Path to the panel prefix")
    parser.add_argument("-t", "--threshold", type=float, help="The threshold to judge the similar variety")
    args = parser.parse_args()

    simulation_and_compare(args.input_data, args.panel_prefix, args.threshold)

if __name__ == "__main__":
    main()