import subprocess
import argparse
import pandas as pd
import numpy as np
import os
import logging
import matplotlib.pyplot as plt
from pathlib import Path
from base_function_multiple import run_command

ADMIXTURE_BIN = Path(
    os.environ.get("ADMIXTURE_BIN", "admixture")
).as_posix()

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

# log
log_file = "pop_structure.log"
setup_logger(log_file)


def admixture_p(study_vcf, panel_prefix):
    """
    Run admixture p
    
    study_vcf (str): Path to the study vcf file.
    ref_Q_file (str): Path to the reference chip file.[fixed by system]
    k (int): The number of clusters.[fixed by system]
    """
    
    if study_vcf is not None:
        logging.info("Reading user's vcf and filtering...")

        # filter user vcf
        command_filter = f"plink --vcf {study_vcf} --snps-only --biallelic-only --keep-allele-order --double-id --recode vcf-iid --out filter"
        command_bgzip = f"bgzip filter.vcf && tabix filter.vcf.gz"
        run_command(command_filter)
        run_command(command_bgzip)
        logging.info(f"Filter the genotype data.")
        

        # read reference panel and make sure the snp in user's vcf equal to panel
        command_get_ref = f"bcftools query -f '%CHROM\t%POS\n' {panel_prefix}/{panel_prefix}.panel.vcf.gz > reference.txt"
        run_command(command_get_ref)
        if not Path(f"reference.txt").exists():
            logging.error("Reference file does not exist.")
            
        command_get_study = f"bcftools view -R reference.txt filter.vcf.gz -Oz -o filter.md.vcf.gz" 
        run_command(command_get_study)
        if not Path("filter.md.vcf.gz").exists():
            logging.error("Study file does not exist.")
            
        command_tabix = f"tabix filter.md.vcf.gz"
        run_command(command_tabix)
        
        command_merge = f"bcftools merge {panel_prefix}/{panel_prefix}.panel.vcf.gz filter.md.vcf.gz --force-samples --missing-to-ref -o merge.vcf.gz -Oz"
        run_command(command_merge)
        if not Path("merge.vcf.gz").exists():
            logging.error("Merge file does not exist.")
            
        command_tabix_merge = f"tabix merge.vcf.gz"
        run_command(command_tabix_merge)
        
        command_get_study_id = f"bcftools query -l filter.vcf.gz > user_id.txt"
        run_command(command_get_study_id)
        if not Path("user_id.txt").exists():
            logging.error("User ID file does not exist.")
        
        command_get_study = f"bcftools view -S user_id.txt merge.vcf.gz -Ov -o study.vcf"
        run_command(command_get_study)
        if not Path("study.vcf").exists():
            logging.error("Study file does not exist.")
        
        logging.info(f"Extract the common snps between uesr's vcf and panel")
        
        # admixture
        logging.info("Analysing population structure with admixture...")
        # make sure the studt.8.P.in file is in the same directory
        # do study
        k_data = pd.read_csv('chip.version.file.txt', sep='\t', header=0)
        k = k_data.loc[k_data['Chip'] == panel_prefix, 'k_value'].iloc[0]
        command_makebed_study = f"plink --vcf study.vcf --vcf-half-call missing --out study"
        command_getid = f"cut -f 1 study.fam > sample_list.txt"
        command_cp_ref = f"cp {panel_prefix}/{panel_prefix}.{k}.P.in study.{k}.P.in"
        command_admix_P = f"{ADMIXTURE_BIN} -P study.bed {k}" # K depend on panel
        command_modify_Q = f"paste sample_list.txt study.{k}.Q > study.{k}.md.Q" # K depend on panel
        command_modify_Q_1 = f"sed -i 's/ /\t/g' study.{k}.md.Q" # K depend on panel
        run_command(command_makebed_study)
        run_command(command_getid)
        run_command(command_cp_ref)
        run_command(command_admix_P)
        run_command(command_modify_Q)
        run_command(command_modify_Q_1)
        logging.info(f"Admixture calculate finished.")
        return f"study.{k}.md.Q" # K depend on panel
    
    else:
        logging.error("Error: No input file provided.")
    

def visualization(study_Q_file, panel_prefix):
    """
    Visualization of the results of admixture p

    study_Q_file (str): Path to the study Q file.
    """
    # read files
    logging.info("Visualizing the results...")
    df = pd.read_csv(study_Q_file, sep='\t', header=None)
    cluster_info = pd.read_csv(f'{panel_prefix}/{panel_prefix}.reference_cluster_name.txt', sep='\t', header=None)

    # extract the cluster information
    cluster_info = cluster_info.iloc[:, 0]
    cluster_info = cluster_info.tolist()
    df.columns = cluster_info

    # write table
    df.to_csv('cluster.info.txt', sep='\t', index=False, header=True, float_format='%.6f')

    # calculate the sum of each column
    k = df.shape[1] - 1
    # add new column
    # The maximum proportion of each sample in the different population
    df['maxval'] = df.iloc[:, 1:(df.shape[1] + 1)].max(axis=1)
    # The column name of the largest proportion of each sample in different populations
    df['matchval'] = df.iloc[:, 1:(df.shape[1] + 1)].idxmax(axis=1)
    # Calculates the sum of maxval column values in each matchval group
    cluster_sums = df.groupby('matchval')['maxval'].sum()
    sorted_clusters = cluster_sums.sort_values(ascending=False).index.tolist()
    df['cluster_rank'] = df['matchval'].map(
        {cluster: rank for rank, cluster in enumerate(sorted_clusters, start=1)})

    # set the order for drawing 
    cluster_counts = df['matchval'].value_counts()
    cluster_counts_sorted_small = cluster_counts.sort_values(ascending=True).index.tolist()
    
    df = df.sort_values(by=['cluster_rank', 'maxval'], ascending=[True, False])

    # create the sample id list
    ids = df.iloc[:, 0].tolist()

    # set the color
    colors = ['#aec7e8', '#2ca02c', '#cf0087', '#e377c2', '#8c564b',
              '#7826B6', '#1f77b4', '#ff7f0e', '#ffbb78', '#98df8a',
              '#17becf', '#bcbd22', '#7f7f7f', '#ff9896', '#c5b0d5']
    
    # sst the bottom
    bottom = np.zeros(len(ids))

    # set the figure
    fig, ax = plt.subplots(figsize=(24, 2))

    # delete the background
    ax.tick_params(which='both', bottom=False, top=False, left=False, right=False, labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)

    # draw plot
    for i, cluster in enumerate(cluster_counts_sorted_small):
        plt.bar(ids, df[cluster], bottom=bottom, color=colors[i], width=1, align='edge')
        bottom += df[cluster]
    
    # x label
    plt.xticks(ids)
    # rotation
    plt.xticks(rotation=90)
    # modify the x_axis
    plt.tick_params(axis='x', labelsize=3)
    # modify the y_axis
    plt.annotate(f'K = {k}', xy=(0, 0.5), xytext=(-20, 0), textcoords='offset points',
                    va='center', ha='right',fontsize=22)
    plt.savefig('cluster.out.pdf', format='pdf', bbox_inches='tight', dpi=300)


def main():
    parser = argparse.ArgumentParser(description="Convert genotype files and perform population structure analysis.")
    parser.add_argument("-i", "--input_vcf", help="Input genotype vcf.")
    parser.add_argument("-p", "--panel_prefix", type=str, help="Path to the panel prefix")
    args = parser.parse_args()

    # step 1: make q file
    study_Q_data = admixture_p(args.input_vcf, args.panel_prefix)

    # step 2: visualization
    visualization(study_Q_data, args.panel_prefix)

if __name__ == "__main__":
    main()

