import argparse
import logging
import pandas as pd
import json
import matplotlib.pyplot as plt
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

def pretreatment_output(input_file, sample_heter, sample_miss, marker_miss, maf):
    """
    Pretreat the genotype data uploaded by users.

    Args:
        input_file (str): Path to the input genotype data file.
        sample_miss (float): The missing rate for filtering samples.
        marker_miss (float): The missing rate for filtering markers.
        maf (float): The minor allele frequency for filtering markers.
    """

    log_file = "pretreatment.log"
    setup_logger(log_file)    

    if not all(isinstance(x, float) for x in [sample_miss, marker_miss]):
        logging.error("Error: All parameters must be of type float.")
        return None

    # check whether between 0 and 1
    if not all(0 <= x <= 1 for x in [sample_miss, marker_miss]):
        logging.error("Error: All parameters must be between 0 and 1.")
        return None

    logging.info("Starting genotype data pretreatment...")

    # Read genotype data
    logging.info("Reading vcf information...")
    command_info = f"plink --vcf {input_file} --double-id --freq --het --missing --out output"
    run_command(command_info)
    logging.info(f"Save sample heterozygous rate to file output.het")
    logging.info(f"Save sample missing rate to file output.imiss")
    logging.info(f"Save marker missing rate to file output.lmiss")
    logging.info(f"Save marker minor allele frequency to file output.frq")
    
    logging.info("Catch the high heterozygous rate samples...")
    het_data = pd.read_csv('output.het', sep='\s+', index_col=False)
    het_data['het_rate'] = 1 - (het_data['O(HOM)'] / het_data['N(NM)'])
    high_het_samples = het_data[het_data['het_rate'] > sample_heter]['IID'].tolist()
    remove_file = "high_het_samples_to_remove.txt"
    with open(remove_file, 'w') as f:
        for sample in high_het_samples:
            f.write(f"{sample} {sample}\n")
    logging.info(f"Catch {len(high_het_samples)} samples with heterozygous rate > {sample_heter}")
    logging.info(f"high heterozygous rate samples: {high_het_samples}")

    # Filter genotype data
    logging.info("Filtering genotype data...")
    command_filter = f"plink --vcf {input_file} \
                             --mind {sample_miss} \
                             --geno {marker_miss} \
                             --maf {maf} \
                             --remove {remove_file} \
                             --snps-only --double-id --recode vcf-iid --biallelic-only --keep-allele-order \
                             --out mind{sample_miss}.geno{marker_miss}.maf{maf}.heter{sample_heter}.filter"
    run_command(command_filter)
    logging.info(f"Filter the genotype data and save to mind{sample_miss}.geno{marker_miss}.maf{maf}.heter{sample_heter}.filter.vcf")
    logging.info(f"------------------------Data pretreatment done-----------------------------")
    logging.info(f"------------------------Begin data information display---------------------")
    logging.info(f"Sample missing rate")
    imiss_data = pd.read_csv('output.imiss', sep='\s+', index_col=False)
    imiss_sample = len(imiss_data)
    imiss_filtered = len(imiss_data[imiss_data['F_MISS'] > sample_miss])
    imiss_keep = imiss_sample - imiss_filtered
    logging.info(f"Sample heterorzygous rate")
    heter_sample = len(het_data)
    heter_filtered = len(high_het_samples)
    heter_keep = heter_sample - heter_filtered
    logging.info(f"Marker missing rate")
    lmiss_data = pd.read_csv('output.lmiss', sep='\s+', index_col=False)
    lmiss_marker = len(lmiss_data)
    lmiss_filtered = len(lmiss_data[lmiss_data['F_MISS'] > marker_miss])
    lmiss_keep = lmiss_marker - lmiss_filtered
    
    json_output_rows = []
    imiss_row = {
        "Type" : "Sample missing rate",
        "Threshold" : sample_miss,
        "Total number" : imiss_sample,
        "Filtered number" : imiss_filtered,
        "Kept number" : imiss_keep,
    }
    heter_row = {
        "Type" : "Sample heterozygous rate",
        "Threshold" : sample_heter,
        "Total number" : heter_sample,
        "Filtered number" : heter_filtered,
        "Kept number" : heter_keep,
    }
    lmiss_row = {
        "Type" : "Marker missing rate",
        "Threshold" : marker_miss,
        "Total number" : lmiss_marker,
        "Filtered number" : lmiss_filtered,
        "Kept number" : lmiss_keep,
    }
    json_output_rows.append(imiss_row)
    json_output_rows.append(heter_row)
    json_output_rows.append(lmiss_row)
    
    json_output = {
        "metadata" : {
            "title" : "Genotype Data Pretreatment Summary",
            "description" : "This table summarizes the filtering results based on sample missing rate, sample heterozygous rate, and marker missing rate.",
        },
        "columns" : [
            {"key" : "Type", "title" : "The pretreatment type"},
            {"key" : "Threshold", "title" : "The filtering threshold"},
            {"key" : "Total number", "title" : "The total number before filtering"},
            {"key" : "Filtered number", "title" : "The number of samples/markers filtered"},
            {"key" : "Kept number", "title" : "The number of samples/markers kept"},
        ],
        "data" : json_output_rows
    }
    with open('output_table.json', 'w') as json_file:
        json.dump(json_output, json_file, indent=2)
    logging.info(f"Save the filtering summary to output_table.json")
    
    logging.info(f"MAF output picture")
    maf_data = pd.read_csv('output.frq', sep='\s+', index_col=False)
    bins = [0, 0.001, 0.01, 0.05, 0.5]
    bin_labels = ['0-0.001', '0.001-0.01', '0.01-0.05', '>0.05']

    maf_categories = pd.cut(maf_data['MAF'], bins=bins, labels=bin_labels, right=False)
    counts = maf_categories.value_counts().sort_index()

    plt.figure(figsize=(12, 8))

    colors = ['#ff9999', '#ffcc99', '#99ccff', '#99ff99']
    bars = plt.bar(bin_labels, counts, color=colors, alpha=0.8, edgecolor='black', linewidth=1.2)
    for i, (bar, count) in enumerate(zip(bars, counts)):
        percentage = count / len(maf_data) * 100
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                f'{count:,}\n({percentage:.1f}%)', 
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    total_variants = len(maf_data)
    ultra_rare_count = counts[bin_labels[0]]
    rare_count = counts[bin_labels[1]]
    low_freq_count = counts[bin_labels[2]]
    common_count = counts[bin_labels[3]]

    stats_text = f'''Total markers: {total_variants:,}
    MAF<0.001: {ultra_rare_count:,} ({ultra_rare_count/total_variants*100:.2f}%)
    0.001≤MAF<0.01: {rare_count:,} ({rare_count/total_variants*100:.2f}%)
    0.01≤MAF<0.05: {low_freq_count:,} ({low_freq_count/total_variants*100:.2f}%)
    MAF≥0.05: {common_count:,} ({common_count/total_variants*100:.2f}%)'''

    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, fontsize=9,
            bbox=dict(boxstyle='round,pad=0.5', alpha=0.9, facecolor='lightyellow'),
            verticalalignment='top')

    plt.xlabel('MAF', fontsize=12, fontweight='bold')
    plt.ylabel('Number of markers', fontsize=12, fontweight='bold')
    plt.title('Marker MAF Distribution Statistics', fontsize=14, fontweight='bold')

    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.xticks(fontsize=9)
    plt.ylim(0, max(counts) * 1.15)

    plt.tight_layout()

    plt.savefig('maf_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Save the MAF distribution figure to maf_distribution.png")
    


def main():
    parser = argparse.ArgumentParser(description="Pretreat the genotype data uploaded by users")
    parser.add_argument("-i", "--input_file", type=str, help="Path to the input genotype data file")
    parser.add_argument("-s", "--sample_heter", type=float, default=0.05, help="The heterozygous rate for filtering samples")
    parser.add_argument("-m", "--missing_rate", type=float, default=0.05, help="The missing rate for filtering samples")
    parser.add_argument("-n", "--marker_miss", type=float, default=0.2, help="The missing rate for filtering markers.")
    parser.add_argument("-a", "--maf", type=float, default=0.05, help="The minor allele frequency for filtering markers.")
    args = parser.parse_args()

    # Pretreat the genotype data and write to new file
    pretreatment_output(args.input_file, args.sample_heter, args.missing_rate, args.marker_miss, args.maf)

if __name__ == "__main__":
    main()