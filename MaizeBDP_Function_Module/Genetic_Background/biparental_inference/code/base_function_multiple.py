import pandas as pd
import numpy as np
import subprocess
from collections import Counter


def run_command(command):
    """
    Run shell command and print
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    if error:
        print(f"Error: {error.decode('utf-8')}")
    else:
        print(output.decode('utf-8'))
        
### read vcf as input
def read_data_v2c(input_vcf, sample_miss, marker_miss, maf):
    
    """
    Convert the vcf file to chip data.
    input_vcf (str): Path to the input vcf file.
    sample_miss (float): The missing rate for filtering samples.
    marker_miss (float): The missing rate for filtering markers.
    maf (float): The minor allele frequency for filtering markers.
    :return chip data
    """
    if sample_miss == 0 or marker_miss == 0 or maf == 0:
        command_filter = f"plink --vcf {input_vcf} --keep-allele-order --snps-only --double-id --biallelic-only --recode vcf-iid --out filtered"
        run_command(command_filter)
    else:
        # filter vcf
        command_filter = f"plink --vcf {input_vcf} --mind {sample_miss} --geno {marker_miss} --maf {maf} --keep-allele-order --snps-only  --double-id --biallelic-only --recode vcf-iid --out filtered"
        run_command(command_filter)
 
    # convert to df
    samples = []
    data = []
    
    with open('filtered.vcf', 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Parse sample names from header
                samples = line.strip().split('\t')[9:]
                continue
            
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                chrom, pos, var_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
                
                # Skip multi-allelic sites
                if ',' in alt:
                    continue
                
                # Process genotypes
                row = {
                    'ID': var_id if var_id else '.',
                    'chrom': chrom,
                    'position': pos,
                    'ref': ref
                }
                
                for i, sample_gt in enumerate(fields[9:]):
                    gt = sample_gt.split(':')[0]  # Extract GT (e.g., "0/0")
                    
                    if gt == '0/0':
                        row[samples[i]] = ref + ref
                    elif gt == '0/1' or gt == '1/0':
                        row[samples[i]] = ref + alt
                    elif gt == '1/1':
                        row[samples[i]] = alt + alt
                    else:
                        row[samples[i]] = 'NA'
                
                data.append(row)
    
    return pd.DataFrame(data)

### impute with major allele
def major_imputation(input_data):
    """
    Imputate the chip data with major allele.
    input_data (str):  Path to the input genotype data file.
    :return imputed chip data
    """
    # read file
    data = input_data
    data_geno = data.iloc[:, 4:]

    def imputate_by_major(row):
        row = [item if isinstance(item, str) else "nan" for item in row]
        genotype = Counter([geno for geno in row if geno != 'nan'])
        most_common_genotype = genotype.most_common(1)[0][0]
        for i in range(len(row)):
            if row[i] == 'nan':
                row[i] = most_common_genotype

        return row

    # extract column name
    columns = data_geno.columns

    # imputation
    data_md = data_geno.apply(lambda row: pd.Series(imputate_by_major(list(row)), index=columns), axis=1)

    result_df = pd.concat([data.iloc[:, :4], data_md.iloc[:, :]], axis=1)
    


    return result_df

### data pretreatment
def pretreatment_out(input_vcf, sample_miss=0.05, marker_miss=0.2, maf=0.05):
    """
    Pretreat the genotype data uploaded by users.

    Args:
        input_vcf (str): Path to the input genotype data file.
        heter_rate: (float): The heterozygous rate for filtering samples.
    """

    # Read genotype data
    df = read_data_v2c(input_vcf, sample_miss, marker_miss, maf)

    if df is not None:
        # reset the index
        df = df.reset_index(drop=True)

        # rename
        col_list = ['ID', 'Chr', 'Pos', 'Ref']
        df.rename(columns=dict(zip(df.columns[:4], col_list)), inplace=True)
        
        # Write the output
        df = major_imputation(df) # imputation
        return df

### convert chip to vcf
def chip2vcf_temp(input_file, output):
    """
    Convert the chip data to vcf

    Args:
        input_file (str):  Path to the input chip data.
    """
    # Read genotype data
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Get the information of samples
    samples = lines[0].strip().split('\t')[4:]

    # Create a dataframe
    data_rows = []
    header_skipped = False
    for line in lines[0:]:
        if not header_skipped:
            header_skipped = True
            continue
        parts = line.strip().split('\t')
        id, chrom, pos, ref = parts[:4]
        alt = 'N'
        genotypes = []
        for genotype in parts[4:]:
            if genotype == 'NA' or 'del' in genotype or 'ins' in genotype:
                genotypes.append('./.')
            elif genotype[0] == genotype[1] == ref:
                genotypes.append('0/0')
            elif genotype[0] == genotype[1] != ref:
                alt = genotype[1]
                genotypes.append('1/1')
            else:
                alt = genotype[1]
                genotypes.append('0/1')

        # create a line
        data_row = [chrom, pos, id, ref, alt, '100', '.', 'AC=154;AF=0.987;AN=156', 'GT:AD:DP:GQ:PL'] + genotypes
        data_rows.append(data_row)

    # convert to dataframe
    df = pd.DataFrame(data_rows,
                      columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples)

    # function
    def process_chrom(chrom):
        if chrom.startswith('chr'):
            return chrom[3:]
        else:
            return chrom

    # apply
    df['#CHROM'] = df['#CHROM'].apply(process_chrom)


    # add the header
    vcf_header = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'
    ]

    # write to file
    with open(output, 'w') as o:
        o.write('\n'.join(vcf_header) + '\n')
        df.to_csv(o, sep='\t', index=False, lineterminator='\n')

# calculate IBS distance with plink
def ibs_distance_matrix(input_file):
    """
    Calculate the IBS distance matrix with plink.
    input_file (str): Path to the input genotype data file.
    return: IBS distance matrix
    """
    #command_filter = 
    command_ibs = f"plink --vcf {input_file} --double-id --keep-allele-order --distance 1-ibs square flat-missing --out ibs"
    run_command(command_ibs)
    
    ibs_matrix = pd.read_csv('ibs.mdist', sep='\t', header=None)
    id_data = pd.read_csv('ibs.mdist.id', sep='\t', header=None)
    
    ibs_matrix.columns = id_data.iloc[:, 0].tolist()
    ibs_matrix.index = id_data.iloc[:, 0].tolist()
    
    return ibs_matrix