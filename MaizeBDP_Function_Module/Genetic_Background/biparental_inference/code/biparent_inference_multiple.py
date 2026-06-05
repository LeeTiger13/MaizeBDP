import pandas as pd
import logging
from itertools import product
import argparse
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
    
def generate_offspring_genotypes(parent1, parent2):
    """
    Given two parental genotypes, generate all possible offspring genotypes.
    """
    # Extract unique alleles from both parents
    alleles = set(parent1 + parent2)

    # Generate all combinations of these alleles for the offspring
    offspring_genotypes = [''.join(pair) for pair in product(alleles, repeat=2)]

    # Filter out invalid combinations that don't respect dominance
    # assuming simple Mendelian inheritance without considering co-dominance complexities
    valid_genotypes = [gen for gen in offspring_genotypes if gen[0] <= gen[1]]

    return valid_genotypes


def check_offspring(parent1, parent2, offspring):
    """
    Judge whether the offspring genotype is formed by parents.
    """
    possible_genotypes = generate_offspring_genotypes(parent1, parent2)

    return int(offspring in possible_genotypes)


def parent_judge(input_parent, input_off):
    """
    Comparing different samples with genome window

    input_parent (str): Path to the input parent data
    input_off (str):  Path to the input offspring data
    """

    # log
    log_file = "parent_judge.log"
    setup_logger(log_file)
    logging.info("Starting parent judge...")

    # read file and pretreat
    parent_geno = pretreatment_out(input_parent, 0, 0, 0)
    off_geno = pretreatment_out(input_off, 0, 0, 0)

    print(parent_geno.head())
    print(off_geno.head())
    # rename
    parent_geno.columns = ['ID', 'Chr', 'Pos', 'Ref', 'parent1', 'parent2']
    col_list = ['ID', 'Chr', 'Pos', 'Ref']
    off_geno.rename(columns=dict(zip(off_geno.columns[:4], col_list)), inplace=True)

    # merge
    df = pd.merge(parent_geno, off_geno, on=['ID', 'Chr', 'Pos', 'Ref'], how='inner')

    # delete NA
    df = df.dropna()

    results = []
    # iterate over each offspring sample column
    for col in off_geno.columns[4:]:
        # convert the genotype to validity code
        df['validity_code'] = df.apply(lambda row: check_offspring(row['parent1'], row['parent2'], row[col]), axis=1)

        # create a new column for modified code with sliding window
        df['validity_code_md'] = df['validity_code']

        for chrom, group in df.groupby('Chr'):
            code_column = group['validity_code'].values
            code_md_column = group['validity_code_md'].values

            # initialize the start and end
            start = 0
            end = start + 20

            # When the end of the window does not exceed the length of the data, continue sliding
            while end <= len(code_column):
                # Get the data in the current window
                window = code_column[start:end]

                # Check whether there are at least 19 ones in the window
                if window.sum() >= 19:
                    # Replace all zeros in the window with ones
                    code_md_column[start:end] = 1

                # move the window
                start += 1
                end = start + 20

            # update the data
            df.loc[group.index, 'validity_code_md'] = code_md_column

        #  make sure the code_md is int
        df['validity_code_md'] = df['validity_code_md'].astype(int)

        # count the rate
        all_rate = sum(df['validity_code_md']) / df.shape[0]
        results.append({'Sample': col, 'Rate': all_rate})

    result_df = pd.DataFrame(results)
    result_df.to_csv('output.biparent.judge.txt', index=False, sep='\t')
    logging.info(f"Results saved.")


def main():
    parser = argparse.ArgumentParser(description="Comparing different samples with sliding window")
    parser.add_argument("-p", "--parent_file", type=str, help="Input the parent chip file")
    parser.add_argument("-i", "--offspring_file", type=str, help="Input the offspring chip file")
    args = parser.parse_args()

    # Call the function with the output file as the last argument
    parent_judge(args.parent_file, args.offspring_file)


if __name__ == "__main__":
    main()