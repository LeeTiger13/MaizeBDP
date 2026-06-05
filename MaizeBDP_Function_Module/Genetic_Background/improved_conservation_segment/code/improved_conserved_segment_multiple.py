import pandas as pd
import argparse
import math
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

def conserved_region_find(input_founder, input_improved, conserved_rate = 0.8, min_bin=3000000):
    """
    Identify conserved region of improved lines and judge it origin

    input_founder (str): Path to the input founder line
    input_improved (str):  Path to the input improved line
    conserved_rate (float): The rate for setting conserved region
    min_bin (int): The minimum bin for each region
    """
    # log
    log_file = "conserved_region_origin.log"
    setup_logger(log_file)
    logging.info("Start to find the origin segment of improved line...")
    
    if not isinstance(min_bin, int):
        logging.error("Error: The parameter 'min_bin' must be of type int.")
        return None

    if not isinstance(conserved_rate, float):
        logging.error("Error: The parameter 'conserved_rate' must be of type float.")
        return None

    # rad file
    data_founder = pretreatment_out(input_founder, 0, 0, 0)
    if data_founder.shape[1] > 5:
        logging.error("Error: The founder line has more than 1 sample.")
        return None
    data_improved = pretreatment_out(input_improved, 0, 0, 0)
    if data_improved.shape[1] <= 5:
        logging.error("Error: The improved line has less than 2 samples.")
        return None

    # combine data
    df = pd.merge(data_founder, data_improved, on=['ID', 'Chr', 'Pos', 'Ref'], how='inner')
    df.columns = (['ID', 'Chr', 'Pos', 'Ref', 'founder_genotype'] +
                  [f'improved_{i}' for i in range(1, data_improved.shape[1] - 3)])

    # Match each improved line and founder line and convert to a value
    improved_colname = df.columns[5:].tolist()
    for i in improved_colname:
        df[f'code_{i}'] = (df['founder_genotype'] == df[f'{i}']).astype(int)

    # The values of each improved line are converted by sliding window method
    code_colname_list = [f'code_{i}' for i in improved_colname]
    
    # parameter
    window_size = 20
    step = 1
    threshold = 19

    for i in code_colname_list:
        df[f'{i}_md'] = df[f'{i}']
        for chrom, group in df.groupby('Chr'):
            code_column = group[f'{i}'].values
            code_md_column = group[f'{i}_md'].values

            # initialize the start and end
            start = 0
            end = start + window_size

            # When the end of the window does not exceed the length of the data, continue sliding
            while end <= len(code_column):
                # Get the data in the current window
                window = code_column[start:end]

                # Check whether there are at least 19 ones in the window
                if window.sum() >= threshold:
                    # Replace all zeros in the window with ones
                    code_md_column[start:end] = 1

                # move the window
                start += step
                end = start + window_size

            # update the data
            df.loc[group.index, f'{i}_md'] = code_md_column

        df[f'{i}_md'] = df[f'{i}_md'].astype(int)

    # Calculate the sum of each material value
    code_md_colname_list = [f'code_{i}_md' for i in improved_colname]
    df['conserved_code'] = df[code_md_colname_list].sum(axis=1)

    # Further conversion of values according to conservative rate
    conserved_value = math.ceil(len(improved_colname) * conserved_rate)
    df['conserved_code_md'] = df['conserved_code'].apply(lambda x: 1 if x >= conserved_value else 0)

    # set region with window size
    # create a new list for result
    result = []

    for chrom, group in df.groupby('Chr'):
        group['is_one'] = group['conserved_code_md'] == 1
        group['segment_id'] = (group['is_one'] != group['is_one'].shift()).cumsum()
        segment_len = group[group['is_one']].groupby('segment_id').size()
        long_segments = segment_len[segment_len >= window_size].index
        result_group = (
            group[group['is_one']]
            .groupby(['segment_id'])
            .agg({'Pos': ['min', 'max']})
            .loc[long_segments]
            .reset_index()
        )

        result_group['Chr'] = chrom
        result_group.columns = ['Segment_ID', 'Start', 'End', 'Chr']
        result.append(result_group)

    result_df = pd.concat(result)
    result_df = result_df.drop(columns=['Segment_ID'])
    result_df = result_df[['Chr', 'Start', 'End']]
    result_df['region_size'] = result_df['End'].astype(int) - result_df['Start'].astype(int)
    result_df_md = result_df[result_df['region_size'] >= min_bin]

    return result_df_md

def main():
    parser = argparse.ArgumentParser(description="Identify conserved region of improved lines and judge it origin")
    parser.add_argument("-f", "--founder_file", type=str, help="Path to the input founder line file")
    parser.add_argument("-i", "--improved_file", type=str, help="Path to the input improved line file")
    parser.add_argument("-c", "--conserved_rate", type=float, default=0.8, help="conserved rate for filter region")
    parser.add_argument("-m", "--min_bin", type=int, default=3000000, help="Minimum bin for region")
    args = parser.parse_args()

    result_df = conserved_region_find(args.founder_file, args.improved_file, args.conserved_rate, args.min_bin)
    result_df.to_csv('output.conserved.segment.txt', index=False, sep='\t')
    logging.info("The conserved segment of improved line has been saved to output.conserved.segment.txt")

if __name__ == "__main__":
    main()
