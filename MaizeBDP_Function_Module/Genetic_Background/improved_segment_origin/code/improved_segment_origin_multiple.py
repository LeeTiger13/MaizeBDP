import pandas as pd
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

def match_same_genotype(input_founder, input_improved, min_bin=3000000):
    """
    Find the segment of improved line which origined from improved line. One on One.

    input_founder (str): Path to the input founder line
    input_improved (str):  Path to the input improved line
    min_bin (int): The minimum bin for each region.
    """
    # log
    log_file = "improved_line_origin.log"    
    setup_logger(log_file)
    logging.info("Start to find the origin segment of improved line...")

    if not isinstance(min_bin, int):
        logging.error("Error: The parameter 'min_bin' must be of type int.")
        return None

    # rad file
    data_founder = pretreatment_out(input_founder, 0, 0, 0)
    if data_founder.shape[1] > 5:
        logging.error("Error: The founder line has more than 1 sample.")
        return None
    data_improved = pretreatment_out(input_improved, 0, 0, 0)
    if data_improved.shape[1] > 5:
        logging.error("Error: The improved line has more than 1 sample.")
        return None

    # combine data
    df = pd.merge(data_founder, data_improved, on=['ID', 'Chr', 'Pos', 'Ref'], how='inner')

    # rename
    df.columns = ['ID', 'Chr', 'Pos', 'Ref', 'founder_genotype', 'improved_genotype']

    # compare two sample, and set the same genotype SNP as 1, otherwise 0
    df['code'] = (df['founder_genotype'] == df['improved_genotype']).astype(int)

    # create a new column for modified code with sliding window
    df['code_md'] = df['code']

    # parameter
    window_size = 20
    step = 1
    threshold = 19

    for chrom, group in df.groupby('Chr'):
        code_column = group['code'].values
        code_md_column = group['code_md'].values

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
        df.loc[group.index, 'code_md'] = code_md_column

    #  make sure the code_md is int
    df['code_md'] = df['code_md'].astype(int)

    # set region with window size
    # create a new list for result
    result = []

    for chrom, group in df.groupby('Chr'):
        group['is_one'] = group['code_md'] == 1
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
    parser = argparse.ArgumentParser(description="Identify which region of the improved line were derived from the founder line")
    parser.add_argument("-f", "--founder_file", type=str, help="Path to the input founder line file")
    parser.add_argument("-i", "--improved_file", type=str, help="Path to the input improved line file")
    parser.add_argument("-m", "--min_bin", type=int, default=3000000, help="Minimum bin for region")
    args = parser.parse_args()

    # find the origin segment
    result_df = match_same_genotype(args.founder_file, args.improved_file, args.min_bin)
    result_df.to_csv('output.segment.origin.txt', index=False, sep='\t')
    logging.info("The origin segment of improved line has been saved to output.segment.origin.txt")

if __name__ == "__main__":
    main()
