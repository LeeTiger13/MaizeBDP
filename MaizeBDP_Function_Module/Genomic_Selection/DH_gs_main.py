#!/usr/bin/env python3
# Author: lihu0628@qq.com

import argparse
from DH_gs_base import pretreat_data, call_rrblup, call_Bayes, call_SVR, call_lgm, call_gblup

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Genomic Selection using different models")
    parser.add_argument("--genotype_train", required=True, help="Path to the genotype data for training")
    parser.add_argument("--phenotype_train", required=True, help="Path to the phenotype data for training")
    parser.add_argument("--genotype_test", help="Path to the genotype data for testing")
    parser.add_argument("--models", choices=["rrBLUP", "GBLUP", "BayesC", "SVR", "lightgbm"], required=True, help="Model for genomic selection")
    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_arguments()
    
    model_map = {
        "rrBLUP": {'func': call_rrblup, 'encoding': 'dominance'},
        "GBLUP": {'func': call_gblup, 'encoding': 'dominance'},
        "BayesC": {'func': call_Bayes, 'encoding': 'dominance'},
        "SVR": {'func': call_SVR, 'encoding': 'additive'},
        "lightgbm": {'func': call_lgm, 'encoding': 'additive'}
    }

    if args.models not in model_map:
        print(f"Error: Model '{args.models}' is not recognized.")
        return

    selected_model = model_map[args.models]
    
    print(f"--- Running Genomic Selection with {args.models} ---")
    
    # Preprocess data based on the model's required encoding
    if args.genotype_test:
        X_train, y_train, X_test = pretreat_data(
            args.genotype_train, 
            args.phenotype_train, 
            selected_model['encoding'], 
            args.genotype_test
        )
        # Call the selected model function
        selected_model['func'](X_train, y_train, X_test)
    else:
        X_train, y_train = pretreat_data(
            args.genotype_train, 
            args.phenotype_train, 
            selected_model['encoding']
        )
        # Call the selected model function
        selected_model['func'](X_train, y_train)

    print(f"--- {args.models} analysis complete. ---")

if __name__ == "__main__":
    main()