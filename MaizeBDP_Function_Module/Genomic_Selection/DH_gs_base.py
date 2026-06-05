#!/usr/bin/env python3
# Author: lihu0628@qq.com

import pandas as pd
import numpy as np
import subprocess
import os
import glob
import lightgbm as lgb
import tempfile
from base_function_multiple import pretreatment_out
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV, RepeatedKFold
from sklearn.metrics import mean_squared_error

RSCRIPT_BIN = os.environ.get("RSCRIPT_BIN", "Rscript")

# --- Refactored Genotype Encoding ---

def encode_genotype(genotype_data, ref_alleles, encoding_scheme='additive'):
    """
    Encodes genotype strings into numeric values based on a specified scheme.
    This function is more robust and easier to understand.

    Args:
        genotype_data (np.array): 3D numpy array of shape (n_markers, n_samples, 2).
        ref_alleles (np.array): 1D array of reference alleles for each marker.
        encoding_scheme (str): 'additive' (0, 1, 2) or 'dominance' (-1, 0, 1).

    Returns:
        np.array: 2D array of encoded genotypes.
    """
    # Count the number of reference alleles for each genotype
    allele1 = genotype_data[:, :, 0]
    allele2 = genotype_data[:, :, 1]
    ref_allele_count = (allele1 == ref_alleles[:, None]) + (allele2 == ref_alleles[:, None])

    if encoding_scheme == 'additive':
        # Homozygous reference (2 ref alleles) -> 0
        # Heterozygous (1 ref allele) -> 1
        # Homozygous alternate (0 ref alleles) -> 2
        return 2 - ref_allele_count
    elif encoding_scheme == 'dominance':
        # Homozygous reference (2 ref alleles) -> -1
        # Heterozygous (1 ref allele) -> 0
        # Homozygous alternate (0 ref alleles) -> 1
        return ref_allele_count - 1
    else:
        raise ValueError("Unknown encoding_scheme. Use 'additive' or 'dominance'.")

# --- Refactored Pretreatment Function ---

def pretreat_data(train_geno_path, train_pheno_path, encoding_scheme, test_geno_path=None):
    """
    A single, unified function for data preprocessing.
    """
    # Genotype filtering
    train_genotype = pretreatment_out(train_geno_path, 0, 0, 0)
    train_markers = train_genotype.iloc[:, 0].values
    train_samples = train_genotype.columns[4:].values

    if test_geno_path:
        test_genotype = pretreatment_out(test_geno_path, 0, 0, 0)
        test_markers = test_genotype.iloc[:, 0].values
        test_samples = test_genotype.columns[4:].values
        
        common_markers = list(set(train_markers) & set(test_markers))
        
        train_geno_df = train_genotype[train_genotype.iloc[:, 0].isin(common_markers)]
        test_geno_df = test_genotype[test_genotype.iloc[:, 0].isin(common_markers)]
        
        train_geno_df = train_geno_df.set_index(train_geno_df.columns[0]).loc[common_markers].reset_index()
        test_geno_df = test_geno_df.set_index(test_geno_df.columns[0]).loc[common_markers].reset_index()
    else:
        train_geno_df = train_genotype
        test_geno_df = None
        test_samples = None

    # Helper function for encoding and transposing
    def process_genotypes(geno_df, samples, ref_alleles, scheme):
        geno_vals = geno_df.loc[:, samples].values
        geno_array_3d = np.array([np.array([list(gt) for gt in row]) for row in geno_vals])
        encoded_numeric = encode_genotype(geno_array_3d, ref_alleles, scheme)
        
        df = pd.DataFrame(encoded_numeric, columns=samples, index=geno_df.iloc[:, 0])
        df_t = df.T
        df_t.index.name = 'Sample'
        return df_t.reset_index()

    # Process training data
    ref_alleles = train_geno_df['Ref'].values
    X_train = process_genotypes(train_geno_df, train_samples, ref_alleles, encoding_scheme)
    
    # Process test data if it exists
    X_test = None
    if test_geno_df is not None:
        X_test = process_genotypes(test_geno_df, test_samples, ref_alleles, encoding_scheme)

    # Phenotype processing
    train_pheno_data = pd.read_csv(train_pheno_path, sep='\t', header=0)
    train_pheno_data.columns = ['Sample', 'Value']
    y_train = train_pheno_data[train_pheno_data['Sample'].isin(X_train['Sample'])]
    
    if test_geno_path:
        return X_train, y_train, X_test
    else:
        return X_train, y_train

# --- Refactored R Script Callers ---

def call_r_script(r_script_path, train_geno, train_pheno, test_geno=None):
    """
    A general function to call an R script with genotype and phenotype data,
    using secure temporary files.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_X_train_path = os.path.join(temp_dir, "X_train.csv")
        temp_y_train_path = os.path.join(temp_dir, "y_train.csv")
        
        train_geno.to_csv(temp_X_train_path, index=False, sep='\t')
        train_pheno.to_csv(temp_y_train_path, index=False, sep='\t')

        command = [RSCRIPT_BIN, r_script_path, "--train_geno", temp_X_train_path, "--train_pheno", temp_y_train_path]

        if test_geno is not None:
            temp_X_test_path = os.path.join(temp_dir, "X_test.csv")
            test_geno.to_csv(temp_X_test_path, index=False, sep='\t')
            command.extend(["--test_geno", temp_X_test_path])
        
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            print("R Script Standard Output:\n", result.stdout)
            if result.stderr:
                print("R Script Standard Error:\n", result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"--- Error running {r_script_path} ---\n", e.stderr, "\n---")

def call_rrblup(train_geno, train_pheno, test_geno=None):
    """Calls the rrblup.R script."""
    call_r_script("rrblup.R", train_geno, train_pheno, test_geno)

def call_gblup(train_geno, train_pheno, test_geno=None):
    """Calls the gblup.R script."""
    call_r_script("GBLUP.R", train_geno, train_pheno, test_geno)

def call_Bayes(train_geno, train_pheno, test_geno=None):
    """Calls the Bayes.R script and cleans up .dat files."""
    call_r_script("Bayes.R", train_geno, train_pheno, test_geno)
    for dat_file in glob.glob('*.dat'):
        if os.path.exists(dat_file):
            os.unlink(dat_file)

# --- Refactored Python Model Callers (SVR & LightGBM) ---

def run_cv_and_predict(model, param_grid, geno_data, pheno_data, test_data=None, output_prefix='model'):
    """
    A single, unified function to run GridSearchCV, report CV metrics, 
    and predict on a test set for scikit-learn compatible models.
    """
    # 1. Prepare data directly from memory
    merged_data = pd.merge(geno_data, pheno_data, on='Sample')
    train_pheno = merged_data['Value'].to_numpy()
    train_geno = merged_data.drop(columns=['Sample', 'Value']).to_numpy()
    
    # 2. Run GridSearchCV to find best parameters
    cv_K = RepeatedKFold(n_splits=5, n_repeats=10, random_state=42)
    grid_search = GridSearchCV(estimator=model, param_grid=param_grid, scoring='neg_mean_squared_error', cv=cv_K, verbose=1)
    grid_search.fit(train_geno, train_pheno)
    
    best_params = grid_search.best_params_
    best_model = grid_search.best_estimator_ # This model is already trained on the full dataset
    print(f"Best parameters for {output_prefix}: {best_params}")

    # 3. Manually loop through CV splits to generate detailed report
    results = []
    for i, (train_index, test_index) in enumerate(cv_K.split(train_geno)):
        X_train, X_test = train_geno[train_index], train_geno[test_index]
        y_train, y_test = train_pheno[train_index], train_pheno[test_index]
        
        # Instantiate and train model with the best found parameters
        eval_model = model.set_params(**best_params)
        eval_model.fit(X_train, y_train)
        y_pred = eval_model.predict(X_test)

        corr = np.corrcoef(y_test, y_pred)[0, 1] if np.std(y_pred) > 0 else 0
        mse = mean_squared_error(y_test, y_pred)
        results.append({'Repeat': i + 1, 'Correlation': corr, 'MSE': mse, 'RMSE': np.sqrt(mse)})

    # 4. Save results
    results_df = pd.DataFrame(results, columns=['Repeat', 'Correlation', 'MSE', 'RMSE'])
    results_df.to_csv(f'output_accuracy.{output_prefix}.txt', index=False, sep='\t')

    # 5. Predict on test set if provided
    if test_data is not None:
        test_samples = test_data['Sample']
        test_geno_array = test_data.drop(columns=['Sample']).to_numpy()
        test_pred = best_model.predict(test_geno_array)
        test_pred_df = pd.DataFrame({'Sample': test_samples, 'Prediction': test_pred})
        test_pred_df.to_csv(f'output_predicted_value.{output_prefix}.txt', index=False, sep='\t')

def call_SVR(geno_data, pheno_data, test_data=None):
    """Wrapper function for SVR."""
    svr = SVR()
    param_grid = {'C': [0.1, 1, 10], 'gamma': [0.001, 0.01, 0.1, 1], 'kernel': ['rbf']}
    run_cv_and_predict(svr, param_grid, geno_data, pheno_data, test_data, output_prefix='svr')

def call_lgm(geno_data, pheno_data, test_data=None):
    """Wrapper function for LightGBM."""
    lgbm = lgb.LGBMRegressor(random_state=42)
    param_grid = {
        'num_leaves': [15, 31], 'max_depth': [3, 5], 'min_child_samples': [20, 50],
        'learning_rate': [0.05, 0.1], 'n_estimators': [100, 200]
    }
    run_cv_and_predict(lgbm, param_grid, geno_data, pheno_data, test_data, output_prefix='lgbm')
