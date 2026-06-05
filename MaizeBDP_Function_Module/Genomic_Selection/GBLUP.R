# Load required libraries
library(rrBLUP)
library(cvTools)
library(optparse)

# Define the function for GBLUP Genomic Selection
Gblup_gs <- function(geno_data, pheno_data, test_data = NULL) {
  #--------------- 1. Read and Prepare Data ---------------#
  # Read genotype data for the training set
  train_geno <- read.table(geno_data,
    header = TRUE, sep = '\t',
    stringsAsFactors = FALSE, row.names = 1
  )
  # Read phenotype data for the training set
  train_pheno <- read.table(pheno_data,
    sep = '\t', header = TRUE,
    stringsAsFactors = FALSE, row.names = 1
  )
  colnames(train_pheno) <- 'Value'

  # Calculate the Genomic Relationship Matrix (Kinship Matrix) K
  # This is the core of GBLUP when using the kin.blup function.
  K <- A.mat(train_geno)

  #--------------- 2. Cross-Validation Setup ---------------#
  # Set parameters for cross-validation
  nfold <- 5 # 5-fold CV
  nrepeat <- 10 # 10 repetitions

  # Create an empty matrix to store the results from each fold/repeat
  output_matrix <- matrix(NA,
    nrow = nfold * nrepeat, ncol = 4,
    dimnames = list(NULL, c("Fold", "Correlation", "MSE", "RMSE"))
  )

  #--------------- 3. Run Cross-Validation ---------------#
  cat("Starting", nrepeat, "repeats of", nfold, "-fold cross-validation...\n")
  ix <- 0 # Initialize a counter for rows in the output matrix
  for (i in 1:nrepeat) {
    # The cvFolds are created INSIDE the repeat loop to ensure each repeat is different.
    set.seed(108 + i)
    folds <- cvFolds(nrow(train_pheno), K = nfold)

    for (j in 1:nfold) {
      ix <- ix + 1 # Increment counter

      # Get the indices for the test samples for the current fold
      test_indices <- folds$subsets[folds$which == j]

      # --- Prepare data for kin.blup ---
      # Copy the original phenotype vector
      y_na <- train_pheno$Value
      # Mask the phenotype of the test set with NA
      y_na[test_indices] <- NA

      # Create a data frame as required by kin.blup
      data_cv <- data.frame(y = y_na, gid = rownames(train_pheno))

      # --- Run the kin.blup model ---
      # The model predicts the NA values based on the kinship matrix
      model <- kin.blup(data = data_cv, geno = "gid", pheno = "y", K = K)

      # --- Evaluate the prediction ---
      # Get the predicted values for the test samples
      y_est_temp <- model$g[test_indices]
      # Get the true, observed values for the test samples
      y_obs_temp <- train_pheno$Value[test_indices]

      # Calculate and store performance metrics
      output_matrix[ix, 1] <- paste0("R", i, "_F", j)
      output_matrix[ix, 2] <- cor(y_est_temp, y_obs_temp)
      output_matrix[ix, 3] <- mean((y_obs_temp - y_est_temp)^2)
      output_matrix[ix, 4] <- sqrt(mean((y_obs_temp - y_est_temp)^2))
    }
  }

  # Write the cross-validation accuracy results to a file
  write.table(output_matrix,
    file = 'output_predict_accuracy.GBLUP.txt',
    quote = FALSE, sep = '\t', row.names = FALSE
  )
  cat("Cross-validation finished. Results saved to 'output_predict_accuracy.txt'.\n")

  #--------------- 4. Predict an Independent Test Set (Optional) ---------------#
  # This part executes only if a test genotype file is provided.
  if (!is.null(test_data)) {
    cat("Predicting phenotypes for the independent test set...\n")
    # Read the genotype data for the external test set
    test_geno <- read.table(test_data,
      header = TRUE, sep = '\t',
      stringsAsFactors = F, row.names = 1
    )

    # --- Final Model Training ---
    # To predict new individuals, we must estimate marker effects (u).
    # We use mixed.solve with the marker matrix (Z) for this.
    # This model is trained on the FULL training dataset.
    final_model <- mixed.solve(train_pheno$Value, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)

    # --- Predict Phenotypes for the Test Set ---
    # Calculate Genomic Estimated Breeding Values (GEBV)
    GEBV_test <- as.matrix(test_geno) %*% final_model$u
    # Add the overall mean (beta) to get the final predicted phenotype
    y_pred_test <- GEBV_test[, 1] + final_model$beta

    # Create a data frame for the output
    predicted_pheno <- data.frame(
      Sample = rownames(test_geno),
      Pre_value = y_pred_test
    )

    # Write the predicted phenotypes to a file
    write.table(predicted_pheno,
      file = 'output_test_predicted_value.GBLUP.txt',
      quote = FALSE, sep = '\t', row.names = FALSE
    )
    cat("Prediction on test set finished. Results saved to 'output_test_predicted_value.txt'.\n")
  }
}

#--------------- 5. Parse Command-Line Arguments and Run ---------------#
# Define the command-line options
option_list <- list(
  make_option(c("--train_geno"), type = "character", help = "File path for training genotype data"),
  make_option(c("--train_pheno"), type = "character", help = "File path for training phenotype data"),
  make_option(c("--test_geno"), type = "character", default = NULL, help = "Optional file path for testing genotype data")
)

# Parse the provided arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if mandatory arguments are provided
if (is.null(opt$train_geno) || is.null(opt$train_pheno)) {
  print_help(opt_parser)
  stop("Training genotype and phenotype files must be supplied.", call. = FALSE)
}

# Call the main function with the parsed arguments
Gblup_gs(opt$train_geno, opt$train_pheno, opt$test_geno)