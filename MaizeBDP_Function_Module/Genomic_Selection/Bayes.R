library(BGLR)
library(cvTools)
library(optparse)

# using BayesC do GS
# def function
bayesC_gs <- function(geno_data, pheno_data, test_data = NULL){
  # load file
  train_geno <- read.table(geno_data, header=TRUE, sep='\t',
                           stringsAsFactors = FALSE, row.names = 1)
  train_pheno <- read.table(pheno_data, sep='\t', header = TRUE, 
                            stringsAsFactors = FALSE, row.names = 1)
  colnames(train_pheno) <- 'Value'
  
  # set the relative parameter
  nfold <- 5
  nrepeat <- 10
  
  # set the output file 
  output_matrix <- matrix(NA, nrow = nfold*nrepeat, ncol = 4, 
                          dimnames = list(c(1:(nfold*nrepeat)), c("Repeat","Correlation", "MSE", "RMSE")))
  
  # begin the prediction
  ix <- 0 # set the row of output
  for (i in 1:nrepeat){
    # cycle for repeat
    set.seed(108+i)
  fold <- cvFolds(nrow(train_pheno), K = nfold)
    for (j in 1:nfold){
      ix <- ix + 1
      test_sample <- fold[[4]][fold[[5]]==j, 1] # get the test sample
      # training
      train_geno_temp <- train_geno[-test_sample, ]
      train_pheno_temp <- train_pheno[-test_sample, ]
      # testing
      test_geno_temp <- train_geno[test_sample, ]
      test_pheno_temp <- train_pheno[test_sample, ]
      # predicting
      model <- BGLR(y=train_pheno_temp, ETA = list(list(X=train_geno_temp, model="BayesC")),
                    nIter=6000, burnIn = 1000, saveAt = 'bc_')
      GEBV <- as.matrix(test_geno_temp) %*% model$ETA[[1]]$b
      y_est_temp <- GEBV[, 1] + model$mu
      y_obs_temp <- test_pheno_temp
      # calculate
      output_matrix[ix, 1] <- ix
      output_matrix[ix, 2] <- cor(y_est_temp, y_obs_temp)
      output_matrix[ix, 3] <- mean((y_obs_temp - y_est_temp)^2)
      output_matrix[ix, 4] <- sqrt(mean((y_obs_temp - y_est_temp)^2))
    }
  }
  write.table(output_matrix, file = paste0('output_predict_accuracy.BayesC.txt'), 
              quote = FALSE, sep = '\t', row.names = FALSE)
  
  # test
  if (!is.null(test_data)){
    # read test_geno data
    test_geno <- read.table(test_data, header=TRUE, sep='\t',
                            stringsAsFactors = F, row.names = 1)
    # final model
    final_model <- BGLR(y=train_pheno$Value, ETA = list(list(X=train_geno, model="BayesC")),
                        nIter=6000, burnIn = 1000, saveAt = 'bc_')
    GEBV_test <- as.matrix(test_geno) %*% final_model$ETA[[1]]$b
    y_pred_test <- GEBV_test[, 1] + final_model$mu
    # create a dataframe
    predicted_pheno <- data.frame(
      Sample = rownames(test_geno),
      Pre_value = y_pred_test
    )
    
    # output
    write.table(predicted_pheno, file = paste0('output_test_predicted_value.BayesC.txt'), 
                quote = FALSE, sep = '\t', row.names = FALSE)
  }
}

# set parameters
option_list <- list(
  make_option(c("--train_geno"), type = "character", help = "Training genotype"),
  make_option(c("--train_pheno"), type = "character", help = "Training phenotype"),
  make_option(c("--test_geno"), type = "character", help = "Testing genotype")
)

# parse the parameters
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# call
bayesC_gs(opt$train_geno, opt$train_pheno, opt$test_geno)

