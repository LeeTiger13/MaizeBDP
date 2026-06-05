library(rrBLUP)
library(cvTools)
library(optparse)

# Using rrblup do GS
# geno data : col=marker; row=Sample
# pheno data: col1=sample,col2=value
# def function
rrblup_gs <- function(geno_data, pheno_data, test_data = NULL){
  # read data
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
    for (j in 1:nfold) {
      ix <- ix + 1
      test_sample <- fold[[4]][fold[[5]]==j, 1] # get the test sample
      # training
      train_geno_temp <- train_geno[-test_sample, ]
      train_pheno_temp <- train_pheno[-test_sample, ]
      # testing
      test_geno_temp <- train_geno[test_sample, ]
      test_pheno_temp <- train_pheno[test_sample, ]
      # predicting
      model <- mixed.solve(train_pheno_temp, Z = train_geno_temp, K = NULL, SE =FALSE, return.Hinv = FALSE)
      GEBV <- as.matrix(test_geno_temp) %*% model$u # calculate the GEBV
      y_est_temp <- GEBV[, 1] + model$beta # predict the phenotype
      y_obs_temp <- test_pheno_temp # true phenotype
      # calculate the corr and mse,rmse
      output_matrix[ix, 1] <- ix
      output_matrix[ix, 2] <- cor(y_est_temp, y_obs_temp)
      output_matrix[ix, 3] <- mean((y_obs_temp - y_est_temp)^2)
      output_matrix[ix, 4] <- sqrt(mean((y_obs_temp - y_est_temp)^2))
  
  write.table(output_matrix, file = paste0('output_predict_accuracy.rrBLUP.txt'), 
              quote = FALSE, sep = '\t', row.names = FALSE)
    }
  }
  
  if (!is.null(test_data)){
    # read test_geno data
    test_geno <- read.table(test_data, header=TRUE, sep='\t',
                            stringsAsFactors = F, row.names = 1)
    
    # final model 
    final_model <- mixed.solve(train_pheno$Value, Z = train_geno, K = NULL, SE =FALSE, return.Hinv = FALSE)
    
    # predict phenotype 
    GEBV_test <- as.matrix(test_geno) %*% final_model$u
    y_pred_test <- GEBV_test[, 1] + final_model$beta
    
    # create a dataframe
    predicted_pheno <- data.frame(
      Sample = rownames(test_geno),
      Pre_value = y_pred_test
    )
    
    # output
    write.table(predicted_pheno, file = paste0('output_test_predicted_value.rrBLUP.txt'), 
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
rrblup_gs(opt$train_geno, opt$train_pheno, opt$test_geno)

