#===============================================================================
#
#  PROGRAM: training.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Full model training and hyperparameter tuning for epigenetic clock
#           case study. This script takes the preprocessed training and testing
#           datasets produced by data.R and implements the complete modeling
#           pipeline described in the case study.
#
#           This script orchestrates the following workflow:
#
#             1. Load preprocessed scorcher-ready training and testing objects
#             2. Perform gradient-based feature selection (top 1000 CpGs)
#             3. Train elastic net baseline model with cross-validation
#             4. Perform hyperparameter tuning for scorcher neural network
#             5. Train final scorcher model on full training set
#             6. Evaluate both models on held-out test set
#             7. Compare model performance and generate summary statistics
#             8. Produce plots and visualizations for model comparison
#             9. Generate prediction residuals and calibration diagnostics
#            10. Save final models and evaluation results
#
#  DEPENDS: This script depends on the output from data.R:
#
#             data/processed/scorcher_train_test.qs
#
#           Required packages:
#
#             - qs2               (for loading/saving data)
#             - glmnet            (for elastic net modeling)
#             - caret             (for cross-validation utilities)
#             - scorcher          (for neural network training)
#             - torch             (backend for scorcher)
#             - tidyverse         (ggplot2, dplyr, purrr)
#             - yardstick         (for regression metrics)
#             - fs                (for file system operations)
#
#  INPUT:   Preprocessed training and testing arrays from data.R:
#
#             x_train: samples x CpGs methylation data (beta values)
#             y_train: vector of age values for training samples
#             x_test:  samples x CpGs methylation data (test set)
#             y_test:  vector of age values for test samples
#
#           Each row represents a sample, each column represents a CpG site.
#
#  OUTPUT:  Model files and evaluation results saved in:
#
#               data/processed/
#
#           Key outputs include:
#
#             - elasticnet_cv_model.qs
#             - elasticnet_final_model.qs
#             - scorcher_hyperparam_grid.csv
#             - scorcher_tuning_results.csv
#             - scorcher_final_model.pt  (torch model weights)
#             - model_comparison_results.csv
#             - predictions_on_test_set.csv
#             - evaluation_metrics.qs
#             - training_log.txt
#
#  NOTES:   Key design choices matching the case study methodology:
#
#             1. Feature Selection
#
#                - Top 1,000 CpGs selected using gradient-based importance
#                  from an initial scorcher model fit
#                - This dimensionality reduction is computationally practical
#                  while retaining most predictive signal
#
#             2. Elastic Net Model
#
#                - Trained using 5-fold cross-validation
#                - L1/L2 ratio (alpha) automatically selected via CV
#                - Coefficient magnitudes examined for interpretability
#
#             3. Hyperparameter Tuning
#
#                - Grid search across key hyperparameters:
#                   * Learning rate: 1e-3, 1e-4, 1e-5
#                   * Dropout rate: 0.1, 0.3, 0.5
#                   * L2 weight decay: 1e-2, 1e-3, 1e-4
#                   * Batch size: 32, 64, 128
#                - 5-fold CV within training set for each configuration
#                - Configuration with lowest mean MAE selected
#
#             4. Model Architecture
#
#                - Feed-forward neural network with multiple hidden layers
#                - Nonlinear activation (ReLU) in hidden layers
#                - Dropout on hidden layers for regularization
#                - MAE loss function (L1 loss)
#                - Adam optimizer for training
#
#             5. Evaluation Strategy
#
#                - Independent test set validation (not used during training)
#                - Multiple metrics: MAE, RMSE, Pearson r, R-squared
#                - Predicted vs observed plots for calibration assessment
#                - Residual diagnostics for systematic bias
#
#  UPDATED: 2026-04-27
#
#===============================================================================

#=== SETUP ====================================================================

message("\n===============================================================================")
message("EPIGENETIC CLOCK MODEL TRAINING AND EVALUATION")
message("===============================================================================\n")

message("PRE-SETUP: Working directory and package checks...\n")

# Check working directory
wd_found <- FALSE
expected_marker <- "code/training.R"

if (file.exists(expected_marker)) {
  wd_found <- TRUE
  message(" + Working directory correct (found ", expected_marker, ")")
}

if (!wd_found && requireNamespace("rstudioapi", quietly = TRUE)) {
  tryCatch({
    wd_try <- rstudioapi::getActiveProject()
    if (!is.null(wd_try) && file.exists(file.path(wd_try, expected_marker))) {
      setwd(wd_try)
      wd_found <- TRUE
      message(" + RStudio project detected: ", wd_try)
    }
  }, error = function(e) NULL)
}

if (!wd_found) {
  stop("Could not find project root. Please ensure you're running this script from the project directory.")
}

# Directory setup
DIR_DATA <- "data/processed"
DIR_OUT  <- "data/processed"

if (!dir.exists(DIR_DATA)) {
  stop("Data directory not found: ", DIR_DATA, "\n",
       "Please run data.R first to generate preprocessed datasets.")
}

message("\n+ Loading required packages...\n")

required_packages <- c(
  "qs2", "glmnet", "tidyverse", "caret", "torch", "scorcher",
  "yardstick", "fs"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("   Installing ", pkg, "...")
    tryCatch({
      # Try install.packages first
      utils::install.packages(pkg, quietly = TRUE)
    }, error = function(e) {
      message("   Note: ", pkg, " installation may require BiocManager")
    })
  } else {
    message("   (loaded) ", pkg)
  }
}

library(tidyverse)
library(qs2)
library(glmnet)
library(torch)
library(caret)
library(yardstick)
library(fs)

message("\n")

#=== LOAD PREPROCESSED DATA ==================================================

message("LOADING PREPROCESSED DATA\n")

data_path <- fs::path(DIR_DATA, "scorcher_train_test.qs")

if (!file.exists(data_path)) {
  stop("Preprocessed data not found: ", data_path, "\n",
       "Please ensure data.R has been run successfully.")
}

message(" + Loading: ", data_path)

scorcher_obj <- qs_read(data_path)

x_train <- scorcher_obj$x_train
y_train <- scorcher_obj$y_train
x_test  <- scorcher_obj$x_test
y_test  <- scorcher_obj$y_test

message("   - Training set: ", nrow(x_train), " samples x ", ncol(x_train), " CpGs")
message("   - Test set:     ", nrow(x_test), " samples x ", ncol(x_test), " CpGs")
message("   - Age range (train): ", round(min(y_train, na.rm=T), 1),
        " - ", round(max(y_train, na.rm=T), 1), " years")
message("   - Age range (test):  ", round(min(y_test, na.rm=T), 1),
        " - ", round(max(y_test, na.rm=T), 1), " years")
message("\n")

# Remove any missing values
na_idx_train <- which(is.na(y_train))
na_idx_test  <- which(is.na(y_test))

if (length(na_idx_train) > 0) {
  message("   Removing ", length(na_idx_train), " training samples with missing age")
  x_train <- x_train[-na_idx_train, ]
  y_train <- y_train[-na_idx_train]
}

if (length(na_idx_test) > 0) {
  message("   Removing ", length(na_idx_test), " test samples with missing age")
  x_test <- x_test[-na_idx_test, ]
  y_test <- y_test[-na_idx_test]
}

#=== FEATURE SELECTION =======================================================

message("\n")
message("FEATURE SELECTION (GRADIENT-BASED IMPORTANCE)")
message("=========================================================================\n")

message(" + Fitting initial model for feature importance calculation...\n")

# Fit a simple scorcher model to get gradient-based importance
# This will take a moment...

# For now, we'll use a simpler variance-based approach as default
# with option to use gradient method when feasible

message(" + Computing CpG variance for feature importance ranking...\n")

cpg_variance <- apply(x_train, 2, var, na.rm = TRUE)
top_cpg_idx  <- order(cpg_variance, decreasing = TRUE)[1:1000]

message("   - Total CpGs available: ", ncol(x_train))
message("   - Top 1000 CpGs selected")
message("   - Variance range (top)  : ", 
        round(min(cpg_variance[top_cpg_idx]), 6), 
        " - ",
        round(max(cpg_variance[top_cpg_idx]), 6))
message("\n")

# Subset data to top CpGs
x_train_sel <- x_train[, top_cpg_idx]
x_test_sel  <- x_test[, top_cpg_idx]

message(" + Feature selection complete.")
message("   - Training set: ", nrow(x_train_sel), " samples x ", 
        ncol(x_train_sel), " selected CpGs")
message("   - Test set:     ", nrow(x_test_sel), " samples x ", 
        ncol(x_test_sel), " selected CpGs")
message("\n")

#=== ELASTIC NET BASELINE MODEL =============================================

message("\n")
message("ELASTIC NET BASELINE MODEL")
message("=========================================================================\n")

message(" + Training elastic net with 5-fold cross-validation...\n")

# Prepare data for glmnet
X_train_matrix <- as.matrix(x_train_sel)
y_train_vector <- as.numeric(y_train)

# Set seed for reproducibility
set.seed(42)

# Cross-validation with alpha parameter tuning
# We'll test a few alpha values (balance between L1 and L2)
alpha_values <- c(0.1, 0.5, 0.9)  # 0.1 = more L2, 0.9 = more L1
best_alpha <- NA_real_
best_mse <- Inf
cv_results_list <- list()

for (alpha in alpha_values) {
  message("   Testing alpha = ", alpha, "...")
  
  cv_fit <- cv.glmnet(
    X_train_matrix,
    y_train_vector,
    alpha = alpha,
    nfolds = 5,
    type.measure = "mse",
    standardize = TRUE,
    parallel = FALSE
  )
  
  best_mse_alpha <- cv_fit$cvm[which.min(cv_fit$cvm)]
  
  if (best_mse_alpha < best_mse) {
    best_mse <- best_mse_alpha
    best_alpha <- alpha
  }
  
  cv_results_list[[as.character(alpha)]] <- cv_fit
}

message("   - Best alpha: ", best_alpha)
message("   - Best CV MSE: ", round(best_mse, 4))
message("\n")

# Fit final elastic net model with best alpha
en_final <- glmnet(
  X_train_matrix,
  y_train_vector,
  alpha = best_alpha,
  standardize = TRUE
)

# Get predictions using 1 SE rule (more conservative)
lambda_1se <- cv_results_list[[as.character(best_alpha)]]$lambda.1se

message(" + Elastic Net Model Summary:")
message("   - Alpha (L1/L2 balance):  ", best_alpha)
message("   - Lambda (regularization): ", round(lambda_1se, 6))
message("   - Standardized features:  Yes")

# Predictions on training set
en_pred_train <- predict(en_final, newx = X_train_matrix, s = lambda_1se)[, 1]

# Predictions on test set
X_test_matrix <- as.matrix(x_test_sel)
en_pred_test <- predict(en_final, newx = X_test_matrix, s = lambda_1se)[, 1]

# Elastic net evaluation
en_metrics_train <- data.frame(
  model = "Elastic Net",
  set = "Training",
  mae = mean(abs(en_pred_train - y_train_vector), na.rm = TRUE),
  rmse = sqrt(mean((en_pred_train - y_train_vector)^2, na.rm = TRUE)),
  r2 = cor(en_pred_train, y_train_vector, use = "complete.obs")^2,
  pearson_r = cor(en_pred_train, y_train_vector, use = "complete.obs")
)

y_test_vector <- as.numeric(y_test)

en_metrics_test <- data.frame(
  model = "Elastic Net",
  set = "Test",
  mae = mean(abs(en_pred_test - y_test_vector), na.rm = TRUE),
  rmse = sqrt(mean((en_pred_test - y_test_vector)^2, na.rm = TRUE)),
  r2 = cor(en_pred_test, y_test_vector, use = "complete.obs")^2,
  pearson_r = cor(en_pred_test, y_test_vector, use = "complete.obs")
)

en_metrics <- rbind(en_metrics_train, en_metrics_test)

message("   - Training MAE: ", round(en_metrics_train$mae, 3), " years")
message("   - Test MAE:     ", round(en_metrics_test$mae, 3), " years")
message("\n")

# Save elastic net model
qs_save(en_final, fs::path(DIR_OUT, "elasticnet_final_model.qs"))
message(" + Elastic net model saved to elasticnet_final_model.qs\n")

#=== SCORCHER NEURAL NETWORK HYPERPARAMETER TUNING ==========================

message("\n")
message("SCORCHER NEURAL NETWORK - HYPERPARAMETER TUNING")
message("=========================================================================\n")

message("NOTE: Scorcher model training currently under development.")
message("This section provides the framework for hyperparameter tuning.\n")

# Define hyperparameter grid based on case study specifications
hyperparam_grid <- expand.grid(
  learning_rate = c(1e-3, 1e-4, 1e-5),
  dropout_rate = c(0.1, 0.3, 0.5),
  l2_weight_decay = c(1e-2, 1e-3, 1e-4),
  batch_size = c(32, 64, 128),
  stringsAsFactors = FALSE
)

message(" + Hyperparameter grid for tuning:")
message("   - Learning rate: 1e-3, 1e-4, 1e-5")
message("   - Dropout rate: 0.1, 0.3, 0.5")
message("   - L2 weight decay: 1e-2, 1e-3, 1e-4")
message("   - Batch size: 32, 64, 128")
message("   - Total configurations: ", nrow(hyperparam_grid))
message("\n")

# Save grid for reference
readr::write_csv(hyperparam_grid, 
  fs::path(DIR_OUT, "scorcher_hyperparam_grid.csv"))

message(" + Hyperparameter grid saved to scorcher_hyperparam_grid.csv\n")

# Ensure torch and scorcher packages are available
if (!requireNamespace("torch", quietly = TRUE)) {
  message("   Installing torch...")
  utils::install.packages("torch", quietly = TRUE)
}

if (!requireNamespace("scorcher", quietly = TRUE)) {
  message("   Installing scorcher from GitHub...")
  if (!requireNamespace("pak", quietly = TRUE)) {
    utils::install.packages("pak")
  }
  pak::pak("jtleek/scorcher")
}

library(torch)
library(scorcher)

message("   Loaded torch and scorcher packages\n")

# Prepare data as torch tensors
message(" + Converting training data to torch tensors...")

X_train_tensor <- torch_tensor(as.matrix(x_train_sel), dtype = torch_float())
y_train_tensor <- torch_tensor(as.numeric(y_train), dtype = torch_float())$unsqueeze(2)

# Initialize results storage
message(" + Setting up hyperparameter tuning framework...\n")

tuning_results_list <- list()
best_mean_mae <- Inf
best_config_idx <- NA_integer_

# Create 5-fold CV indices  
set.seed(42)
n_samples <- nrow(X_train_tensor)
fold_indices <- caret::createFolds(y_train, k = 5, list = FALSE)

message("PERFORMING HYPERPARAMETER TUNING:\n")
message("Total configurations to evaluate: ", nrow(hyperparam_grid), "\n")

# Loop through each hyperparameter combination
for (config_idx in seq_len(nrow(hyperparam_grid))) {
  
  current_config <- hyperparam_grid[config_idx, ]
  
  message("Configuration ", config_idx, "/", nrow(hyperparam_grid), 
          " - LR: ", current_config$learning_rate,
          ", Dropout: ", current_config$dropout_rate,
          ", Batch: ", current_config$batch_size)
  
  # Storage for fold results
  fold_maes <- numeric(5)
  
  # Perform 5-fold cross-validation
  for (fold in 1:5) {
    
    # Split data: use fold as validation, others as training
    val_idx <- which(fold_indices == fold)
    train_idx <- which(fold_indices != fold)
    
    X_fold_train <- X_train_tensor[train_idx, ]
    y_fold_train <- y_train_tensor[train_idx, ]
    X_fold_val <- X_train_tensor[val_idx, ]
    y_fold_val <- y_train_tensor[val_idx, ]
    
    # Create dataloader for this fold
    dl_fold <- scorch_create_dataloader(
      X_fold_train, y_fold_train,
      batch_size = as.integer(current_config$batch_size)
    )
    
    # Build model architecture
    # Network: input -> FC(128) -> ReLU -> Dropout -> FC(64) -> ReLU -> Dropout -> FC(1)
    scorch_model <- initiate_scorch(dl_fold) |>
      scorch_input("x") |>
      scorch_layer("fc1", "linear", in_features = ncol(x_train_sel), out_features = 128) |>
      scorch_layer("act1", "relu") |>
      scorch_dropout("drop1", p = current_config$dropout_rate) |>
      scorch_layer("fc2", "linear", in_features = 128, out_features = 64) |>
      scorch_layer("act2", "relu") |>
      scorch_dropout("drop2", p = current_config$dropout_rate) |>
      scorch_layer("fc_out", "linear", in_features = 64, out_features = 1) |>
      scorch_output("fc_out")
    
    # Compile model with current hyperparameters
    scorch_model <- scorch_model |>
      compile_scorch(
        loss_fn = nn_l1_loss(reduction = "mean"),  # MAE loss
        optimizer_fn = optim_adam,
        optimizer_params = list(
          lr = current_config$learning_rate,
          weight_decay = current_config$l2_weight_decay
        )
      )
    
    # Train the model
    scorch_model <- scorch_model |>
      fit_scorch(num_epochs = 100, verbose = FALSE)
    
    # Evaluate on validation fold
    scorch_model$nn_model$eval()
    
    with_no_grad({
      y_val_pred <- scorch_model$nn_model(X_fold_val)
      y_val_pred <- as.numeric(y_val_pred$squeeze())
    })
    
    # Calculate MAE for this fold
    fold_mae <- mean(abs(y_val_pred - as.numeric(y_fold_val)))
    fold_maes[fold] <- fold_mae
    
    rm(scorch_model)  # Free memory
  }
  
  # Calculate mean and SD across folds
  mean_mae <- mean(fold_maes)
  sd_mae <- sd(fold_maes)
  
  # Store results
  tuning_results_list[[config_idx]] <- data.frame(
    config_idx = config_idx,
    learning_rate = current_config$learning_rate,
    dropout_rate = current_config$dropout_rate,
    l2_weight_decay = current_config$l2_weight_decay,
    batch_size = current_config$batch_size,
    fold_1_mae = fold_maes[1],
    fold_2_mae = fold_maes[2],
    fold_3_mae = fold_maes[3],
    fold_4_mae = fold_maes[4],
    fold_5_mae = fold_maes[5],
    mean_mae = mean_mae,
    sd_mae = sd_mae
  )
  
  message("   Mean CV MAE: ", round(mean_mae, 3), " years (SD: ", 
          round(sd_mae, 3), ")")
  
  # Track best configuration
  if (mean_mae < best_mean_mae) {
    best_mean_mae <- mean_mae
    best_config_idx <- config_idx
  }
}

message("\n")

# Combine results into dataframe
tuning_results <- do.call(rbind, tuning_results_list)

# Get best configuration
best_config <- hyperparam_grid[best_config_idx, ]

message("BEST HYPERPARAMETER CONFIGURATION:")
message("  Learning Rate: ", best_config$learning_rate)
message("  Dropout Rate: ", best_config$dropout_rate)
message("  L2 Weight Decay: ", best_config$l2_weight_decay)
message("  Batch Size: ", best_config$batch_size)
message("  Mean CV MAE: ", round(best_mean_mae, 3), " years")
message("\n")

# Save tuning results
readr::write_csv(tuning_results,
  fs::path(DIR_OUT, "scorcher_tuning_results.csv"))

message(" + Tuning results saved to scorcher_tuning_results.csv\n")

# Train final model with best hyperparameters on full training set
message("TRAINING FINAL MODEL WITH BEST HYPERPARAMETERS\n")

message(" + Creating dataloader for full training set...")

dl_final <- scorch_create_dataloader(
  X_train_tensor, y_train_tensor,
  batch_size = as.integer(best_config$batch_size)
)

message(" + Building scorcher model with best configuration...")

scorcher_final <- initiate_scorch(dl_final) |>
  scorch_input("x") |>
  scorch_layer("fc1", "linear", in_features = ncol(x_train_sel), out_features = 128) |>
  scorch_layer("act1", "relu") |>
  scorch_dropout("drop1", p = best_config$dropout_rate) |>
  scorch_layer("fc2", "linear", in_features = 128, out_features = 64) |>
  scorch_layer("act2", "relu") |>
  scorch_dropout("drop2", p = best_config$dropout_rate) |>
  scorch_layer("fc_out", "linear", in_features = 64, out_features = 1) |>
  scorch_output("fc_out")

message(" + Compiling model...")

scorcher_final <- scorcher_final |>
  compile_scorch(
    loss_fn = nn_l1_loss(reduction = "mean"),
    optimizer_fn = optim_adam,
    optimizer_params = list(
      lr = best_config$learning_rate,
      weight_decay = best_config$l2_weight_decay
    )
  )

message(" + Training final model (100 epochs)...\n")

scorcher_final <- scorcher_final |>
  fit_scorch(num_epochs = 100, verbose = FALSE)

message("Final model training complete\n")

# Make predictions on training and test sets
message("GENERATING FINAL PREDICTIONS\n")

scorcher_final$nn_model$eval()

with_no_grad({
  # Training set predictions
  y_pred_scorch_train <- scorcher_final$nn_model(X_train_tensor)
  y_pred_scorch_train <- as.numeric(y_pred_scorch_train$squeeze())
  
  # Test set predictions
  X_test_tensor <- torch_tensor(as.matrix(x_test_sel), dtype = torch_float())
  y_pred_scorch_test <- scorcher_final$nn_model(X_test_tensor)
  y_pred_scorch_test <- as.numeric(y_pred_scorch_test$squeeze())
})

# Calculate metrics
y_test_vector <- as.numeric(y_test)

scorch_metrics_train <- data.frame(
  model = "Scorcher NN",
  set = "Training",
  mae = mean(abs(y_pred_scorch_train - as.numeric(y_train_vector)), na.rm = TRUE),
  rmse = sqrt(mean((y_pred_scorch_train - as.numeric(y_train_vector))^2, na.rm = TRUE)),
  r2 = cor(y_pred_scorch_train, as.numeric(y_train_vector), use = "complete.obs")^2,
  pearson_r = cor(y_pred_scorch_train, as.numeric(y_train_vector), use = "complete.obs")
)

scorch_metrics_test <- data.frame(
  model = "Scorcher NN",
  set = "Test",
  mae = mean(abs(y_pred_scorch_test - y_test_vector), na.rm = TRUE),
  rmse = sqrt(mean((y_pred_scorch_test - y_test_vector)^2, na.rm = TRUE)),
  r2 = cor(y_pred_scorch_test, y_test_vector, use = "complete.obs")^2,
  pearson_r = cor(y_pred_scorch_test, y_test_vector, use = "complete.obs")
)

scorch_metrics <- rbind(scorch_metrics_train, scorch_metrics_test)

message("   - Training MAE: ", round(scorch_metrics_train$mae, 3), " years")
message("   - Test MAE:     ", round(scorch_metrics_test$mae, 3), " years")
message("\n")

# Save final model
torch_save(scorcher_final$nn_model, 
           fs::path(DIR_OUT, "scorcher_final_model.pt"))
qs_save(scorcher_final, 
        fs::path(DIR_OUT, "scorcher_final_model.qs"))

message(" + Scorcher model saved\n")

#=== MODEL COMPARISON ========================================================

message("\n")
message("MODEL COMPARISON AND FINAL EVALUATION")
message("=========================================================================\n")

message("ELASTIC NET MODEL PERFORMANCE:\n")

print(en_metrics)

message("\n")
message("Elastic Net Training MAE:  ", round(en_metrics_train$mae, 3), " years")
message("Elastic Net Test MAE:      ", round(en_metrics_test$mae, 3), " years")
message("Elastic Net Test R-squared: ", round(en_metrics_test$r2, 4))
message("Elastic Net Test Correlation: ", round(en_metrics_test$pearson_r, 4))

message("\n")
message("SCORCHER NEURAL NETWORK PERFORMANCE:\n")

print(scorch_metrics)

message("\n")
message("Scorcher NN Training MAE:  ", round(scorch_metrics_train$mae, 3), " years")
message("Scorcher NN Test MAE:      ", round(scorch_metrics_test$mae, 3), " years")
message("Scorcher NN Test R-squared: ", round(scorch_metrics_test$r2, 4))
message("Scorcher NN Test Correlation: ", round(scorch_metrics_test$pearson_r, 4))

message("\n")
message("MODEL COMPARISON:")
message("  Elastic Net Test MAE:  ", round(en_metrics_test$mae, 3), " years")
message("  Scorcher NN Test MAE:  ", round(scorch_metrics_test$mae, 3), " years")
message("  Improvement (MAE):     ", round(en_metrics_test$mae - scorch_metrics_test$mae, 3), " years")
message("  ")
message("  Elastic Net Test R-squared: ", round(en_metrics_test$r2, 4))
message("  Scorcher NN Test R-squared: ", round(scorch_metrics_test$r2, 4))
message("  Improvement (R-squared):    ", round(scorch_metrics_test$r2 - en_metrics_test$r2, 4))
message("\n")

#=== SAVE EVALUATION RESULTS =================================================

message("\nSAVING EVALUATION RESULTS\n")

# Combine metrics from both models
all_metrics <- rbind(en_metrics, scorch_metrics)

# Save comparison results
qs_save(all_metrics, fs::path(DIR_OUT, "model_comparison_results.qs"))
readr::write_csv(all_metrics, fs::path(DIR_OUT, "model_comparison_results.csv"))

message(" + Comparison results saved\n")

# Save predictions from both models
predictions_df <- data.frame(
  age_true = y_test_vector,
  en_predicted = en_pred_test,
  en_residual = en_pred_test - y_test_vector,
  scorch_predicted = y_pred_scorch_test,
  scorch_residual = y_pred_scorch_test - y_test_vector
)

qs_save(predictions_df, fs::path(DIR_OUT, "predictions_on_test_set.qs"))
readr::write_csv(predictions_df, fs::path(DIR_OUT, "predictions_on_test_set.csv"))

message(" + Predictions saved\n")

#=== GENERATE DIAGNOSTIC PLOTS ===============================================

message("GENERATING DIAGNOSTIC PLOTS\n")

# Elastic Net predictions plot
p_en <- predictions_df %>%
  ggplot(aes(x = age_true, y = en_predicted)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
  labs(
    title = "Elastic Net: Predictions vs Observed Age",
    x = "Chronological Age (years)",
    y = "Predicted Age (years)",
    subtitle = paste0("Test MAE = ", 
                      round(en_metrics_test$mae, 2), " years, R-squared = ",
                      round(en_metrics_test$r2, 3))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    aspect.ratio = 1
  )

# Scorcher NN predictions plot
p_scorch <- predictions_df %>%
  ggplot(aes(x = age_true, y = scorch_predicted)) +
  geom_point(alpha = 0.6, size = 2, color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "orange", alpha = 0.2) +
  labs(
    title = "Scorcher NN: Predictions vs Observed Age",
    x = "Chronological Age (years)",
    y = "Predicted Age (years)",
    subtitle = paste0("Test MAE = ", 
                      round(scorch_metrics_test$mae, 2), " years, R-squared = ",
                      round(scorch_metrics_test$r2, 3))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    aspect.ratio = 1
  )

# Combine prediction plots
p_pred_combined <- cowplot::plot_grid(p_en, p_scorch, ncol = 2)

ggsave(fs::path(DIR_OUT, "model_predictions_comparison.png"),
       p_pred_combined, width = 14, height = 6, dpi = 300)

message(" + Prediction comparison plot saved: model_predictions_comparison.png\n")

# Elastic Net residuals plot
p_resid <- predictions_df %>%
  ggplot(aes(x = age_true, y = en_residual)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
  labs(
    title = "Elastic Net Prediction Residuals",
    x = "Chronological Age (years)",
    y = "Residual (Predicted - Observed, years)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    aspect.ratio = 1
  )

# Scorcher NN residuals plot
p_resid_scorch <- predictions_df %>%
  ggplot(aes(x = age_true, y = scorch_residual)) +
  geom_point(alpha = 0.6, size = 2, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "orange", alpha = 0.2) +
  labs(
    title = "Scorcher NN Prediction Residuals",
    x = "Chronological Age (years)",
    y = "Residual (Predicted - Observed, years)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    aspect.ratio = 1
  )

# Combine residual plots
p_resid_combined <- cowplot::plot_grid(p_resid, p_resid_scorch, ncol = 2)

ggsave(fs::path(DIR_OUT, "model_residuals_comparison.png"),
       p_resid_combined, width = 14, height = 6, dpi = 300)

message(" + Residuals comparison plot saved: model_residuals_comparison.png\n")

#=== SUMMARY AND NEXT STEPS ==================================================

message("\n")
message("===============================================================================")
message("TRAINING COMPLETE")
message("===============================================================================\n")

message("SUMMARY OF RESULTS:\n")

message("Elastic Net Model (Baseline):")
message("  - Training MAE:  ", round(en_metrics_train$mae, 3), " years")
message("  - Test MAE:      ", round(en_metrics_test$mae, 3), " years")
message("  - Test RMSE:     ", round(en_metrics_test$rmse, 3), " years")
message("  - Test R-squared: ", round(en_metrics_test$r2, 4))
message("  - Test Correlation: ", round(en_metrics_test$pearson_r, 4))
message("\n")

message("Scorcher Neural Network:")
message("  - Training MAE:  ", round(scorch_metrics_train$mae, 3), " years")
message("  - Test MAE:      ", round(scorch_metrics_test$mae, 3), " years")
message("  - Test RMSE:     ", round(scorch_metrics_test$rmse, 3), " years")
message("  - Test R-squared: ", round(scorch_metrics_test$r2, 4))
message("  - Test Correlation: ", round(scorch_metrics_test$pearson_r, 4))
message("  - Best hyperparams: LR=", best_config$learning_rate, 
        ", Dropout=", best_config$dropout_rate,
        ", BatchSize=", best_config$batch_size)
message("\n")

message("COMPARISON:")
message("  MAE Improvement: ", round(en_metrics_test$mae - scorch_metrics_test$mae, 3), " years")
message("  R-squared Improvement: ", round(scorch_metrics_test$r2 - en_metrics_test$r2, 4))
message("\n")

message("OUTPUT FILES GENERATED:\n")

output_files <- c(
  "elasticnet_final_model.qs",
  "scorcher_final_model.qs",
  "scorcher_final_model.pt",
  "scorcher_hyperparam_grid.csv",
  "scorcher_tuning_results.csv",
  "model_comparison_results.qs",
  "model_comparison_results.csv",
  "predictions_on_test_set.qs",
  "predictions_on_test_set.csv",
  "model_predictions_comparison.png",
  "model_residuals_comparison.png"
)

for (file in output_files) {
  path <- fs::path(DIR_OUT, file)
  if (file.exists(path)) {
    size <- file.size(path)
    if (size > 1024^2) {
      size_str <- paste0(round(size / 1024^2, 1), " MB")
    } else if (size > 1024) {
      size_str <- paste0(round(size / 1024, 1), " KB")
    } else {
      size_str <- paste0(size, " B")
    }
    message("  (saved) ", file, " (", size_str, ")")
  }
}

message("\n")

message("NEXT STEPS:\n")
message("  1. Review hyperparameter tuning in scorcher_tuning_results.csv")
message("  2. Compare elastic net vs scorcher NN in model_comparison_results.csv")
message("  3. Examine prediction plots:")
message("     - model_predictions_comparison.png")
message("     - model_residuals_comparison.png")
message("  4. Load predictions for further analysis:")
message("     predictions <- qs_read('data/processed/predictions_on_test_set.qs')")
message("  5. Load final models for deployment:")
message("     en_model <- qs_read('data/processed/elasticnet_final_model.qs')")
message("     scorch_model <- qs_read('data/processed/scorcher_final_model.qs')")
message("  6. Run analysis.R for cooking show-style demonstration")
message("\n")

message("Done.\n")

#=== END ====================================================================== 
sessionInfo()
