#===============================================================================
#
#  PROGRAM: training.R
#
#  PURPOSE: Train a DeepMAge-style epigenetic aging clock from the sesame-based
#           preprocessing output produced by code/data.R.
#
#  This is written as a pedagogical script. The scorcher graph construction is
#  shown directly in each modeling section instead of hidden behind helper
#  functions.
#
#===============================================================================

message("\nDeepMAge-style scorcher training\n")

set.seed(1)

DIR_DATA <- "data/processed"
DIR_OUT <- "data/processed"

data_path <- file.path(DIR_DATA, "scorcher_train_test.qs")

if (!file.exists(data_path)) {
  stop("Missing ", data_path, ". Run code/data.R first.")
}

cran_pkgs <- c("qs2", "fs", "glmnet", "ggplot2", "readr", "tibble")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("torch", quietly = TRUE)) {
  install.packages("torch")
}

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

if (!requireNamespace("scorcher", quietly = TRUE)) {
  pak::pak("jtleek/scorcher@ssupdates")
}

library(qs2)
library(torch)
library(scorcher)

# Paper-reported final model settings.
n_features <- 1000L
paper_hidden_units <- 512L
paper_n_hidden_layers <- 4L
paper_activation <- "elu"
paper_optimizer <- "adam"
paper_dropout <- 0.30
paper_l2 <- 1e-3
paper_learning_rate <- 1e-4

# Runtime settings. The full DeepMAge grid is large, so demo mode evaluates a
# representative subset by default while still writing the complete grid.
batch_size <- 64L
feature_epochs <- 50L
tuning_epochs <- 50L
final_epochs <- 200L
cv_folds <- 5L
tuning_mode <- tolower(Sys.getenv("DEEPMAGE_TUNING_MODE", "demo"))
tuning_max_configs <- as.integer(Sys.getenv("DEEPMAGE_TUNING_MAX_CONFIGS", "12"))
if (is.na(tuning_max_configs)) tuning_max_configs <- 12L

#=== LOAD DATA =================================================================

message("Loading ", data_path)
obj <- qs2::qs_read(data_path)

x_train <- obj$x_train
x_test <- obj$x_test
y_train <- as.numeric(obj$y_train)
y_test <- as.numeric(obj$y_test)
pheno_train <- obj$pheno_train
pheno_test <- obj$pheno_test

if (!"is_control" %in% names(pheno_train)) {
  pheno_train$is_control <- NA
}

if (!"is_control" %in% names(pheno_test)) {
  pheno_test$is_control <- NA
}

keep_train <- !is.na(y_train)
keep_test <- !is.na(y_test)

x_train <- x_train[keep_train, , drop = FALSE]
y_train <- y_train[keep_train]
pheno_train <- pheno_train[keep_train, , drop = FALSE]

x_test <- x_test[keep_test, , drop = FALSE]
y_test <- y_test[keep_test]
pheno_test <- pheno_test[keep_test, , drop = FALSE]

message("Training data: ", nrow(x_train), " samples x ", ncol(x_train), " CpGs")
message("Test data:     ", nrow(x_test), " samples x ", ncol(x_test), " CpGs")

#=== STANDARDIZE FEATURES ======================================================

message("\nStandardizing CpGs using training-set means and SDs")

train_center <- colMeans(x_train, na.rm = TRUE)
train_scale <- apply(x_train, 2, stats::sd, na.rm = TRUE)
train_scale[is.na(train_scale) | train_scale == 0] <- 1

x_train_z <- sweep(x_train, 2, train_center, "-")
x_train_z <- sweep(x_train_z, 2, train_scale, "/")

x_test_z <- sweep(x_test, 2, train_center, "-")
x_test_z <- sweep(x_test_z, 2, train_scale, "/")

#=== FEATURE SELECTION MODEL ===================================================

message("\nSelecting 1,000 CpGs by input-gradient magnitude")
message("Fitting the initial scorcher model used only for feature ranking")

x_feature_tensor <- torch_tensor(as.matrix(x_train_z), dtype = torch_float())
y_feature_tensor <- torch_tensor(y_train, dtype = torch_float())$unsqueeze(2)

dl_feature <- scorch_create_dataloader(
  x_feature_tensor,
  y_feature_tensor,
  batch_size = batch_size
)

feature_fit <- initiate_scorch(dl_feature) |>
  scorch_input(features) |>
  scorch_layer(linear, in_features = ncol(x_train_z), out_features = 512) |>
  scorch_layer(elu) |>
  scorch_dropout(p = 0.30) |>
  scorch_layer(linear, in_features = 512, out_features = 512) |>
  scorch_layer(elu) |>
  scorch_dropout(p = 0.30) |>
  scorch_layer(linear, in_features = 512, out_features = 512) |>
  scorch_layer(elu) |>
  scorch_dropout(p = 0.30) |>
  scorch_layer(linear, in_features = 512, out_features = 512) |>
  scorch_layer(elu) |>
  scorch_dropout(p = 0.30) |>
  scorch_layer(linear, in_features = 512, out_features = 1, .name = prediction) |>
  scorch_output(prediction) |>
  compile_scorch(
    loss_fn = nn_l1_loss(reduction = "mean"),
    optimizer_fn = optim_adam,
    optimizer_params = list(
      lr = paper_learning_rate,
      weight_decay = paper_l2
    )
  ) |>
  fit_scorch(num_epochs = feature_epochs, seed = 1L)

feature_fit$nn_model$eval()
x_gradient_tensor <- torch_tensor(as.matrix(x_train_z), dtype = torch_float())
x_gradient_tensor$requires_grad_(TRUE)

gradient_predictions <- feature_fit$nn_model(x_gradient_tensor)
gradient_predictions$sum()$backward()

gradient_matrix <- as_array(x_gradient_tensor$grad)
gradient_importance <- colMeans(abs(gradient_matrix), na.rm = TRUE)
feature_rank <- order(gradient_importance, decreasing = TRUE, na.last = NA)
feature_idx <- feature_rank[seq_len(min(n_features, length(feature_rank)))]
selected_features <- colnames(x_train_z)[feature_idx]

feature_tbl <- tibble::tibble(
  cpg = selected_features,
  gradient_importance = gradient_importance[feature_idx],
  rank = seq_along(selected_features)
)

readr::write_csv(feature_tbl, fs::path(DIR_OUT, "deepmage_selected_cpgs.csv"))
qs2::qs_save(feature_tbl, fs::path(DIR_OUT, "deepmage_selected_cpgs.qs"))

x_train_sel <- x_train_z[, selected_features, drop = FALSE]
x_test_sel <- x_test_z[, selected_features, drop = FALSE]

message("Selected CpGs: ", ncol(x_train_sel))

#=== HYPERPARAMETER GRID =======================================================

message("\nBuilding DeepMAge paper-range hyperparameter grid")

paper_grid <- expand.grid(
  n_hidden_layers = 2:5,
  hidden_units = seq(128L, 1024L, by = 128L),
  activation = c("elu", "relu", "selu"),
  optimizer = c("adam", "amsgrad", "nadam"),
  dropout = seq(0.15, 0.50, by = 0.05),
  l2 = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
  learning_rate = paper_learning_rate,
  stringsAsFactors = FALSE
)

if (tuning_mode == "skip") {
  tuning_grid <- paper_grid[0, , drop = FALSE]
} else if (tuning_mode == "full") {
  tuning_grid <- paper_grid
} else {
  best_default <- which(
    paper_grid$n_hidden_layers == paper_n_hidden_layers &
      paper_grid$hidden_units == paper_hidden_units &
      paper_grid$activation == paper_activation &
      paper_grid$optimizer == paper_optimizer &
      abs(paper_grid$dropout - paper_dropout) < 1e-8 &
      abs(paper_grid$l2 - paper_l2) < 1e-12
  )

  n_demo <- min(max(tuning_max_configs, 1L), nrow(paper_grid))
  set.seed(1)
  demo_idx <- unique(c(best_default[1], sample(seq_len(nrow(paper_grid)), n_demo - 1L)))
  tuning_grid <- paper_grid[demo_idx, , drop = FALSE]
}

readr::write_csv(paper_grid, fs::path(DIR_OUT, "deepmage_full_hyperparameter_grid.csv"))
readr::write_csv(tuning_grid, fs::path(DIR_OUT, "deepmage_tuning_grid.csv"))

message("Full grid size: ", nrow(paper_grid), " configurations")
message("Tuning mode: ", tuning_mode)
message("Configurations to evaluate now: ", nrow(tuning_grid))

#=== HYPERPARAMETER TUNING =====================================================

tuning_results <- tuning_grid

if (nrow(tuning_grid) > 0) {
  tuning_results$n_folds <- NA_integer_
  tuning_results$mean_mae <- NA_real_
  tuning_results$mean_medae <- NA_real_
  tuning_results$status <- NA_character_

  set.seed(1)
  k <- min(cv_folds, nrow(x_train_sel))
  fold_id <- sample(rep(seq_len(k), length.out = nrow(x_train_sel)))

  for (config_i in seq_len(nrow(tuning_grid))) {
    current <- tuning_grid[config_i, , drop = FALSE]

    current_n_hidden_layers <- as.integer(current$n_hidden_layers)
    current_hidden_units <- as.integer(current$hidden_units)
    current_activation <- as.character(current$activation)
    current_optimizer <- as.character(current$optimizer)
    current_dropout <- as.numeric(current$dropout)
    current_l2 <- as.numeric(current$l2)
    current_learning_rate <- as.numeric(current$learning_rate)

    message(
      "  config ", config_i, "/", nrow(tuning_grid),
      ": layers=", current_n_hidden_layers,
      ", units=", current_hidden_units,
      ", activation=", current_activation,
      ", optimizer=", current_optimizer,
      ", dropout=", current_dropout,
      ", l2=", current_l2
    )

    current_optimizer_params <- list(
      lr = current_learning_rate,
      weight_decay = current_l2
    )

    current_optimizer_fn <- optim_adam
    current_status <- "ok"

    if (current_optimizer == "amsgrad") {
      current_optimizer_params$amsgrad <- TRUE
    }

    if (current_optimizer == "nadam") {
      if (exists("optim_nadam", envir = asNamespace("torch"))) {
        current_optimizer_fn <- get("optim_nadam", envir = asNamespace("torch"))
      } else {
        current_status <- "skipped_optimizer_unavailable"
      }
    }

    fold_mae <- rep(NA_real_, k)
    fold_medae <- rep(NA_real_, k)

    if (current_status == "ok") {
      for (fold in seq_len(k)) {
        val_idx <- which(fold_id == fold)
        trn_idx <- which(fold_id != fold)

        x_fold_train_tensor <- torch_tensor(
          as.matrix(x_train_sel[trn_idx, , drop = FALSE]),
          dtype = torch_float()
        )
        y_fold_train_tensor <- torch_tensor(
          y_train[trn_idx],
          dtype = torch_float()
        )$unsqueeze(2)

        dl_fold <- scorch_create_dataloader(
          x_fold_train_tensor,
          y_fold_train_tensor,
          batch_size = batch_size
        )

        fold_model <- initiate_scorch(dl_fold) |>
          scorch_input(features)

        in_features <- ncol(x_train_sel)

        for (layer_i in seq_len(current_n_hidden_layers)) {
          fold_model <- fold_model |>
            scorch_layer(
              linear,
              in_features = in_features,
              out_features = current_hidden_units
            )

          if (current_activation == "elu") {
            fold_model <- fold_model |>
              scorch_layer(elu)
          } else if (current_activation == "relu") {
            fold_model <- fold_model |>
              scorch_layer(relu)
          } else if (current_activation == "selu") {
            fold_model <- fold_model |>
              scorch_layer(selu)
          }

          fold_model <- fold_model |>
            scorch_dropout(p = current_dropout)

          in_features <- current_hidden_units
        }

        fold_fit <- tryCatch({
          fold_model |>
            scorch_layer(
              linear,
              in_features = in_features,
              out_features = 1,
              .name = prediction
            ) |>
            scorch_output(prediction) |>
            compile_scorch(
              loss_fn = nn_l1_loss(reduction = "mean"),
              optimizer_fn = current_optimizer_fn,
              optimizer_params = current_optimizer_params
            ) |>
            fit_scorch(num_epochs = tuning_epochs, seed = 1L)
        }, error = function(e) {
          current_status <<- conditionMessage(e)
          NULL
        })

        if (is.null(fold_fit)) {
          break
        }

        fold_fit$nn_model$eval()
        x_fold_val_tensor <- torch_tensor(
          as.matrix(x_train_sel[val_idx, , drop = FALSE]),
          dtype = torch_float()
        )

        with_no_grad({
          fold_pred <- as.numeric(fold_fit$nn_model(x_fold_val_tensor)$squeeze())
        })

        fold_abs_error <- abs(fold_pred - y_train[val_idx])
        fold_mae[fold] <- mean(fold_abs_error, na.rm = TRUE)
        fold_medae[fold] <- stats::median(fold_abs_error, na.rm = TRUE)

        rm(fold_fit)
        gc()
      }
    }

    tuning_results$n_folds[config_i] <- sum(!is.na(fold_mae))
    tuning_results$mean_mae[config_i] <- mean(fold_mae, na.rm = TRUE)
    tuning_results$mean_medae[config_i] <- mean(fold_medae, na.rm = TRUE)
    tuning_results$status[config_i] <- current_status
  }
}

readr::write_csv(tuning_results, fs::path(DIR_OUT, "deepmage_tuning_results.csv"))
qs2::qs_save(tuning_results, fs::path(DIR_OUT, "deepmage_tuning_results.qs"))

final_n_hidden_layers <- paper_n_hidden_layers
final_hidden_units <- paper_hidden_units
final_activation <- paper_activation
final_optimizer <- paper_optimizer
final_dropout <- paper_dropout
final_l2 <- paper_l2
final_learning_rate <- paper_learning_rate

if ("status" %in% names(tuning_results)) {
  usable_tuning <- tuning_results[
    tuning_results$status == "ok" & is.finite(tuning_results$mean_mae),
    ,
    drop = FALSE
  ]
} else {
  usable_tuning <- tuning_results[0, , drop = FALSE]
}

if (nrow(usable_tuning) > 0) {
  best_i <- which.min(usable_tuning$mean_mae)
  best_row <- usable_tuning[best_i, , drop = FALSE]

  final_n_hidden_layers <- as.integer(best_row$n_hidden_layers)
  final_hidden_units <- as.integer(best_row$hidden_units)
  final_activation <- as.character(best_row$activation)
  final_optimizer <- as.character(best_row$optimizer)
  final_dropout <- as.numeric(best_row$dropout)
  final_l2 <- as.numeric(best_row$l2)
  final_learning_rate <- as.numeric(best_row$learning_rate)

  message(
    "Best tuned config: layers=", final_n_hidden_layers,
    ", units=", final_hidden_units,
    ", activation=", final_activation,
    ", optimizer=", final_optimizer,
    ", dropout=", final_dropout,
    ", l2=", final_l2,
    ", CV MAE=", round(best_row$mean_mae, 3)
  )
} else {
  message("No usable tuning result; using paper's final reported configuration.")
}

#=== ELASTIC-NET BASELINE ======================================================

message("\nTraining elastic-net baseline on selected CpGs")

set.seed(1)
cv_en <- glmnet::cv.glmnet(
  as.matrix(x_train_sel),
  y_train,
  alpha = 0.5,
  nfolds = cv_folds,
  type.measure = "mae",
  standardize = FALSE
)

en_pred_train <- as.numeric(predict(cv_en, as.matrix(x_train_sel), s = "lambda.1se"))
en_pred_test <- as.numeric(predict(cv_en, as.matrix(x_test_sel), s = "lambda.1se"))

#=== FINAL SCORCHER MODEL ======================================================

message("\nTraining final scorcher DeepMAge-style model")

final_optimizer_params <- list(
  lr = final_learning_rate,
  weight_decay = final_l2
)

final_optimizer_fn <- optim_adam

if (final_optimizer == "amsgrad") {
  final_optimizer_params$amsgrad <- TRUE
}

if (final_optimizer == "nadam" && exists("optim_nadam", envir = asNamespace("torch"))) {
  final_optimizer_fn <- get("optim_nadam", envir = asNamespace("torch"))
}

x_final_train_tensor <- torch_tensor(as.matrix(x_train_sel), dtype = torch_float())
y_final_train_tensor <- torch_tensor(y_train, dtype = torch_float())$unsqueeze(2)

dl_final <- scorch_create_dataloader(
  x_final_train_tensor,
  y_final_train_tensor,
  batch_size = batch_size
)

deepmage_model <- initiate_scorch(dl_final) |>
  scorch_input(features)

in_features <- ncol(x_train_sel)

for (layer_i in seq_len(final_n_hidden_layers)) {
  deepmage_model <- deepmage_model |>
    scorch_layer(
      linear,
      in_features = in_features,
      out_features = final_hidden_units
    )

  if (final_activation == "elu") {
    deepmage_model <- deepmage_model |>
      scorch_layer(elu)
  } else if (final_activation == "relu") {
    deepmage_model <- deepmage_model |>
      scorch_layer(relu)
  } else if (final_activation == "selu") {
    deepmage_model <- deepmage_model |>
      scorch_layer(selu)
  }

  deepmage_model <- deepmage_model |>
    scorch_dropout(p = final_dropout)

  in_features <- final_hidden_units
}

deepmage_fit <- deepmage_model |>
  scorch_layer(linear, in_features = in_features, out_features = 1, .name = prediction) |>
  scorch_output(prediction) |>
  compile_scorch(
    loss_fn = nn_l1_loss(reduction = "mean"),
    optimizer_fn = final_optimizer_fn,
    optimizer_params = final_optimizer_params
  ) |>
  fit_scorch(num_epochs = final_epochs, seed = 1L)

deepmage_fit$nn_model$eval()

with_no_grad({
  dm_pred_train <- as.numeric(deepmage_fit$nn_model(x_final_train_tensor)$squeeze())

  x_final_test_tensor <- torch_tensor(as.matrix(x_test_sel), dtype = torch_float())
  dm_pred_test <- as.numeric(deepmage_fit$nn_model(x_final_test_tensor)$squeeze())
})

#=== METRICS AND OUTPUTS =======================================================

metrics <- tibble::tibble()

metric_sets <- list(
  list(model = "Elastic Net", set = "train", truth = y_train, pred = en_pred_train, pheno = pheno_train),
  list(model = "Elastic Net", set = "test", truth = y_test, pred = en_pred_test, pheno = pheno_test),
  list(model = "Scorcher DeepMAge", set = "train", truth = y_train, pred = dm_pred_train, pheno = pheno_train),
  list(model = "Scorcher DeepMAge", set = "test", truth = y_test, pred = dm_pred_test, pheno = pheno_test)
)

for (metric_i in seq_along(metric_sets)) {
  current_metric <- metric_sets[[metric_i]]

  idx <- seq_along(current_metric$truth)
  truth_i <- current_metric$truth[idx]
  pred_i <- current_metric$pred[idx]

  metrics <- rbind(
    metrics,
    tibble::tibble(
      model = current_metric$model,
      set = current_metric$set,
      cohort = "all",
      n = length(truth_i),
      medae = stats::median(abs(pred_i - truth_i), na.rm = TRUE),
      mae = mean(abs(pred_i - truth_i), na.rm = TRUE),
      rmse = sqrt(mean((pred_i - truth_i)^2, na.rm = TRUE)),
      pearson_r = stats::cor(pred_i, truth_i, use = "complete.obs"),
      r2 = 1 - sum((truth_i - pred_i)^2, na.rm = TRUE) /
        sum((truth_i - mean(truth_i, na.rm = TRUE))^2, na.rm = TRUE)
    )
  )

  idx <- which(current_metric$pheno$is_control %in% TRUE)

  if (length(idx) > 1) {
    truth_i <- current_metric$truth[idx]
    pred_i <- current_metric$pred[idx]

    metrics <- rbind(
      metrics,
      tibble::tibble(
        model = current_metric$model,
        set = current_metric$set,
        cohort = "controls",
        n = length(truth_i),
        medae = stats::median(abs(pred_i - truth_i), na.rm = TRUE),
        mae = mean(abs(pred_i - truth_i), na.rm = TRUE),
        rmse = sqrt(mean((pred_i - truth_i)^2, na.rm = TRUE)),
        pearson_r = stats::cor(pred_i, truth_i, use = "complete.obs"),
        r2 = 1 - sum((truth_i - pred_i)^2, na.rm = TRUE) /
          sum((truth_i - mean(truth_i, na.rm = TRUE))^2, na.rm = TRUE)
      )
    )
  }
}

predictions <- rbind(
  data.frame(
    model = "Elastic Net",
    set = "train",
    sample_id = pheno_train$sample_id,
    age = y_train,
    predicted_age = en_pred_train,
    is_control = pheno_train$is_control
  ),
  data.frame(
    model = "Elastic Net",
    set = "test",
    sample_id = pheno_test$sample_id,
    age = y_test,
    predicted_age = en_pred_test,
    is_control = pheno_test$is_control
  ),
  data.frame(
    model = "Scorcher DeepMAge",
    set = "train",
    sample_id = pheno_train$sample_id,
    age = y_train,
    predicted_age = dm_pred_train,
    is_control = pheno_train$is_control
  ),
  data.frame(
    model = "Scorcher DeepMAge",
    set = "test",
    sample_id = pheno_test$sample_id,
    age = y_test,
    predicted_age = dm_pred_test,
    is_control = pheno_test$is_control
  )
)

predictions$residual <- predictions$predicted_age - predictions$age

final_config <- list(
  n_hidden_layers = final_n_hidden_layers,
  hidden_units = final_hidden_units,
  activation = final_activation,
  optimizer = final_optimizer,
  dropout = final_dropout,
  l2 = final_l2,
  learning_rate = final_learning_rate,
  batch_size = batch_size,
  feature_epochs = feature_epochs,
  tuning_epochs = tuning_epochs,
  final_epochs = final_epochs,
  cv_folds = cv_folds
)

training_object <- list(
  config = final_config,
  full_hyperparameter_grid = paper_grid,
  tuning_grid = tuning_grid,
  tuning_results = tuning_results,
  feature_model = feature_fit,
  final_model = deepmage_fit,
  elastic_net = cv_en,
  selected_features = selected_features,
  feature_importance = feature_tbl,
  train_center = train_center,
  train_scale = train_scale,
  metrics = metrics,
  predictions = predictions
)

qs2::qs_save(training_object, fs::path(DIR_OUT, "deepmage_scorcher_training.qs"))
qs2::qs_save(cv_en, fs::path(DIR_OUT, "elasticnet_final_model.qs"))
qs2::qs_save(deepmage_fit, fs::path(DIR_OUT, "scorcher_deepmage_model.qs"))
torch::torch_save(deepmage_fit$nn_model, fs::path(DIR_OUT, "scorcher_deepmage_model.pt"))

readr::write_csv(metrics, fs::path(DIR_OUT, "model_comparison_results.csv"))
readr::write_csv(predictions, fs::path(DIR_OUT, "predictions_all_samples.csv"))
qs2::qs_save(metrics, fs::path(DIR_OUT, "model_comparison_results.qs"))
qs2::qs_save(predictions, fs::path(DIR_OUT, "predictions_all_samples.qs"))

plot_df <- predictions[predictions$set == "test", , drop = FALSE]

p <- ggplot2::ggplot(plot_df, ggplot2::aes(age, predicted_age, color = model)) +
  ggplot2::geom_point(alpha = 0.65, size = 1.8) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Chronological age",
    y = "Predicted age",
    color = NULL
  )

ggplot2::ggsave(
  fs::path(DIR_OUT, "test_predictions_deepmage_scorcher.png"),
  p,
  width = 7,
  height = 6,
  dpi = 300
)

message("\nHeadline controls-only metrics")
print(metrics[metrics$cohort == "controls", , drop = FALSE])

message("\nDone.")
sessionInfo()
