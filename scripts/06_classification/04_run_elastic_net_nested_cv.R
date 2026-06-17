# ============================================================
# Elastic Net nested CV classifier for cfDNA-MS feature sets
#
# Preserved from original workflow:
#   - Same 5 binary tasks
#   - Outer 5-fold CV
#   - Inner 3-fold CV
#   - alpha/lambda tuning inside inner CV
#   - threshold selection from inner CV
#   - final test evaluation only on outer test fold
#
# Usage:
#   Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R bioanalyzer_10bp_length
#   Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R regional_mds_167
#   Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R global_endmotif_256
# ============================================================

required_packages <- c(
  "glmnet",
  "caret",
  "pROC",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "stringr",
  "readr"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(pROC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(readr)
})

set.seed(111)
options(stringsAsFactors = FALSE)
options(mc.cores = 1)

# ============================================================
# 0. Input and output
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Please provide a feature set name.\n",
    "Example:\n",
    "Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R bioanalyzer_10bp_length"
  )
}

feature_set_name <- args[1]

feature_file <- file.path(
  "results/intermediate/classification/feature_sets",
  paste0(feature_set_name, ".tsv")
)

if (!file.exists(feature_file)) {
  stop("Feature file not found: ", feature_file)
}

table_out_dir <- file.path(
  "results/tables/classification/elastic_net",
  feature_set_name
)

figure_out_dir <- file.path(
  "results/figures/classification/elastic_net",
  feature_set_name
)

dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("Elastic Net nested CV classification\n")
cat("Feature set:", feature_set_name, "\n")
cat("Feature file:", feature_file, "\n")
cat("Table output:", table_out_dir, "\n")
cat("Figure output:", figure_out_dir, "\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load feature matrix
# ============================================================

feature_df <- read_tsv(feature_file, show_col_types = FALSE)

if (!all(c("sample", "group") %in% colnames(feature_df))) {
  stop("Feature file must contain columns: sample and group")
}

feature_df <- feature_df %>%
  mutate(
    sample = as.character(sample),
    group = factor(
      as.character(group),
      levels = c("Healthy", "Remission", "Relapse")
    )
  ) %>%
  arrange(group, sample)

feature_cols <- setdiff(colnames(feature_df), c("sample", "group"))

feature_matrix <- feature_df %>%
  select(all_of(feature_cols)) %>%
  as.data.frame()

feature_matrix[] <- lapply(feature_matrix, as.numeric)
colnames(feature_matrix) <- make.names(colnames(feature_matrix), unique = TRUE)
rownames(feature_matrix) <- feature_df$sample

groups <- feature_df$group

cat("Feature matrix dimensions:\n")
print(dim(feature_matrix))

cat("\nGroup counts:\n")
print(table(groups))

cat("\nMissing values in full feature matrix:\n")
print(sum(is.na(feature_matrix)))

write_csv(
  feature_df %>% select(sample, group),
  file.path(table_out_dir, "sample_group_assignment_check.csv")
)

# ============================================================
# 2. Tasks
# ============================================================

tasks <- list(
  list(name = "Healthy vs MS",          case = c("Remission", "Relapse"), control = "Healthy"),
  list(name = "Healthy vs Remission",   case = "Remission",               control = "Healthy"),
  list(name = "Healthy vs Relapse",     case = "Relapse",                 control = "Healthy"),
  list(name = "NonRelapse vs Relapse",  case = "Relapse",                 control = c("Healthy", "Remission")),
  list(name = "Remission vs Relapse",   case = "Relapse",                 control = "Remission")
)

method_label <- "Elastic Net"
method_color <- c("Elastic Net" = "#E63946")

# ============================================================
# 3. Settings
# ============================================================

alpha_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
lambda_grid <- 10^seq(-4, 1, length.out = 80)

n_outer_folds <- 5
n_inner_folds <- 3

use_class_weights <- FALSE

# ============================================================
# 4. Helper functions
# ============================================================

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

filter_train_test_features <- function(X_train, X_test,
                                       max_missing_fraction = 0.10,
                                       min_sd = 1e-8) {

  n_train <- nrow(X_train)

  keep_missing <- colSums(is.na(X_train)) <= floor(max_missing_fraction * n_train)

  X_train_f <- X_train[, keep_missing, drop = FALSE]
  X_test_f  <- X_test[,  keep_missing, drop = FALSE]

  if (ncol(X_train_f) == 0) {
    stop("No features left after missingness filtering.")
  }

  feature_sd <- apply(X_train_f, 2, sd, na.rm = TRUE)
  keep_sd <- !is.na(feature_sd) & feature_sd > min_sd

  X_train_f <- X_train_f[, keep_sd, drop = FALSE]
  X_test_f  <- X_test_f[,  keep_sd, drop = FALSE]

  if (ncol(X_train_f) == 0) {
    stop("No features left after variance filtering.")
  }

  list(
    X_train = X_train_f,
    X_test = X_test_f,
    kept_features = colnames(X_train_f)
  )
}

impute_train_test_median <- function(X_train, X_test) {

  medians <- apply(X_train, 2, function(x) median(x, na.rm = TRUE))
  medians[is.na(medians)] <- 0

  for (j in seq_along(medians)) {
    X_train[is.na(X_train[, j]), j] <- medians[j]
    X_test[is.na(X_test[, j]), j]  <- medians[j]
  }

  list(
    X_train = X_train,
    X_test = X_test,
    medians = medians
  )
}

make_class_weights <- function(y_bin) {
  n_case <- sum(y_bin == 1)
  n_control <- sum(y_bin == 0)
  n_total <- length(y_bin)

  w_case <- n_total / (2 * n_case)
  w_control <- n_total / (2 * n_control)

  ifelse(y_bin == 1, w_case, w_control)
}

safe_roc_auc_numeric <- function(response, predictor) {

  if (length(unique(response[!is.na(response)])) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(roc_obj$auc)
}

safe_roc_obj_numeric <- function(response, predictor) {

  pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

# ============================================================
# 5. Inner CV: tune alpha, lambda, and threshold
# ============================================================

tune_elastic_net_inner_cv <- function(X_outer_train_raw,
                                      y_outer_train,
                                      n_inner_folds = 3,
                                      alpha_grid,
                                      lambda_grid,
                                      use_class_weights = FALSE) {

  y_outer_train <- factor(y_outer_train, levels = c("control", "case"))

  inner_folds <- createFolds(
    y_outer_train,
    k = n_inner_folds,
    returnTrain = TRUE
  )

  tuning_rows <- list()
  prediction_store <- list()

  for (alpha_value in alpha_grid) {

    for (lambda_value in lambda_grid) {

      inner_predictions <- list()

      for (j in seq_along(inner_folds)) {

        inner_train_idx <- inner_folds[[j]]
        inner_valid_idx <- setdiff(seq_len(nrow(X_outer_train_raw)), inner_train_idx)

        X_inner_train_raw <- X_outer_train_raw[inner_train_idx, , drop = FALSE]
        X_inner_valid_raw <- X_outer_train_raw[inner_valid_idx, , drop = FALSE]

        y_inner_train <- y_outer_train[inner_train_idx]
        y_inner_valid <- y_outer_train[inner_valid_idx]

        y_inner_train_bin <- ifelse(y_inner_train == "case", 1, 0)
        y_inner_valid_bin <- ifelse(y_inner_valid == "case", 1, 0)

        imp_inner <- impute_train_test_median(
          X_inner_train_raw,
          X_inner_valid_raw
        )

        X_inner_train <- as.matrix(imp_inner$X_train)
        X_inner_valid <- as.matrix(imp_inner$X_test)

        weights_inner <- NULL
        if (use_class_weights) {
          weights_inner <- make_class_weights(y_inner_train_bin)
        }

        enet_model <- tryCatch(
          glmnet(
            x = X_inner_train,
            y = y_inner_train_bin,
            family = "binomial",
            alpha = alpha_value,
            lambda = lambda_value,
            standardize = TRUE,
            weights = weights_inner
          ),
          error = function(e) NULL
        )

        if (is.null(enet_model)) {
          next
        }

        pred_prob <- as.numeric(
          predict(
            enet_model,
            newx = X_inner_valid,
            type = "response",
            s = lambda_value
          )
        )

        inner_predictions[[length(inner_predictions) + 1]] <- data.frame(
          Inner_Fold = j,
          Alpha = alpha_value,
          Lambda = lambda_value,
          True_Binary = y_inner_valid_bin,
          True_Class = as.character(y_inner_valid),
          Predicted_Probability_Case = pred_prob,
          stringsAsFactors = FALSE
        )
      }

      pred_param <- bind_rows(inner_predictions)

      auc_param <- NA_real_

      if (nrow(pred_param) > 0 &&
          length(unique(pred_param$True_Binary)) == 2) {

        auc_param <- safe_roc_auc_numeric(
          response = pred_param$True_Binary,
          predictor = pred_param$Predicted_Probability_Case
        )
      }

      key <- paste(alpha_value, lambda_value, sep = "_")
      prediction_store[[key]] <- pred_param

      tuning_rows[[length(tuning_rows) + 1]] <- data.frame(
        Alpha = alpha_value,
        Lambda = lambda_value,
        Inner_AUC = auc_param,
        stringsAsFactors = FALSE
      )
    }
  }

  tuning_table <- bind_rows(tuning_rows)

  if (all(is.na(tuning_table$Inner_AUC))) {
    return(list(
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      best_threshold = NA_real_,
      tuning_table = tuning_table,
      inner_predictions_best = data.frame()
    ))
  }

  best_row <- tuning_table %>%
    arrange(desc(Inner_AUC), desc(Alpha), desc(Lambda)) %>%
    slice(1)

  best_alpha <- best_row$Alpha[1]
  best_lambda <- best_row$Lambda[1]

  best_key <- paste(best_alpha, best_lambda, sep = "_")
  pred_best <- prediction_store[[best_key]]

  inner_roc <- safe_roc_obj_numeric(
    response = pred_best$True_Binary,
    predictor = pred_best$Predicted_Probability_Case
  )

  best_threshold <- coords(
    inner_roc,
    x = "best",
    best.method = "youden",
    transpose = FALSE
  )$threshold[1]

  list(
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_threshold = as.numeric(best_threshold),
    tuning_table = tuning_table,
    inner_predictions_best = pred_best
  )
}

# ============================================================
# 6. Outer CV
# ============================================================

run_nested_cv_elastic_net <- function(data,
                                      y_binary,
                                      task_name,
                                      n_outer_folds = 5,
                                      n_inner_folds = 3,
                                      alpha_grid,
                                      lambda_grid,
                                      use_class_weights = FALSE) {

  y_binary <- factor(y_binary, levels = c("control", "case"))

  outer_folds <- createFolds(
    y_binary,
    k = n_outer_folds,
    returnTrain = TRUE
  )

  fold_metrics <- list()
  prediction_rows <- list()
  selected_feature_rows <- list()
  tuning_rows_all <- list()

  for (i in seq_along(outer_folds)) {

    cat(sprintf("      Outer fold %d/%d\n", i, length(outer_folds)))

    outer_train_idx <- outer_folds[[i]]
    outer_test_idx  <- setdiff(seq_len(nrow(data)), outer_train_idx)

    X_outer_train_raw0 <- data[outer_train_idx, , drop = FALSE]
    X_outer_test_raw0  <- data[outer_test_idx,  , drop = FALSE]

    y_outer_train <- y_binary[outer_train_idx]
    y_outer_test  <- y_binary[outer_test_idx]

    y_outer_train_bin <- ifelse(y_outer_train == "case", 1, 0)
    y_outer_test_bin  <- ifelse(y_outer_test  == "case", 1, 0)

    filt <- filter_train_test_features(
      X_train = X_outer_train_raw0,
      X_test = X_outer_test_raw0,
      max_missing_fraction = 0.10,
      min_sd = 1e-8
    )

    X_outer_train_raw <- filt$X_train
    X_outer_test_raw  <- filt$X_test

    tune_res <- tune_elastic_net_inner_cv(
      X_outer_train_raw = X_outer_train_raw,
      y_outer_train = y_outer_train,
      n_inner_folds = n_inner_folds,
      alpha_grid = alpha_grid,
      lambda_grid = lambda_grid,
      use_class_weights = use_class_weights
    )

    best_alpha <- tune_res$best_alpha
    best_lambda <- tune_res$best_lambda
    best_threshold <- tune_res$best_threshold

    tuning_table <- tune_res$tuning_table %>%
      mutate(
        Feature_Set = feature_set_name,
        Task = task_name,
        Fold = i,
        Method = method_label
      )

    tuning_rows_all[[length(tuning_rows_all) + 1]] <- tuning_table

    if (is.na(best_alpha) || is.na(best_lambda)) {

      fold_metrics[[length(fold_metrics) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Task = task_name,
        Fold = i,
        Method = method_label,
        N_Train = length(outer_train_idx),
        N_Test = length(outer_test_idx),
        N_Train_Control = sum(y_outer_train == "control"),
        N_Train_Case = sum(y_outer_train == "case"),
        N_Test_Control = sum(y_outer_test == "control"),
        N_Test_Case = sum(y_outer_test == "case"),
        N_Features_After_Filtering = ncol(X_outer_train_raw),
        Best_Alpha = NA_real_,
        Best_Lambda = NA_real_,
        Threshold = NA_real_,
        Accuracy = NA_real_,
        Balanced_Accuracy = NA_real_,
        Sensitivity = NA_real_,
        Specificity = NA_real_,
        PPV = NA_real_,
        NPV = NA_real_,
        F1 = NA_real_,
        AUC = NA_real_,
        N_Selected_Features = NA_real_,
        stringsAsFactors = FALSE
      )

      next
    }

    imp_outer <- impute_train_test_median(
      X_train = X_outer_train_raw,
      X_test = X_outer_test_raw
    )

    X_outer_train <- as.matrix(imp_outer$X_train)
    X_outer_test  <- as.matrix(imp_outer$X_test)

    weights_outer <- NULL
    if (use_class_weights) {
      weights_outer <- make_class_weights(y_outer_train_bin)
    }

    final_model <- glmnet(
      x = X_outer_train,
      y = y_outer_train_bin,
      family = "binomial",
      alpha = best_alpha,
      lambda = best_lambda,
      standardize = TRUE,
      weights = weights_outer
    )

    pred_prob <- as.numeric(
      predict(
        final_model,
        newx = X_outer_test,
        type = "response",
        s = best_lambda
      )
    )

    pred_class <- factor(
      ifelse(pred_prob >= best_threshold, "case", "control"),
      levels = c("case", "control")
    )

    y_test_fac <- factor(
      y_outer_test,
      levels = c("case", "control")
    )

    cm <- confusionMatrix(
      data = pred_class,
      reference = y_test_fac,
      positive = "case"
    )

    outer_auc <- safe_roc_auc_numeric(
      response = y_outer_test_bin,
      predictor = pred_prob
    )

    coef_mat <- coef(final_model, s = best_lambda)
    coef_vec <- as.numeric(coef_mat)
    coef_names <- rownames(coef_mat)

    selected_idx <- which(coef_vec != 0 & coef_names != "(Intercept)")
    selected_features <- coef_names[selected_idx]
    selected_coefs <- coef_vec[selected_idx]

    fold_metrics[[length(fold_metrics) + 1]] <- data.frame(
      Feature_Set = feature_set_name,
      Task = task_name,
      Fold = i,
      Method = method_label,
      N_Train = length(outer_train_idx),
      N_Test = length(outer_test_idx),
      N_Train_Control = sum(y_outer_train == "control"),
      N_Train_Case = sum(y_outer_train == "case"),
      N_Test_Control = sum(y_outer_test == "control"),
      N_Test_Case = sum(y_outer_test == "case"),
      N_Features_After_Filtering = ncol(X_outer_train),
      Best_Alpha = best_alpha,
      Best_Lambda = best_lambda,
      Threshold = best_threshold,
      Accuracy = unname(cm$overall["Accuracy"]),
      Balanced_Accuracy = unname(cm$byClass["Balanced Accuracy"]),
      Sensitivity = unname(cm$byClass["Sensitivity"]),
      Specificity = unname(cm$byClass["Specificity"]),
      PPV = unname(cm$byClass["Pos Pred Value"]),
      NPV = unname(cm$byClass["Neg Pred Value"]),
      F1 = unname(cm$byClass["F1"]),
      AUC = outer_auc,
      N_Selected_Features = length(selected_features),
      stringsAsFactors = FALSE
    )

    prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
      Feature_Set = feature_set_name,
      Task = task_name,
      Fold = i,
      Method = method_label,
      Sample = rownames(data)[outer_test_idx],
      True_Class = as.character(y_outer_test),
      True_Binary = y_outer_test_bin,
      Predicted_Probability_Case = pred_prob,
      Predicted_Class = as.character(pred_class),
      stringsAsFactors = FALSE
    )

    if (length(selected_features) > 0) {
      selected_feature_rows[[length(selected_feature_rows) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Task = task_name,
        Fold = i,
        Method = method_label,
        Feature = selected_features,
        Coefficient = selected_coefs,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    fold_metrics = bind_rows(fold_metrics),
    predictions = bind_rows(prediction_rows),
    selected_features = bind_rows(selected_feature_rows),
    tuning_table = bind_rows(tuning_rows_all)
  )
}

# ============================================================
# 7. Run all tasks
# ============================================================

cat("\n============================================================\n")
cat("RUNNING ELASTIC NET NESTED CV\n")
cat("Feature set:", feature_set_name, "\n")
cat("============================================================\n\n")

all_task_results <- list()

for (task in tasks) {

  cat(sprintf(">> Task: %s\n", task$name))

  control_idx <- which(groups %in% task$control)
  case_idx <- which(groups %in% task$case)

  keep_idx <- c(control_idx, case_idx)

  task_data <- feature_matrix[keep_idx, , drop = FALSE]

  task_y_binary <- factor(
    c(rep("control", length(control_idx)), rep("case", length(case_idx))),
    levels = c("control", "case")
  )

  rownames(task_data) <- rownames(feature_matrix)[keep_idx]

  res <- run_nested_cv_elastic_net(
    data = task_data,
    y_binary = task_y_binary,
    task_name = task$name,
    n_outer_folds = n_outer_folds,
    n_inner_folds = n_inner_folds,
    alpha_grid = alpha_grid,
    lambda_grid = lambda_grid,
    use_class_weights = use_class_weights
  )

  all_task_results[[task$name]] <- list(
    n_control = length(control_idx),
    n_case = length(case_idx),
    results = res
  )

  cat("\n")
}

# ============================================================
# 8. Combine and save detailed results
# ============================================================

all_fold_metrics <- bind_rows(
  lapply(names(all_task_results), function(task_name) {
    all_task_results[[task_name]]$results$fold_metrics
  })
)

all_predictions <- bind_rows(
  lapply(names(all_task_results), function(task_name) {
    all_task_results[[task_name]]$results$predictions
  })
)

all_selected_features <- bind_rows(
  lapply(names(all_task_results), function(task_name) {
    all_task_results[[task_name]]$results$selected_features
  })
)

all_tuning_table <- bind_rows(
  lapply(names(all_task_results), function(task_name) {
    all_task_results[[task_name]]$results$tuning_table
  })
)

write_csv(
  all_fold_metrics,
  file.path(table_out_dir, "elastic_net_fold_metrics.csv")
)

write_csv(
  all_predictions,
  file.path(table_out_dir, "elastic_net_outer_fold_predictions.csv")
)

write_csv(
  all_selected_features,
  file.path(table_out_dir, "elastic_net_selected_features_long.csv")
)

write_csv(
  all_tuning_table,
  file.path(table_out_dir, "elastic_net_inner_tuning_results.csv")
)

# ============================================================
# 9. Summary table with pooled AUC
# ============================================================

pooled_auc_table <- all_predictions %>%
  group_by(Feature_Set, Task, Method) %>%
  summarise(
    Pooled_AUC = safe_roc_auc_numeric(
      response = True_Binary,
      predictor = Predicted_Probability_Case
    ),
    .groups = "drop"
  )

summary_table <- all_fold_metrics %>%
  group_by(Feature_Set, Task, Method) %>%
  summarise(
    N_Control = first(N_Test_Control + N_Train_Control),
    N_Case = first(N_Test_Case + N_Train_Case),

    Fold_AUC_mean = mean(AUC, na.rm = TRUE),
    Fold_AUC_sd = sd(AUC, na.rm = TRUE),

    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_sd = sd(Accuracy, na.rm = TRUE),

    Balanced_Accuracy_mean = mean(Balanced_Accuracy, na.rm = TRUE),
    Balanced_Accuracy_sd = sd(Balanced_Accuracy, na.rm = TRUE),

    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Sensitivity_sd = sd(Sensitivity, na.rm = TRUE),

    Specificity_mean = mean(Specificity, na.rm = TRUE),
    Specificity_sd = sd(Specificity, na.rm = TRUE),

    PPV_mean = mean(PPV, na.rm = TRUE),
    PPV_sd = sd(PPV, na.rm = TRUE),

    NPV_mean = mean(NPV, na.rm = TRUE),
    NPV_sd = sd(NPV, na.rm = TRUE),

    F1_mean = mean(F1, na.rm = TRUE),
    F1_sd = sd(F1, na.rm = TRUE),

    Threshold_mean = mean(Threshold, na.rm = TRUE),
    Threshold_sd = sd(Threshold, na.rm = TRUE),

    Alpha_mean = mean(Best_Alpha, na.rm = TRUE),
    Lambda_mean = mean(Best_Lambda, na.rm = TRUE),

    N_Features_Filtered_mean = mean(N_Features_After_Filtering, na.rm = TRUE),
    N_Selected_Features_mean = mean(N_Selected_Features, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  left_join(
    pooled_auc_table,
    by = c("Feature_Set", "Task", "Method")
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  )

write_csv(
  summary_table,
  file.path(table_out_dir, "elastic_net_results_summary.csv")
)

print_summary <- summary_table %>%
  mutate(
    Fold_AUC = sprintf("%.3f ± %.3f", Fold_AUC_mean, Fold_AUC_sd),
    Balanced_Accuracy = sprintf("%.3f ± %.3f", Balanced_Accuracy_mean, Balanced_Accuracy_sd),
    Sensitivity = sprintf("%.3f ± %.3f", Sensitivity_mean, Sensitivity_sd),
    Specificity = sprintf("%.3f ± %.3f", Specificity_mean, Specificity_sd),
    F1 = sprintf("%.3f ± %.3f", F1_mean, F1_sd),
    Alpha = sprintf("%.2f", Alpha_mean),
    Lambda = sprintf("%.5f", Lambda_mean),
    Selected_Features = sprintf("%.1f", N_Selected_Features_mean)
  ) %>%
  select(
    Feature_Set,
    Task,
    Method,
    N_Control,
    N_Case,
    Pooled_AUC,
    Fold_AUC,
    Balanced_Accuracy,
    Sensitivity,
    Specificity,
    F1,
    Alpha,
    Lambda,
    Selected_Features
  )

write_csv(
  print_summary,
  file.path(table_out_dir, "elastic_net_results_summary_printable.csv")
)

cat("\n============================================================\n")
cat("ELASTIC NET RESULTS SUMMARY\n")
cat("============================================================\n\n")
print(print_summary, row.names = FALSE)

# ============================================================
# 10. ROC curves
# ============================================================

for (task_name in unique(all_predictions$Task)) {

  df_task <- all_predictions %>%
    filter(Task == task_name)

  if (nrow(df_task) == 0 ||
      length(unique(df_task$True_Binary)) < 2) {
    next
  }

  roc_obj <- safe_roc_obj_numeric(
    response = df_task$True_Binary,
    predictor = df_task$Predicted_Probability_Case
  )

  roc_plot_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Method = "Elastic Net"
  )

  pooled_auc <- as.numeric(roc_obj$auc)

  n_info <- summary_table %>%
    filter(Task == task_name) %>%
    slice(1)

  p_roc <- ggplot(
    roc_plot_df,
    aes(x = FPR, y = TPR, colour = Method)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.4
    ) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = method_color, drop = FALSE) +
    annotate(
      "text",
      x = 0.56,
      y = 0.24,
      label = paste0("Elastic Net: AUC = ", sprintf("%.3f", pooled_auc)),
      hjust = 0,
      size = 3.1,
      colour = "#1D3557",
      fontface = "italic"
    ) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = paste0(task_name, " - ", feature_set_name),
      subtitle = sprintf(
        "Elastic Net | n = %d controls / %d cases",
        n_info$N_Control[1],
        n_info$N_Case[1]
      ),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      colour = "Model"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, colour = "#1D3557"),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      axis.title = element_text(size = 10),
      panel.grid.major = element_line(colour = "grey92"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )

  ggsave(
    filename = file.path(
      figure_out_dir,
      paste0("ROC_", clean_filename(task_name), "_elastic_net.png")
    ),
    plot = p_roc,
    width = 7,
    height = 6,
    units = "in",
    dpi = 300
  )

  ggsave(
    filename = file.path(
      figure_out_dir,
      paste0("ROC_", clean_filename(task_name), "_elastic_net.pdf")
    ),
    plot = p_roc,
    width = 7,
    height = 6,
    units = "in",
    dpi = 300,
    device = "pdf"
  )
}

# ============================================================
# 11. Feature selection summary
# ============================================================

if (nrow(all_selected_features) > 0) {

  feature_frequency <- all_selected_features %>%
    group_by(Feature_Set, Task, Feature) %>%
    summarise(
      Times_Selected = n(),
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      SD_Coefficient = sd(Coefficient, na.rm = TRUE),
      Mean_Abs_Coefficient = mean(abs(Coefficient), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Task, desc(Times_Selected), desc(Mean_Abs_Coefficient))

} else {

  feature_frequency <- tibble(
    Feature_Set = character(),
    Task = character(),
    Feature = character(),
    Times_Selected = integer(),
    Mean_Coefficient = numeric(),
    SD_Coefficient = numeric(),
    Mean_Abs_Coefficient = numeric()
  )
}

write_csv(
  feature_frequency,
  file.path(table_out_dir, "elastic_net_feature_selection_frequency.csv")
)

# ============================================================
# 12. ROC diagnostics
# ============================================================

roc_diagnostics <- all_predictions %>%
  group_by(Feature_Set, Task, Method) %>%
  summarise(
    N = n(),
    N_Control = sum(True_Class == "control"),
    N_Case = sum(True_Class == "case"),
    Unique_Probabilities = n_distinct(round(Predicted_Probability_Case, 8)),
    Min_Probability = min(Predicted_Probability_Case, na.rm = TRUE),
    Max_Probability = max(Predicted_Probability_Case, na.rm = TRUE),
    Mean_Probability_Control = mean(Predicted_Probability_Case[True_Class == "control"], na.rm = TRUE),
    Mean_Probability_Case = mean(Predicted_Probability_Case[True_Class == "case"], na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  roc_diagnostics,
  file.path(table_out_dir, "elastic_net_ROC_diagnostics.csv")
)

# ============================================================
# 13. Session info
# ============================================================

sink(file.path(table_out_dir, "sessionInfo.txt"))
sessionInfo()
sink()

# ============================================================
# 14. Final message
# ============================================================

cat("\n============================================================\n")
cat("Elastic Net pipeline complete.\n")
cat("Feature set:", feature_set_name, "\n")
cat("Tables saved in:\n")
cat(table_out_dir, "\n")
cat("Figures saved in:\n")
cat(figure_out_dir, "\n")
cat("============================================================\n")
