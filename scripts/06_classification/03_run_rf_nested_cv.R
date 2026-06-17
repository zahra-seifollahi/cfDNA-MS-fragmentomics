# ============================================================
# Random Forest nested CV classifier for cfDNA-MS feature sets
#
# Preserved from original workflow:
#   - Same 5 binary tasks
#   - Outer 5-fold CV
#   - Inner 3-fold CV
#   - Feature selection inside CV
#   - mtry tuning inside inner CV
#   - Threshold selection from inner CV
#   - Same 5 feature-selection methods:
#       1) RFE
#       2) RF importance
#       3) Limma FDR < 0.05
#       4) Top50 Limma
#       5) Top50 Jonckheere-Terpstra trend
#
# Usage:
#   Rscript scripts/07_classification/03_run_rf_nested_cv.R bioanalyzer_10bp_length
#   Rscript scripts/07_classification/03_run_rf_nested_cv.R regional_mds_167
#   Rscript scripts/07_classification/03_run_rf_nested_cv.R global_endmotif_256
# ============================================================

required_packages <- c(
  "caret",
  "randomForest",
  "limma",
  "pROC",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "stringr",
  "clinfun",
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
  library(caret)
  library(randomForest)
  library(limma)
  library(pROC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(clinfun)
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
    "Rscript scripts/07_classification/03_run_rf_nested_cv.R bioanalyzer_10bp_length"
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
  "results/tables/classification/random_forest",
  feature_set_name
)

figure_out_dir <- file.path(
  "results/figures/classification/random_forest",
  feature_set_name
)

dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("Random Forest nested CV classification\n")
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
# 2. Tasks and method labels
# ============================================================

tasks <- list(
  list(name = "Healthy vs MS",          case = c("Remission", "Relapse"), control = "Healthy"),
  list(name = "Healthy vs Remission",   case = "Remission",               control = "Healthy"),
  list(name = "Healthy vs Relapse",     case = "Relapse",                 control = "Healthy"),
  list(name = "NonRelapse vs Relapse",  case = "Relapse",                 control = c("Healthy", "Remission")),
  list(name = "Remission vs Relapse",   case = "Relapse",                 control = "Remission")
)

method_order <- c(
  "rfe",
  "rf_importance",
  "limma_sig",
  "top50_limma",
  "top50_trend"
)

method_labels <- c(
  rfe = "RFE",
  rf_importance = "RF importance",
  limma_sig = "Limma FDR<0.05",
  top50_limma = "Top50 Limma",
  top50_trend = "Top50 trend"
)

method_colors <- c(
  "RFE" = "#E63946",
  "RF importance" = "#F77F00",
  "Limma FDR<0.05" = "#E9C46A",
  "Top50 Limma" = "#2A9D8F",
  "Top50 trend" = "#457B9D"
)

# ============================================================
# 3. Helper functions
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

safe_roc_auc <- function(response, predictor) {

  response <- factor(response, levels = c("control", "case"))

  if (length(unique(response[!is.na(response)])) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c("control", "case"),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(roc_obj$auc)
}

safe_roc_obj <- function(response, predictor) {

  response <- factor(response, levels = c("control", "case"))

  pROC::roc(
    response = response,
    predictor = predictor,
    levels = c("control", "case"),
    direction = "<",
    quiet = TRUE
  )
}

train_rf_fixed_mtry <- function(X_train,
                                y_train,
                                mtry_value,
                                ntree = 1000) {

  y_train <- factor(y_train, levels = c("control", "case"))

  randomForest(
    x = X_train,
    y = y_train,
    ntree = ntree,
    mtry = mtry_value,
    importance = TRUE
  )
}

# ============================================================
# 4. Feature selection functions
# ============================================================

select_rfe <- function(X_train, y_train, n_inner_folds = 3) {

  y_train <- factor(y_train, levels = c("control", "case"))

  sizes <- c(10, 20, 50, 100)
  sizes <- sizes[sizes <= ncol(X_train)]

  if (length(sizes) == 0) {
    sizes <- min(10, ncol(X_train))
  }

  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "cv",
    number = n_inner_folds,
    allowParallel = FALSE
  )

  rfe_result <- tryCatch(
    rfe(
      x = X_train,
      y = y_train,
      sizes = sizes,
      rfeControl = ctrl,
      metric = "Accuracy"
    ),
    error = function(e) NULL
  )

  if (is.null(rfe_result)) {
    return(character(0))
  }

  predictors(rfe_result)
}

select_rf_importance <- function(X_train, y_train) {

  y_train <- factor(y_train, levels = c("control", "case"))

  rf_model <- tryCatch(
    randomForest(
      x = X_train,
      y = y_train,
      ntree = 500,
      importance = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(rf_model)) {
    return(character(0))
  }

  importance_scores <- importance(rf_model, type = 1)

  importance_df <- data.frame(
    Feature = rownames(importance_scores),
    Importance = importance_scores[, 1],
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Importance))

  selected <- importance_df$Feature[importance_df$Importance > 0]

  if (length(selected) < 2) {
    selected <- importance_df$Feature[1:min(50, nrow(importance_df))]
  }

  selected
}

select_limma_sig <- function(X_train, y_train, adj_cutoff = 0.05) {

  y_lm <- factor(y_train, levels = c("control", "case"))

  design <- model.matrix(~ 0 + y_lm)
  colnames(design) <- c("control", "case")

  fit <- lmFit(t(X_train), design)

  cont.matrix <- makeContrasts(
    case_vs_control = case - control,
    levels = design
  )

  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  tt <- topTable(
    fit2,
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
  )

  rownames(tt)[tt$adj.P.Val < adj_cutoff]
}

select_top50_limma <- function(X_train, y_train, n_top = 50) {

  y_lm <- factor(y_train, levels = c("control", "case"))

  design <- model.matrix(~ 0 + y_lm)
  colnames(design) <- c("control", "case")

  fit <- lmFit(t(X_train), design)

  cont.matrix <- makeContrasts(
    case_vs_control = case - control,
    levels = design
  )

  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  tt <- topTable(
    fit2,
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
  )

  rownames(tt)[1:min(n_top, nrow(tt))]
}

select_top50_trend <- function(X_train, group_train_original, n_top = 50) {

  group_train_original <- factor(
    as.character(group_train_original),
    levels = c("Healthy", "Remission", "Relapse"),
    ordered = TRUE
  )

  g_num <- as.numeric(group_train_original)

  trend_results <- lapply(colnames(X_train), function(m) {

    x <- X_train[, m]

    p_inc <- tryCatch(
      clinfun::jonckheere.test(x, g_num, alternative = "increasing")$p.value,
      error = function(e) NA_real_
    )

    p_dec <- tryCatch(
      clinfun::jonckheere.test(x, g_num, alternative = "decreasing")$p.value,
      error = function(e) NA_real_
    )

    data.frame(
      Feature = m,
      p_increasing = p_inc,
      p_decreasing = p_dec,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  trend_results <- trend_results %>%
    mutate(
      padj_increasing = p.adjust(p_increasing, method = "BH"),
      padj_decreasing = p.adjust(p_decreasing, method = "BH"),
      best_padj = pmin(padj_increasing, padj_decreasing, na.rm = TRUE),
      direction = ifelse(
        padj_increasing <= padj_decreasing,
        "Increasing",
        "Decreasing"
      )
    ) %>%
    arrange(best_padj)

  trend_results$Feature[1:min(n_top, nrow(trend_results))]
}

select_features <- function(method_code,
                            X_train,
                            y_train_binary,
                            group_train_original,
                            n_inner_folds_for_rfe = 3,
                            n_top = 50) {

  if (method_code == "rfe") {
    return(select_rfe(
      X_train,
      y_train_binary,
      n_inner_folds = n_inner_folds_for_rfe
    ))
  }

  if (method_code == "rf_importance") {
    return(select_rf_importance(
      X_train,
      y_train_binary
    ))
  }

  if (method_code == "limma_sig") {
    return(select_limma_sig(
      X_train,
      y_train_binary,
      adj_cutoff = 0.05
    ))
  }

  if (method_code == "top50_limma") {
    return(select_top50_limma(
      X_train,
      y_train_binary,
      n_top = n_top
    ))
  }

  if (method_code == "top50_trend") {
    return(select_top50_trend(
      X_train,
      group_train_original,
      n_top = n_top
    ))
  }

  stop("Unknown method_code: ", method_code)
}

# ============================================================
# 5. Inner CV: feature selection + mtry tuning
# ============================================================

tune_mtry_inner_cv <- function(X_outer_train_raw,
                               y_outer_train,
                               group_outer_train_original,
                               method_code,
                               n_inner_folds = 3,
                               n_top = 50,
                               ntree_inner = 500) {

  y_outer_train <- factor(y_outer_train, levels = c("control", "case"))

  inner_folds <- createFolds(
    y_outer_train,
    k = n_inner_folds,
    returnTrain = TRUE
  )

  prediction_rows <- list()

  for (j in seq_along(inner_folds)) {

    inner_train_idx <- inner_folds[[j]]
    inner_valid_idx <- setdiff(seq_len(nrow(X_outer_train_raw)), inner_train_idx)

    X_inner_train_raw <- X_outer_train_raw[inner_train_idx, , drop = FALSE]
    X_inner_valid_raw <- X_outer_train_raw[inner_valid_idx, , drop = FALSE]

    y_inner_train <- y_outer_train[inner_train_idx]
    y_inner_valid <- y_outer_train[inner_valid_idx]

    g_inner_train <- group_outer_train_original[inner_train_idx]

    imp_inner <- impute_train_test_median(
      X_inner_train_raw,
      X_inner_valid_raw
    )

    X_inner_train <- imp_inner$X_train
    X_inner_valid <- imp_inner$X_test

    selected_vars <- select_features(
      method_code = method_code,
      X_train = X_inner_train,
      y_train_binary = y_inner_train,
      group_train_original = g_inner_train,
      n_inner_folds_for_rfe = n_inner_folds,
      n_top = n_top
    )

    selected_vars <- intersect(selected_vars, colnames(X_inner_train))

    if (length(selected_vars) < 2) {
      next
    }

    X_train_fs <- X_inner_train[, selected_vars, drop = FALSE]
    X_valid_fs <- X_inner_valid[, selected_vars, drop = FALSE]

    candidate_mtry <- unique(round(c(
      1,
      sqrt(ncol(X_train_fs)),
      ncol(X_train_fs) / 4,
      ncol(X_train_fs) / 2
    )))

    candidate_mtry <- candidate_mtry[candidate_mtry >= 1]
    candidate_mtry <- candidate_mtry[candidate_mtry <= ncol(X_train_fs)]
    candidate_mtry <- sort(unique(candidate_mtry))

    for (mtry_value in candidate_mtry) {

      rf_model <- tryCatch(
        train_rf_fixed_mtry(
          X_train = X_train_fs,
          y_train = y_inner_train,
          mtry_value = mtry_value,
          ntree = ntree_inner
        ),
        error = function(e) NULL
      )

      if (is.null(rf_model)) {
        next
      }

      pred_prob <- predict(
        rf_model,
        X_valid_fs,
        type = "prob"
      )[, "case"]

      prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
        Inner_Fold = j,
        mtry = mtry_value,
        True_Class = as.character(y_inner_valid),
        Predicted_Probability_Case = as.numeric(pred_prob),
        N_Selected_Features = length(selected_vars),
        stringsAsFactors = FALSE
      )
    }
  }

  pred_all <- bind_rows(prediction_rows)

  if (nrow(pred_all) == 0) {
    return(list(
      best_mtry = NA_real_,
      best_threshold = NA_real_,
      tuning_table = data.frame(),
      inner_predictions_best = data.frame()
    ))
  }

  tuning_table <- pred_all %>%
    group_by(mtry) %>%
    summarise(
      Inner_AUC = safe_roc_auc(
        response = True_Class,
        predictor = Predicted_Probability_Case
      ),
      Mean_N_Selected_Features = mean(N_Selected_Features, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Inner_AUC))

  if (all(is.na(tuning_table$Inner_AUC))) {
    return(list(
      best_mtry = NA_real_,
      best_threshold = NA_real_,
      tuning_table = tuning_table,
      inner_predictions_best = data.frame()
    ))
  }

  best_mtry <- tuning_table$mtry[1]

  pred_best <- pred_all %>%
    filter(mtry == best_mtry)

  inner_roc <- safe_roc_obj(
    response = pred_best$True_Class,
    predictor = pred_best$Predicted_Probability_Case
  )

  best_threshold <- coords(
    inner_roc,
    x = "best",
    best.method = "youden",
    transpose = FALSE
  )$threshold[1]

  list(
    best_mtry = best_mtry,
    best_threshold = as.numeric(best_threshold),
    tuning_table = tuning_table,
    inner_predictions_best = pred_best
  )
}

# ============================================================
# 6. Outer CV
# ============================================================

run_nested_cv_rf_5methods <- function(data,
                                      y_binary,
                                      group_original,
                                      task_name,
                                      n_outer_folds = 5,
                                      n_inner_folds = 3,
                                      n_top = 50,
                                      ntree_inner = 500,
                                      ntree_outer = 1000) {

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
  rf_importance_rows <- list()

  for (i in seq_along(outer_folds)) {

    cat(sprintf("      Outer fold %d/%d\n", i, length(outer_folds)))

    outer_train_idx <- outer_folds[[i]]
    outer_test_idx  <- setdiff(seq_len(nrow(data)), outer_train_idx)

    X_outer_train_raw0 <- data[outer_train_idx, , drop = FALSE]
    X_outer_test_raw0  <- data[outer_test_idx,  , drop = FALSE]

    y_outer_train <- y_binary[outer_train_idx]
    y_outer_test  <- y_binary[outer_test_idx]

    g_outer_train <- group_original[outer_train_idx]

    filt <- filter_train_test_features(
      X_train = X_outer_train_raw0,
      X_test = X_outer_test_raw0,
      max_missing_fraction = 0.10,
      min_sd = 1e-8
    )

    X_outer_train_raw <- filt$X_train
    X_outer_test_raw  <- filt$X_test

    for (method_code in method_order) {

      method_label <- method_labels[[method_code]]

      cat(sprintf("        Method: %s\n", method_label))

      tune_res <- tune_mtry_inner_cv(
        X_outer_train_raw = X_outer_train_raw,
        y_outer_train = y_outer_train,
        group_outer_train_original = g_outer_train,
        method_code = method_code,
        n_inner_folds = n_inner_folds,
        n_top = n_top,
        ntree_inner = ntree_inner
      )

      best_mtry <- tune_res$best_mtry
      best_threshold <- tune_res$best_threshold

      tuning_table <- tune_res$tuning_table %>%
        mutate(
          Feature_Set = feature_set_name,
          Task = task_name,
          Fold = i,
          Method = method_label,
          Method_Code = method_code
        )

      tuning_rows_all[[length(tuning_rows_all) + 1]] <- tuning_table

      if (is.na(best_mtry)) {

        fold_metrics[[length(fold_metrics) + 1]] <- data.frame(
          Feature_Set = feature_set_name,
          Task = task_name,
          Fold = i,
          Method = method_label,
          Method_Code = method_code,
          N_Train = length(outer_train_idx),
          N_Test = length(outer_test_idx),
          N_Train_Control = sum(y_outer_train == "control"),
          N_Train_Case = sum(y_outer_train == "case"),
          N_Test_Control = sum(y_outer_test == "control"),
          N_Test_Case = sum(y_outer_test == "case"),
          N_Features_After_Filtering = ncol(X_outer_train_raw),
          N_Selected_Features = 0,
          Best_mtry = NA_real_,
          Threshold = NA_real_,
          Accuracy = NA_real_,
          Balanced_Accuracy = NA_real_,
          Sensitivity = NA_real_,
          Specificity = NA_real_,
          PPV = NA_real_,
          NPV = NA_real_,
          F1 = NA_real_,
          AUC = NA_real_,
          stringsAsFactors = FALSE
        )

        next
      }

      imp_outer <- impute_train_test_median(
        X_train = X_outer_train_raw,
        X_test = X_outer_test_raw
      )

      X_outer_train <- imp_outer$X_train
      X_outer_test  <- imp_outer$X_test

      selected_vars <- select_features(
        method_code = method_code,
        X_train = X_outer_train,
        y_train_binary = y_outer_train,
        group_train_original = g_outer_train,
        n_inner_folds_for_rfe = n_inner_folds,
        n_top = n_top
      )

      selected_vars <- intersect(selected_vars, colnames(X_outer_train))

      if (length(selected_vars) < 2) {

        fold_metrics[[length(fold_metrics) + 1]] <- data.frame(
          Feature_Set = feature_set_name,
          Task = task_name,
          Fold = i,
          Method = method_label,
          Method_Code = method_code,
          N_Train = length(outer_train_idx),
          N_Test = length(outer_test_idx),
          N_Train_Control = sum(y_outer_train == "control"),
          N_Train_Case = sum(y_outer_train == "case"),
          N_Test_Control = sum(y_outer_test == "control"),
          N_Test_Case = sum(y_outer_test == "case"),
          N_Features_After_Filtering = ncol(X_outer_train),
          N_Selected_Features = length(selected_vars),
          Best_mtry = best_mtry,
          Threshold = best_threshold,
          Accuracy = NA_real_,
          Balanced_Accuracy = NA_real_,
          Sensitivity = NA_real_,
          Specificity = NA_real_,
          PPV = NA_real_,
          NPV = NA_real_,
          F1 = NA_real_,
          AUC = NA_real_,
          stringsAsFactors = FALSE
        )

        next
      }

      X_train_fs <- X_outer_train[, selected_vars, drop = FALSE]
      X_test_fs  <- X_outer_test[,  selected_vars, drop = FALSE]

      best_mtry <- min(best_mtry, ncol(X_train_fs))
      best_mtry <- max(1, best_mtry)

      rf_model <- tryCatch(
        train_rf_fixed_mtry(
          X_train = X_train_fs,
          y_train = y_outer_train,
          mtry_value = best_mtry,
          ntree = ntree_outer
        ),
        error = function(e) NULL
      )

      if (is.null(rf_model)) {
        next
      }

      pred_prob <- predict(
        rf_model,
        X_test_fs,
        type = "prob"
      )[, "case"]

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

      outer_auc <- safe_roc_auc(
        response = as.character(y_outer_test),
        predictor = pred_prob
      )

      fold_metrics[[length(fold_metrics) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Task = task_name,
        Fold = i,
        Method = method_label,
        Method_Code = method_code,
        N_Train = length(outer_train_idx),
        N_Test = length(outer_test_idx),
        N_Train_Control = sum(y_outer_train == "control"),
        N_Train_Case = sum(y_outer_train == "case"),
        N_Test_Control = sum(y_outer_test == "control"),
        N_Test_Case = sum(y_outer_test == "case"),
        N_Features_After_Filtering = ncol(X_outer_train),
        N_Selected_Features = length(selected_vars),
        Best_mtry = best_mtry,
        Threshold = best_threshold,
        Accuracy = unname(cm$overall["Accuracy"]),
        Balanced_Accuracy = unname(cm$byClass["Balanced Accuracy"]),
        Sensitivity = unname(cm$byClass["Sensitivity"]),
        Specificity = unname(cm$byClass["Specificity"]),
        PPV = unname(cm$byClass["Pos Pred Value"]),
        NPV = unname(cm$byClass["Neg Pred Value"]),
        F1 = unname(cm$byClass["F1"]),
        AUC = outer_auc,
        stringsAsFactors = FALSE
      )

      prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Classifier = "Random Forest",
        Task = task_name,
        Fold = i,
        Method = method_label,
        Method_Code = method_code,
        Sample = rownames(data)[outer_test_idx],
        True_Class = as.character(y_outer_test),
        Predicted_Probability_Case = as.numeric(pred_prob),
        Predicted_Class = as.character(pred_class),
        stringsAsFactors = FALSE
      )

      selected_feature_rows[[length(selected_feature_rows) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Classifier = "Random Forest",
        Task = task_name,
        Fold = i,
        Method = method_label,
        Method_Code = method_code,
        Feature = selected_vars,
        stringsAsFactors = FALSE
      )

      rf_imp <- importance(rf_model, type = 1)

      rf_importance_rows[[length(rf_importance_rows) + 1]] <- data.frame(
        Feature_Set = feature_set_name,
        Classifier = "Random Forest",
        Task = task_name,
        Fold = i,
        Method = method_label,
        Method_Code = method_code,
        Feature = rownames(rf_imp),
        RF_MeanDecreaseAccuracy = rf_imp[, 1],
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    fold_metrics = bind_rows(fold_metrics),
    predictions = bind_rows(prediction_rows),
    selected_features = bind_rows(selected_feature_rows),
    tuning_table = bind_rows(tuning_rows_all),
    rf_importance = bind_rows(rf_importance_rows)
  )
}

# ============================================================
# 7. Run all tasks
# ============================================================

cat("\n============================================================\n")
cat("RUNNING RANDOM FOREST NESTED CV\n")
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

  task_y_original <- factor(
    as.character(groups[keep_idx]),
    levels = c("Healthy", "Remission", "Relapse"),
    ordered = TRUE
  )

  rownames(task_data) <- rownames(feature_matrix)[keep_idx]

  res <- run_nested_cv_rf_5methods(
    data = task_data,
    y_binary = task_y_binary,
    group_original = task_y_original,
    task_name = task$name,
    n_outer_folds = 5,
    n_inner_folds = 3,
    n_top = 50,
    ntree_inner = 500,
    ntree_outer = 1000
  )

  all_task_results[[task$name]] <- list(
    n_control = length(control_idx),
    n_case = length(case_idx),
    results = res
  )

  cat("\n")
}

# ============================================================
# 8. Combine and save outputs
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

all_rf_importance <- bind_rows(
  lapply(names(all_task_results), function(task_name) {
    all_task_results[[task_name]]$results$rf_importance
  })
)

write_csv(
  all_fold_metrics,
  file.path(table_out_dir, "rf_5methods_fold_metrics.csv")
)

write_csv(
  all_predictions,
  file.path(table_out_dir, "rf_5methods_outer_fold_predictions.csv")
)

write_csv(
  all_selected_features,
  file.path(table_out_dir, "rf_5methods_selected_features_long.csv")
)

write_csv(
  all_tuning_table,
  file.path(table_out_dir, "rf_5methods_inner_tuning_results.csv")
)

write_csv(
  all_rf_importance,
  file.path(table_out_dir, "rf_5methods_importance_long.csv")
)

# ============================================================
# 9. Summary table
# ============================================================

pooled_auc_table <- all_predictions %>%
  group_by(Feature_Set, Task, Method, Method_Code) %>%
  summarise(
    Pooled_AUC = safe_roc_auc(
      response = True_Class,
      predictor = Predicted_Probability_Case
    ),
    .groups = "drop"
  )

summary_table <- all_fold_metrics %>%
  group_by(Feature_Set, Task, Method, Method_Code) %>%
  summarise(
    Classifier = "Random Forest",

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

    Best_mtry_mean = mean(Best_mtry, na.rm = TRUE),

    N_Features_Filtered_mean = mean(N_Features_After_Filtering, na.rm = TRUE),
    N_Selected_Features_mean = mean(N_Selected_Features, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  left_join(
    pooled_auc_table,
    by = c("Feature_Set", "Task", "Method", "Method_Code")
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  ) %>%
  arrange(Task, factor(Method_Code, levels = method_order))

write_csv(
  summary_table,
  file.path(table_out_dir, "rf_5methods_results_summary.csv")
)

print_summary <- summary_table %>%
  mutate(
    Fold_AUC = sprintf("%.3f ± %.3f", Fold_AUC_mean, Fold_AUC_sd),
    Balanced_Accuracy = sprintf("%.3f ± %.3f", Balanced_Accuracy_mean, Balanced_Accuracy_sd),
    Sensitivity = sprintf("%.3f ± %.3f", Sensitivity_mean, Sensitivity_sd),
    Specificity = sprintf("%.3f ± %.3f", Specificity_mean, Specificity_sd),
    F1 = sprintf("%.3f ± %.3f", F1_mean, F1_sd),
    Selected_Features = sprintf("%.1f", N_Selected_Features_mean)
  ) %>%
  select(
    Feature_Set,
    Classifier,
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
    Selected_Features,
    Best_mtry_mean
  )

write_csv(
  print_summary,
  file.path(table_out_dir, "rf_5methods_results_summary_printable.csv")
)

cat("\n============================================================\n")
cat("RANDOM FOREST RESULTS SUMMARY\n")
cat("============================================================\n\n")
print(print_summary, row.names = FALSE)

# ============================================================
# 10. ROC curves
# ============================================================

make_roc_df_for_task <- function(task_name) {

  df_task <- all_predictions %>%
    filter(Task == task_name)

  roc_list <- lapply(unique(df_task$Method), function(m) {

    df_m <- df_task %>%
      filter(Method == m)

    if (nrow(df_m) == 0 || length(unique(df_m$True_Class)) < 2) {
      return(NULL)
    }

    roc_obj <- safe_roc_obj(
      response = df_m$True_Class,
      predictor = df_m$Predicted_Probability_Case
    )

    data.frame(
      Feature_Set = feature_set_name,
      Task = task_name,
      Method = m,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Pooled_AUC = as.numeric(roc_obj$auc),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(roc_list)
}

for (task_name in unique(all_predictions$Task)) {

  roc_plot_df <- make_roc_df_for_task(task_name)

  if (nrow(roc_plot_df) == 0) {
    next
  }

  auc_labels <- roc_plot_df %>%
    group_by(Method) %>%
    summarise(
      Pooled_AUC = first(Pooled_AUC),
      .groups = "drop"
    ) %>%
    arrange(factor(Method, levels = method_labels[method_order])) %>%
    mutate(
      label = paste0(Method, ": AUC = ", sprintf("%.3f", Pooled_AUC))
    )

  label_text <- paste(auc_labels$label, collapse = "\n")

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
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = method_colors) +
    annotate(
      "text",
      x = 0.56,
      y = 0.24,
      label = label_text,
      hjust = 0,
      size = 3.1,
      colour = "#1D3557",
      fontface = "italic"
    ) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = paste0(task_name, " — ", feature_set_name),
      subtitle = sprintf(
        "Random Forest | n = %d controls / %d cases",
        n_info$N_Control[1],
        n_info$N_Case[1]
      ),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      colour = "Method"
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
      paste0("ROC_", clean_filename(task_name), "_RF_5methods.png")
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
      paste0("ROC_", clean_filename(task_name), "_RF_5methods.pdf")
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
# 11. Feature selection and importance summaries
# ============================================================

top_feature_frequency <- all_selected_features %>%
  group_by(Feature_Set, Task, Method, Method_Code, Feature) %>%
  summarise(
    Times_Selected = n(),
    .groups = "drop"
  ) %>%
  arrange(Task, Method, desc(Times_Selected))

write_csv(
  top_feature_frequency,
  file.path(table_out_dir, "rf_5methods_feature_selection_frequency.csv")
)

rf_importance_summary <- all_rf_importance %>%
  group_by(Feature_Set, Task, Method, Method_Code, Feature) %>%
  summarise(
    Mean_RF_Importance = mean(RF_MeanDecreaseAccuracy, na.rm = TRUE),
    SD_RF_Importance = sd(RF_MeanDecreaseAccuracy, na.rm = TRUE),
    Times_In_Model = n(),
    .groups = "drop"
  ) %>%
  arrange(Task, Method, desc(Mean_RF_Importance))

write_csv(
  rf_importance_summary,
  file.path(table_out_dir, "rf_5methods_importance_summary.csv")
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
  file.path(table_out_dir, "rf_5methods_ROC_diagnostics.csv")
)

# ============================================================
# 13. Final message
# ============================================================

cat("\n============================================================\n")
cat("Random Forest pipeline complete.\n")
cat("Feature set:", feature_set_name, "\n")
cat("Tables saved in:\n")
cat(table_out_dir, "\n")
cat("Figures saved in:\n")
cat(figure_out_dir, "\n")
cat("============================================================\n")
