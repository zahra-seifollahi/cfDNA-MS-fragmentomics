# ============================================================
# Consensus Top50 Random Forest using all feature sets
#
# Input feature sets:
#   results/intermediate/classification/feature_sets/global_endmotif_256.tsv
#   results/intermediate/classification/feature_sets/regional_mds_167.tsv
#   results/intermediate/classification/feature_sets/bioanalyzer_10bp_length.tsv
#
# Input consensus feature lists:
#   results/tables/classification/consensus_features/04_top_50_free_features.csv
#   results/tables/classification/consensus_features/05_top_50_balanced_features.csv
#
# Output:
#   results/tables/classification/consensus_top50_rf/
#   results/figures/classification/consensus_top50_rf/
#
# Important:
#   - This script does NOT rebuild raw matrices.
#   - It uses already prepared classification feature sets.
#   - It evaluates consensus-selected Top50 feature lists.
#   - Imputation is done inside each CV fold using training data only.
#
# Usage:
#   Rscript scripts/07_classification/06_run_consensus_top50_rf.R
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(caret)
  library(pROC)
  library(randomForest)
  library(ggplot2)
  library(forcats)
  library(scales)
})

set.seed(111)
options(stringsAsFactors = FALSE)

# ============================================================
# 0. Paths
# ============================================================

feature_dir <- "results/intermediate/classification/feature_sets"

global_endmotif_file <- file.path(
  feature_dir,
  "global_endmotif_256.tsv"
)

regional_mds_file <- file.path(
  feature_dir,
  "regional_mds_167.tsv"
)

bioanalyzer_file <- file.path(
  feature_dir,
  "bioanalyzer_10bp_length.tsv"
)

consensus_dir <- "results/tables/classification/consensus_features"

top50_free_file <- file.path(
  consensus_dir,
  "04_top_50_free_features.csv"
)

top50_balanced_file <- file.path(
  consensus_dir,
  "05_top_50_balanced_features.csv"
)

out_dir <- "results/tables/classification/consensus_top50_rf"
fig_dir <- "results/figures/classification/consensus_top50_rf"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Helper functions
# ============================================================

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

save_plot <- function(plot_object, filename, width, height, dpi = 300) {

  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = "pdf"
  )
}

normalize_feature_key <- function(feature, feature_set) {

  feature <- as.character(feature)
  feature_set <- as.character(feature_set)

  out <- feature

  out <- str_replace(out, "^End motif frequency::", "")
  out <- str_replace(out, "^Regional MDS::", "")
  out <- str_replace(out, "^Bioanalyzer 10-bp::", "")
  out <- str_replace(out, "^Bioanalyzer length bins::", "")
  out <- str_trim(out)

  result <- out

  idx_motif <- feature_set %in% c(
    "End motif frequency",
    "global_endmotif_256"
  )

  if (any(idx_motif, na.rm = TRUE)) {
    tmp <- toupper(out[idx_motif])
    motif_hit <- str_extract(tmp, "[ACGT]{4}")
    result[idx_motif] <- ifelse(!is.na(motif_hit), motif_hit, tmp)
  }

  idx_mds <- feature_set %in% c(
    "Regional MDS",
    "regional_mds_167"
  )

  if (any(idx_mds, na.rm = TRUE)) {
    tmp <- out[idx_mds]
    tmp <- str_replace(tmp, "^regional_MDS_", "")
    tmp <- str_replace(tmp, "^(chr[^:]+):(\\d+)-(\\d+)$", "\\1_\\2_\\3")
    tmp <- str_replace_all(tmp, "[:-]", "_")
    tmp <- str_replace_all(tmp, "__+", "_")
    tmp <- tolower(tmp)
    result[idx_mds] <- tmp
  }

  idx_bio <- feature_set %in% c(
    "Bioanalyzer 10-bp",
    "Bioanalyzer length bins",
    "Fragment length-bin frequency",
    "bioanalyzer_10bp_length"
  )

  if (any(idx_bio, na.rm = TRUE)) {
    tmp <- out[idx_bio]
    tmp <- str_replace_all(tmp, "^length_", "")
    tmp <- str_replace_all(tmp, "Bioanalyzer_10bp_", "")
    tmp <- str_replace_all(tmp, "bioanalyzer_10bp_", "")
    tmp <- str_replace_all(tmp, "Bioanalyzer_", "")
    tmp <- str_replace_all(tmp, "bioanalyzer_", "")
    tmp <- str_replace_all(tmp, "bin_", "")
    tmp <- str_replace_all(tmp, "bp", "")
    tmp <- str_replace_all(tmp, "–", "-")
    tmp <- str_replace_all(tmp, "—", "-")
    tmp <- str_replace_all(tmp, "[^0-9A-Za-z]+", "_")
    tmp <- str_replace_all(tmp, "_+", "_")
    tmp <- str_replace_all(tmp, "^_|_$", "")
    tmp <- tolower(tmp)
    result[idx_bio] <- tmp
  }

  idx_other <- !(idx_motif | idx_mds | idx_bio)

  if (any(idx_other, na.rm = TRUE)) {
    result[idx_other] <- tolower(result[idx_other])
  }

  result
}

standardize_feature_set_name <- function(x) {

  x <- as.character(x)

  case_when(
    x %in% c(
      "End motif frequency",
      "Motif",
      "Motif frequency",
      "global_endmotif_256"
    ) ~ "End motif frequency",

    x %in% c(
      "Regional MDS",
      "MDS",
      "Motif diversity score",
      "regional_mds_167"
    ) ~ "Regional MDS",

    x %in% c(
      "Bioanalyzer 10-bp",
      "Bioanalyzer length bins",
      "Fragment length-bin frequency",
      "bioanalyzer_10bp_length"
    ) ~ "Bioanalyzer length bins",

    TRUE ~ x
  )
}

format_region <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^regional_MDS_", "")
  ifelse(
    str_detect(x, "^chr[^_]+_\\d+_\\d+$"),
    str_replace(x, "^(chr[^_]+)_(\\d+)_(\\d+)$", "\\1:\\2-\\3"),
    x
  )
}

format_feature_label <- function(feature, feature_set) {

  feature <- as.character(feature)
  feature_set <- as.character(feature_set)

  case_when(
    feature_set == "Regional MDS" ~ format_region(feature),

    feature_set == "Bioanalyzer length bins" ~ feature %>%
      str_replace("^length_", "") %>%
      str_replace("bp$", " bp") %>%
      str_replace_all("_", "-"),

    TRUE ~ feature
  )
}

safe_auc <- function(y_true_binary, prob_case) {

  if (length(unique(y_true_binary[!is.na(y_true_binary)])) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = y_true_binary,
    predictor = prob_case,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(roc_obj$auc)
}

safe_roc_obj <- function(y_true_binary, prob_case) {

  pROC::roc(
    response = y_true_binary,
    predictor = prob_case,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

impute_train_test_median <- function(X_train, X_test) {

  medians <- apply(X_train, 2, function(x) median(x, na.rm = TRUE))
  medians[is.na(medians)] <- 0

  for (j in seq_along(medians)) {
    X_train[is.na(X_train[, j]), j] <- medians[j]
    X_test[is.na(X_test[, j]), j] <- medians[j]
  }

  list(
    X_train = X_train,
    X_test = X_test,
    medians = medians
  )
}

# ============================================================
# 2. Check input files
# ============================================================

check_file(global_endmotif_file)
check_file(regional_mds_file)
check_file(bioanalyzer_file)
check_file(top50_free_file)
check_file(top50_balanced_file)

# ============================================================
# 3. Read prepared feature sets
# ============================================================

read_feature_set <- function(path, feature_set_label) {

  df <- read_tsv(path, show_col_types = FALSE)

  if (!all(c("sample", "group") %in% colnames(df))) {
    stop("Feature file must contain sample and group columns: ", path)
  }

  df <- df %>%
    mutate(
      sample = as.character(sample),
      group = as.character(group)
    )

  feature_cols <- setdiff(colnames(df), c("sample", "group"))

  df_prefixed <- df %>%
    select(sample, group, all_of(feature_cols))

  prefixed_names <- paste(feature_set_label, feature_cols, sep = "::")

  colnames(df_prefixed)[match(feature_cols, colnames(df_prefixed))] <- prefixed_names

  df_prefixed
}

global_df <- read_feature_set(
  global_endmotif_file,
  "End motif frequency"
)

regional_df <- read_feature_set(
  regional_mds_file,
  "Regional MDS"
)

bio_df <- read_feature_set(
  bioanalyzer_file,
  "Bioanalyzer length bins"
)

common_samples <- Reduce(
  intersect,
  list(
    global_df$sample,
    regional_df$sample,
    bio_df$sample
  )
)

common_samples <- sort(common_samples)

if (length(common_samples) == 0) {
  stop("No common samples across the three feature sets.")
}

sample_info <- global_df %>%
  select(sample, group) %>%
  filter(sample %in% common_samples) %>%
  arrange(sample) %>%
  mutate(
    group = factor(group, levels = c("Healthy", "Remission", "Relapse"))
  )

combined_matrix <- sample_info %>%
  left_join(
    global_df %>% select(-group),
    by = "sample"
  ) %>%
  left_join(
    regional_df %>% select(-group),
    by = "sample"
  ) %>%
  left_join(
    bio_df %>% select(-group),
    by = "sample"
  ) %>%
  arrange(group, sample)

cat("\nCombined all-feature matrix:\n")
print(dim(combined_matrix))
print(table(combined_matrix$group))

write_csv(
  combined_matrix,
  file.path(out_dir, "combined_feature_matrix_all_feature_sets.csv")
)

write_csv(
  sample_info,
  file.path(out_dir, "sample_info_common_samples.csv")
)

# ============================================================
# 4. Feature map
# ============================================================

make_feature_map <- function(df, feature_set_label) {

  feature_cols <- setdiff(colnames(df), c("sample", "group"))

  tibble(
    Feature_Set = feature_set_label,
    Matrix_Column_Raw = str_replace(
      feature_cols,
      paste0("^", feature_set_label, "::"),
      ""
    ),
    Matrix_Column_Prefixed = feature_cols,
    Feature_Key = normalize_feature_key(
      str_replace(feature_cols, paste0("^", feature_set_label, "::"), ""),
      feature_set_label
    )
  )
}

feature_map_all <- bind_rows(
  make_feature_map(global_df, "End motif frequency"),
  make_feature_map(regional_df, "Regional MDS"),
  make_feature_map(bio_df, "Bioanalyzer length bins")
)

write_csv(
  feature_map_all,
  file.path(out_dir, "feature_name_mapping.csv")
)

# ============================================================
# 5. Read and map consensus Top50 feature lists
# ============================================================

read_consensus_feature_list <- function(path) {

  df <- read_csv(path, show_col_types = FALSE)

  if (!"Feature_Set" %in% colnames(df)) {
    stop("Consensus feature list must contain Feature_Set column: ", path)
  }

  if (!"Feature" %in% colnames(df)) {
    possible_feature_cols <- c(
      "Feature_ID",
      "Feature_Name",
      "Selected_Feature",
      "Matrix_Column_Raw",
      "Matrix_Column_Prefixed"
    )

    hit <- intersect(possible_feature_cols, colnames(df))

    if (length(hit) == 0) {
      stop("Consensus feature list must contain Feature column: ", path)
    }

    df <- df %>%
      rename(Feature = all_of(hit[1]))
  }

  if (!"Weighted_Consensus_Score" %in% colnames(df)) {
    df$Weighted_Consensus_Score <- NA_real_
  }

  df %>%
    mutate(
      Feature_Set = standardize_feature_set_name(Feature_Set),
      Feature = as.character(Feature),
      Feature_Key = normalize_feature_key(Feature, Feature_Set)
    )
}

map_feature_list_to_matrix <- function(feature_list, list_name) {

  mapped <- feature_list %>%
    left_join(
      feature_map_all,
      by = c("Feature_Set", "Feature_Key")
    )

  matched <- mapped %>%
    filter(!is.na(Matrix_Column_Prefixed))

  missing <- mapped %>%
    filter(is.na(Matrix_Column_Prefixed)) %>%
    select(
      Feature_Set,
      Feature,
      Feature_Key,
      Weighted_Consensus_Score
    )

  cat("\nFeature matching for", list_name, "\n")
  cat("Requested:", nrow(feature_list), "\n")
  cat("Matched:", nrow(matched), "\n")
  cat("Missing:", nrow(missing), "\n")
  print(matched %>% count(Feature_Set))

  if (nrow(missing) > 0) {
    write_csv(
      missing,
      file.path(out_dir, paste0("missing_features_", list_name, ".csv"))
    )
  }

  write_csv(
    matched,
    file.path(out_dir, paste0("matched_features_", list_name, ".csv"))
  )

  matched
}

top50_free <- read_consensus_feature_list(top50_free_file)
top50_balanced <- read_consensus_feature_list(top50_balanced_file)

top50_free_mapped <- map_feature_list_to_matrix(
  top50_free,
  "top50_free"
)

top50_balanced_mapped <- map_feature_list_to_matrix(
  top50_balanced,
  "top50_balanced"
)

# ============================================================
# 6. RF task settings
# ============================================================

tasks <- list(
  list(name = "Healthy vs MS",          case = c("Remission", "Relapse"), control = "Healthy"),
  list(name = "Healthy vs Remission",   case = "Remission",               control = "Healthy"),
  list(name = "Healthy vs Relapse",     case = "Relapse",                 control = "Healthy"),
  list(name = "NonRelapse vs Relapse",  case = "Relapse",                 control = c("Healthy", "Remission")),
  list(name = "Remission vs Relapse",   case = "Relapse",                 control = "Remission")
)

task_order <- c(
  "Healthy vs MS",
  "Healthy vs Relapse",
  "Healthy vs Remission",
  "NonRelapse vs Relapse",
  "Remission vs Relapse"
)

feature_colors <- c(
  "End motif frequency" = "#264653",
  "Regional MDS" = "#E9C46A",
  "Bioanalyzer length bins" = "#8ECAE6",
  "Unknown" = "grey70"
)

list_colors <- c(
  "Top50 free" = "#1D3557",
  "Top50 balanced" = "#E76F51"
)

model_colors <- c(
  "Elastic Net" = "#1D3557",
  "SVM" = "#E76F51",
  "Random Forest" = "#2A9D8F"
)

# ============================================================
# 7. Random Forest CV
# ============================================================

prepare_rf_input <- function(mapped_features, list_name) {

  feature_cols <- unique(mapped_features$Matrix_Column_Prefixed)

  if (length(feature_cols) < 2) {
    stop("Too few matched features for ", list_name, ": ", length(feature_cols))
  }

  rf_df <- combined_matrix %>%
    select(sample, group, all_of(feature_cols)) %>%
    filter(!is.na(group))

  cat("\nRF input for", list_name, "\n")
  cat("Samples:", nrow(rf_df), "\n")
  cat("Features:", length(feature_cols), "\n")
  print(table(rf_df$group))

  write_csv(
    rf_df,
    file.path(out_dir, paste0("RF_input_", list_name, ".csv"))
  )

  list(
    data = rf_df,
    feature_cols = feature_cols
  )
}

run_rf_for_task <- function(task,
                            rf_df,
                            feature_cols,
                            list_name,
                            n_folds = 5,
                            ntree = 1000) {

  task_df <- rf_df %>%
    filter(group %in% c(task$control, task$case)) %>%
    mutate(
      class = ifelse(group %in% task$case, "case", "control"),
      class = factor(class, levels = c("control", "case"))
    )

  cat("\nRunning task:", task$name, "|", list_name, "\n")
  print(table(task_df$class))

  if (nrow(task_df) < 10 || length(unique(task_df$class)) < 2) {
    warning("Skipping task: ", task$name)
    return(NULL)
  }

  folds <- caret::createFolds(
    task_df$class,
    k = n_folds,
    returnTrain = TRUE
  )

  pred_list <- list()
  metric_list <- list()
  importance_list <- list()

  for (i in seq_along(folds)) {

    train_idx <- folds[[i]]
    test_idx <- setdiff(seq_len(nrow(task_df)), train_idx)

    train_df <- task_df[train_idx, , drop = FALSE]
    test_df <- task_df[test_idx, , drop = FALSE]

    x_train_raw <- train_df[, feature_cols, drop = FALSE]
    x_test_raw <- test_df[, feature_cols, drop = FALSE]

    y_train <- train_df$class
    y_test <- test_df$class

    imp <- impute_train_test_median(
      X_train = x_train_raw,
      X_test = x_test_raw
    )

    x_train <- imp$X_train
    x_test <- imp$X_test

    mtry_value <- max(1, floor(sqrt(length(feature_cols))))

    rf_fit <- randomForest(
      x = x_train,
      y = y_train,
      ntree = ntree,
      mtry = mtry_value,
      importance = TRUE
    )

    prob_case <- predict(
      rf_fit,
      newdata = x_test,
      type = "prob"
    )[, "case"]

    pred_class <- factor(
      ifelse(prob_case >= 0.5, "case", "control"),
      levels = c("control", "case")
    )

    cm <- caret::confusionMatrix(
      data = pred_class,
      reference = y_test,
      positive = "case"
    )

    auc_val <- safe_auc(
      y_true_binary = ifelse(y_test == "case", 1, 0),
      prob_case = prob_case
    )

    pred_list[[i]] <- data.frame(
      Feature_List = list_name,
      Task = task$name,
      Fold = i,
      Sample = test_df$sample,
      True_Class = as.character(y_test),
      Predicted_Class = as.character(pred_class),
      Predicted_Probability_Case = prob_case,
      stringsAsFactors = FALSE
    )

    metric_list[[i]] <- data.frame(
      Feature_List = list_name,
      Task = task$name,
      Fold = i,
      N_Test = length(test_idx),
      N_Features = length(feature_cols),
      mtry = mtry_value,
      Accuracy = unname(cm$overall["Accuracy"]),
      Balanced_Accuracy = unname(cm$byClass["Balanced Accuracy"]),
      Sensitivity = unname(cm$byClass["Sensitivity"]),
      Specificity = unname(cm$byClass["Specificity"]),
      PPV = unname(cm$byClass["Pos Pred Value"]),
      NPV = unname(cm$byClass["Neg Pred Value"]),
      F1 = unname(cm$byClass["F1"]),
      AUC = auc_val,
      stringsAsFactors = FALSE
    )

    imp_rf <- importance(rf_fit)

    importance_list[[i]] <- data.frame(
      Feature_List = list_name,
      Task = task$name,
      Fold = i,
      Feature_ID = rownames(imp_rf),
      MeanDecreaseAccuracy = imp_rf[, "MeanDecreaseAccuracy"],
      MeanDecreaseGini = imp_rf[, "MeanDecreaseGini"],
      stringsAsFactors = FALSE
    )
  }

  list(
    predictions = bind_rows(pred_list),
    fold_metrics = bind_rows(metric_list),
    importance = bind_rows(importance_list)
  )
}

run_rf_feature_list <- function(mapped_features, list_name) {

  prepared <- prepare_rf_input(
    mapped_features = mapped_features,
    list_name = list_name
  )

  rf_df <- prepared$data
  feature_cols <- prepared$feature_cols

  rf_results <- lapply(
    tasks,
    run_rf_for_task,
    rf_df = rf_df,
    feature_cols = feature_cols,
    list_name = list_name,
    n_folds = 5,
    ntree = 1000
  )

  rf_results <- rf_results[!sapply(rf_results, is.null)]

  if (length(rf_results) == 0) {
    stop("No RF tasks completed for ", list_name)
  }

  all_predictions <- bind_rows(lapply(rf_results, `[[`, "predictions"))
  all_fold_metrics <- bind_rows(lapply(rf_results, `[[`, "fold_metrics"))
  all_importance <- bind_rows(lapply(rf_results, `[[`, "importance"))

  pooled_auc_table <- all_predictions %>%
    group_by(Feature_List, Task) %>%
    summarise(
      Pooled_AUC = safe_auc(
        y_true_binary = ifelse(True_Class == "case", 1, 0),
        prob_case = Predicted_Probability_Case
      ),
      .groups = "drop"
    )

  rf_summary <- all_fold_metrics %>%
    group_by(Feature_List, Task) %>%
    summarise(
      N_Features = first(N_Features),
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
      Fold_AUC_mean = mean(AUC, na.rm = TRUE),
      Fold_AUC_sd = sd(AUC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(pooled_auc_table, by = c("Feature_List", "Task")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  importance_summary <- all_importance %>%
    group_by(Feature_List, Feature_ID) %>%
    summarise(
      MeanDecreaseAccuracy_mean = mean(MeanDecreaseAccuracy, na.rm = TRUE),
      MeanDecreaseAccuracy_sd = sd(MeanDecreaseAccuracy, na.rm = TRUE),
      MeanDecreaseGini_mean = mean(MeanDecreaseGini, na.rm = TRUE),
      Times_Used = n(),
      .groups = "drop"
    ) %>%
    arrange(Feature_List, desc(MeanDecreaseAccuracy_mean))

  write_csv(
    all_predictions,
    file.path(out_dir, paste0("RF_predictions_", clean_filename(list_name), ".csv"))
  )

  write_csv(
    all_fold_metrics,
    file.path(out_dir, paste0("RF_fold_metrics_", clean_filename(list_name), ".csv"))
  )

  write_csv(
    rf_summary,
    file.path(out_dir, paste0("RF_results_summary_", clean_filename(list_name), ".csv"))
  )

  write_csv(
    all_importance,
    file.path(out_dir, paste0("RF_importance_long_", clean_filename(list_name), ".csv"))
  )

  write_csv(
    importance_summary,
    file.path(out_dir, paste0("RF_importance_summary_", clean_filename(list_name), ".csv"))
  )

  list(
    summary = rf_summary,
    predictions = all_predictions,
    fold_metrics = all_fold_metrics,
    importance_long = all_importance,
    importance_summary = importance_summary
  )
}

# ============================================================
# 8. Run RF for Top50 free and Top50 balanced
# ============================================================

rf_free <- run_rf_feature_list(
  mapped_features = top50_free_mapped,
  list_name = "Top50 free"
)

rf_balanced <- run_rf_feature_list(
  mapped_features = top50_balanced_mapped,
  list_name = "Top50 balanced"
)

# ============================================================
# 9. Compare Top50 free vs balanced
# ============================================================

rf_comparison <- bind_rows(
  rf_free$summary,
  rf_balanced$summary
) %>%
  mutate(
    Feature_List = factor(
      Feature_List,
      levels = c("Top50 free", "Top50 balanced")
    ),
    Task = factor(Task, levels = task_order)
  ) %>%
  select(
    Feature_List,
    Task,
    N_Features,
    Pooled_AUC,
    Fold_AUC_mean,
    Fold_AUC_sd,
    Balanced_Accuracy_mean,
    Balanced_Accuracy_sd,
    Sensitivity_mean,
    Specificity_mean,
    F1_mean
  )

write_csv(
  rf_comparison,
  file.path(out_dir, "RF_free_vs_balanced_summary.csv")
)

print(rf_comparison)

p_rf_compare <- ggplot(
  rf_comparison,
  aes(
    x = Task,
    y = Pooled_AUC,
    fill = Feature_List
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Pooled_AUC)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.2
  ) +
  scale_fill_manual(values = list_colors) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Random Forest performance using consensus-selected features",
    subtitle = "Comparison of free and balanced top-50 feature lists",
    x = "Classification task",
    y = "Pooled AUC",
    fill = "Feature list"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

save_plot(
  p_rf_compare,
  "01_RF_top50_free_vs_balanced_pooled_AUC",
  width = 9,
  height = 5.5
)

# ============================================================
# 10. Annotate RF importance by feature set
# ============================================================

importance_long <- bind_rows(
  rf_free$importance_long,
  rf_balanced$importance_long
) %>%
  mutate(
    Feature_Set = case_when(
      str_detect(Feature_ID, "^End motif frequency::") ~ "End motif frequency",
      str_detect(Feature_ID, "^Regional MDS::") ~ "Regional MDS",
      str_detect(Feature_ID, "^Bioanalyzer length bins::") ~ "Bioanalyzer length bins",
      TRUE ~ "Unknown"
    ),
    Feature = Feature_ID %>%
      str_replace("^End motif frequency::", "") %>%
      str_replace("^Regional MDS::", "") %>%
      str_replace("^Bioanalyzer length bins::", ""),
    Feature_Display = format_feature_label(Feature, Feature_Set),
    MeanDecreaseAccuracy_pos = pmax(MeanDecreaseAccuracy, 0),
    Task = factor(Task, levels = task_order),
    Feature_List = factor(
      Feature_List,
      levels = c("Top50 free", "Top50 balanced")
    )
  )

importance_summary <- importance_long %>%
  group_by(Feature_List, Feature_ID, Feature_Set, Feature, Feature_Display) %>%
  summarise(
    MeanDecreaseAccuracy_mean = mean(MeanDecreaseAccuracy, na.rm = TRUE),
    MeanDecreaseAccuracy_pos_mean = mean(MeanDecreaseAccuracy_pos, na.rm = TRUE),
    MeanDecreaseAccuracy_sd = sd(MeanDecreaseAccuracy, na.rm = TRUE),
    Positive_Importance_Count = sum(MeanDecreaseAccuracy > 0, na.rm = TRUE),
    Total_Count = n(),
    Positive_Importance_Fraction = Positive_Importance_Count / Total_Count,
    .groups = "drop"
  ) %>%
  arrange(Feature_List, desc(MeanDecreaseAccuracy_mean))

write_csv(
  importance_long,
  file.path(out_dir, "RF_importance_long_combined.csv")
)

write_csv(
  importance_summary,
  file.path(out_dir, "RF_importance_summary_for_plots.csv")
)

# ============================================================
# 11. Top important features
# ============================================================

plot_top_features_bar <- function(df, feature_list_name, top_n = 25) {

  plot_df <- df %>%
    filter(Feature_List == feature_list_name) %>%
    slice_max(MeanDecreaseAccuracy_mean, n = top_n, with_ties = FALSE) %>%
    mutate(
      Feature_Display = ifelse(
        nchar(Feature_Display) > 55,
        paste0(substr(Feature_Display, 1, 52), "..."),
        Feature_Display
      ),
      Feature_Display = factor(
        Feature_Display,
        levels = rev(Feature_Display)
      )
    )

  ggplot(
    plot_df,
    aes(
      x = Feature_Display,
      y = MeanDecreaseAccuracy_mean,
      fill = Feature_Set
    )
  ) +
    geom_col(color = "black", linewidth = 0.25) +
    coord_flip() +
    scale_fill_manual(values = feature_colors, drop = FALSE) +
    labs(
      title = paste0("Top important features in ", feature_list_name, " Random Forest"),
      subtitle = "Features ranked by mean decrease accuracy across tasks and folds",
      x = "Feature",
      y = "Mean decrease accuracy",
      fill = "Feature set"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(color = "black"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

p_top_free <- plot_top_features_bar(
  importance_summary,
  "Top50 free",
  top_n = 25
)

save_plot(
  p_top_free,
  "02_top25_important_features_top50_free",
  width = 9,
  height = 8
)

p_top_balanced <- plot_top_features_bar(
  importance_summary,
  "Top50 balanced",
  top_n = 25
)

save_plot(
  p_top_balanced,
  "03_top25_important_features_top50_balanced",
  width = 9,
  height = 8
)

# ============================================================
# 12. Importance distribution by feature set
# ============================================================

p_importance_by_set <- importance_long %>%
  filter(Feature_Set != "Unknown") %>%
  ggplot(
    aes(
      x = Feature_Set,
      y = MeanDecreaseAccuracy,
      fill = Feature_Set
    )
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.55,
    color = "black",
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.05,
    size = 0.5,
    alpha = 0.3,
    color = "grey25"
  ) +
  facet_wrap(~ Feature_List, nrow = 1) +
  scale_fill_manual(values = feature_colors, drop = FALSE) +
  labs(
    title = "Distribution of Random Forest feature importance by feature set",
    subtitle = "Mean decrease accuracy values across tasks and folds",
    x = "Feature set",
    y = "Mean decrease accuracy",
    fill = "Feature set"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_plot(
  p_importance_by_set,
  "04_importance_distribution_by_feature_set",
  width = 10,
  height = 5.5
)

# ============================================================
# 13. Importance heatmaps
# ============================================================

make_importance_heatmap <- function(feature_list_name, top_n = 50) {

  top_features <- importance_summary %>%
    filter(Feature_List == feature_list_name) %>%
    slice_max(MeanDecreaseAccuracy_mean, n = top_n, with_ties = FALSE) %>%
    pull(Feature_ID)

  heat_df <- importance_long %>%
    filter(
      Feature_List == feature_list_name,
      Feature_ID %in% top_features
    ) %>%
    group_by(Feature_List, Task, Feature_ID, Feature_Set, Feature_Display) %>%
    summarise(
      Task_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Feature_Display = factor(
        Feature_Display,
        levels = rev(
          importance_summary %>%
            filter(
              Feature_List == feature_list_name,
              Feature_ID %in% top_features
            ) %>%
            arrange(desc(MeanDecreaseAccuracy_mean)) %>%
            pull(Feature_Display)
        )
      )
    )

  ggplot(
    heat_df,
    aes(
      x = Task,
      y = Feature_Display,
      fill = Task_MeanDecreaseAccuracy
    )
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "white",
      high = "#1a9850",
      midpoint = 0,
      name = "Mean decrease\naccuracy"
    ) +
    labs(
      title = paste0("Feature importance heatmap: ", feature_list_name),
      subtitle = "Top features ranked by overall Random Forest importance",
      x = "Classification task",
      y = "Feature"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
      panel.grid = element_blank()
    )
}

p_heat_free <- make_importance_heatmap(
  "Top50 free",
  top_n = 50
)

save_plot(
  p_heat_free,
  "05_importance_heatmap_top50_free",
  width = 5,
  height = 9
)

p_heat_balanced <- make_importance_heatmap(
  "Top50 balanced",
  top_n = 50
)

save_plot(
  p_heat_balanced,
  "06_importance_heatmap_top50_balanced",
  width = 5,
  height = 9
)

# ============================================================
# 14. Stability vs importance scatter
# ============================================================

p_stability <- importance_summary %>%
  filter(Feature_Set != "Unknown") %>%
  ggplot(
    aes(
      x = Positive_Importance_Fraction,
      y = MeanDecreaseAccuracy_mean,
      color = Feature_Set,
      shape = Feature_List
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.4
  ) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(values = feature_colors, drop = FALSE) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Feature importance stability in final Random Forest models",
    subtitle = "Stable features have both high importance and frequent positive contribution",
    x = "Fraction of task-fold combinations with positive importance",
    y = "Mean decrease accuracy",
    color = "Feature set",
    shape = "Feature list"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

save_plot(
  p_stability,
  "07_feature_importance_stability_scatter",
  width = 8,
  height = 6
)

# ============================================================
# 15. Final candidate feature tables
# ============================================================

final_candidate_features_free <- importance_summary %>%
  filter(Feature_List == "Top50 free") %>%
  arrange(desc(MeanDecreaseAccuracy_mean))

final_candidate_features_balanced <- importance_summary %>%
  filter(Feature_List == "Top50 balanced") %>%
  arrange(desc(MeanDecreaseAccuracy_mean))

write_csv(
  final_candidate_features_free,
  file.path(out_dir, "final_candidate_features_ranked_top50_free.csv")
)

write_csv(
  final_candidate_features_balanced,
  file.path(out_dir, "final_candidate_features_ranked_top50_balanced.csv")
)

# ============================================================
# 16. Session info and final message
# ============================================================

sink(file.path(out_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\n============================================================\n")
cat("Consensus Top50 Random Forest completed.\n")
cat("Tables saved in:\n")
cat(out_dir, "\n")
cat("Figures saved in:\n")
cat(fig_dir, "\n")
cat("============================================================\n")
