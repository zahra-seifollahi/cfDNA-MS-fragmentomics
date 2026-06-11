# ============================================================
# Rebuild top_50_free and top_50_balanced from standardized
# classifier feature-selection summaries
#
# This is the GitHub-safe version of the original rebuilding code.
#
# It first creates:
#   feature_selection_main_summary_table.csv
#
# Then rebuilds:
#   04_top_50_free_features.csv
#   05_top_50_balanced_features.csv
#
# Inputs are the standardized outputs from:
#   02_run_svm_nested_cv.R
#   03_run_rf_nested_cv.R
#   04_run_elastic_net_nested_cv.R
#   05_compare_classifiers.R
#
# Usage:
#   Rscript scripts/07_classification/06_create_consensus_feature_lists.R
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

# ============================================================
# 0. Paths
# ============================================================

out_dir <- "results/tables/classification/consensus_features"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

performance_file <- "results/tables/classification/classifier_comparison/all_classifier_results_combined.csv"

feature_sets <- c(
  "global_endmotif_256",
  "regional_mds_167",
  "bioanalyzer_10bp_length"
)

feature_labels <- c(
  global_endmotif_256 = "End motif frequency",
  regional_mds_167 = "Regional MDS",
  bioanalyzer_10bp_length = "Bioanalyzer length bins"
)

# ============================================================
# 1. Helper functions
# ============================================================

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }

  read_csv(path, show_col_types = FALSE)
}

scale_0_1 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0

  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)

  if (!is.finite(mn) || !is.finite(mx) || mx == mn) {
    return(rep(0, length(x)))
  }

  (x - mn) / (mx - mn)
}

standardize_feature_set <- function(x) {
  x <- as.character(x)

  case_when(
    x %in% c(
      "global_endmotif_256",
      "End motif frequency",
      "Motif",
      "Motif frequency"
    ) ~ "End motif frequency",

    x %in% c(
      "regional_mds_167",
      "Regional MDS",
      "MDS",
      "Motif diversity score"
    ) ~ "Regional MDS",

    x %in% c(
      "bioanalyzer_10bp_length",
      "Bioanalyzer 10-bp",
      "Bioanalyzer length bins",
      "Fragment length-bin frequency"
    ) ~ "Bioanalyzer length bins",

    TRUE ~ x
  )
}

clean_feature_key <- function(feature, feature_set) {
  feature <- as.character(feature)
  feature_set <- standardize_feature_set(feature_set)

  feature <- str_replace(feature, "^End motif frequency::", "")
  feature <- str_replace(feature, "^Regional MDS::", "")
  feature <- str_replace(feature, "^Bioanalyzer 10-bp::", "")
  feature <- str_replace(feature, "^Bioanalyzer length bins::", "")
  feature <- str_trim(feature)

  if (feature_set == "End motif frequency") {
    motif_hit <- str_extract(toupper(feature), "[ACGT]{4}")
    return(ifelse(!is.na(motif_hit), motif_hit, toupper(feature)))
  }

  if (feature_set == "Regional MDS") {
    feature <- str_replace(feature, "^regional_MDS_", "")
    feature <- str_replace(feature, "^(chr[^:]+):(\\d+)-(\\d+)$", "\\1_\\2_\\3")
    feature <- str_replace_all(feature, "[:-]", "_")
    feature <- str_replace_all(feature, "__+", "_")
    return(tolower(feature))
  }

  if (feature_set == "Bioanalyzer length bins") {
    feature <- str_replace(feature, "^length_", "")
    feature <- str_replace(feature, "bp$", "")
    feature <- str_replace_all(feature, "Bioanalyzer_10bp_", "")
    feature <- str_replace_all(feature, "bioanalyzer_10bp_", "")
    feature <- str_replace_all(feature, "Bioanalyzer_", "")
    feature <- str_replace_all(feature, "bioanalyzer_", "")
    feature <- str_replace_all(feature, "bin_", "")
    feature <- str_replace_all(feature, "–", "-")
    feature <- str_replace_all(feature, "—", "-")
    feature <- str_replace_all(feature, "[^0-9A-Za-z]+", "_")
    feature <- str_replace_all(feature, "_+", "_")
    feature <- str_replace_all(feature, "^_|_$", "")
    return(tolower(feature))
  }

  tolower(feature)
}

make_display_name <- function(feature, feature_set) {
  feature <- as.character(feature)
  feature_set <- standardize_feature_set(feature_set)

  out <- feature

  idx_regional <- feature_set == "Regional MDS"
  if (any(idx_regional, na.rm = TRUE)) {
    x <- feature[idx_regional]
    x <- str_replace(x, "^regional_MDS_", "")
    out[idx_regional] <- ifelse(
      str_detect(x, "^chr[^_]+_\\d+_\\d+$"),
      str_replace(x, "^(chr[^_]+)_(\\d+)_(\\d+)$", "\\1:\\2-\\3"),
      x
    )
  }

  idx_bio <- feature_set == "Bioanalyzer length bins"
  if (any(idx_bio, na.rm = TRUE)) {
    x <- feature[idx_bio]
    x <- str_replace(x, "^length_", "")
    x <- str_replace(x, "bp$", "")
    x <- str_replace_all(x, "_", "-")
    out[idx_bio] <- paste0(x, " bp")
  }

  out
}

# ============================================================
# 2. Read Elastic Net selected features
# ============================================================

read_en_features <- function(feature_set_id) {
  feature_set_label <- feature_labels[[feature_set_id]]

  path <- file.path(
    "results/tables/classification/elastic_net",
    feature_set_id,
    "elastic_net_feature_selection_frequency.csv"
  )

  df <- safe_read_csv(path)

  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  if (!"Feature" %in% colnames(df)) {
    warning("Elastic Net feature file has no Feature column: ", path)
    return(NULL)
  }

  if (!"Times_Selected" %in% colnames(df)) {
    df$Times_Selected <- 1
  }

  if (!"Mean_Coefficient" %in% colnames(df)) {
    df$Mean_Coefficient <- NA_real_
  }

  if (!"Mean_Abs_Coefficient" %in% colnames(df)) {
    df$Mean_Abs_Coefficient <- abs(df$Mean_Coefficient)
  }

  df %>%
    mutate(
      Feature_Set = feature_set_label,
      Feature = as.character(Feature),
      Feature_Key = clean_feature_key(Feature, feature_set_label)
    ) %>%
    group_by(Task, Feature_Set, Feature_Key) %>%
    summarise(
      Feature = first(Feature),
      EN_Times_Selected = sum(Times_Selected, na.rm = TRUE),
      EN_Mean_Coefficient = mean(Mean_Coefficient, na.rm = TRUE),
      EN_Mean_Abs_Coefficient = mean(Mean_Abs_Coefficient, na.rm = TRUE),
      .groups = "drop"
    )
}

# ============================================================
# 3. Read SVM selected features
# ============================================================

read_svm_features <- function(feature_set_id) {
  feature_set_label <- feature_labels[[feature_set_id]]

  path <- file.path(
    "results/tables/classification/svm",
    feature_set_id,
    "svm_5methods_feature_selection_frequency.csv"
  )

  df <- safe_read_csv(path)

  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  if (!"Feature" %in% colnames(df)) {
    warning("SVM feature file has no Feature column: ", path)
    return(NULL)
  }

  if (!"Times_Selected" %in% colnames(df)) {
    df$Times_Selected <- 1
  }

  if (!"Method" %in% colnames(df)) {
    df$Method <- NA_character_
  }

  df %>%
    mutate(
      Feature_Set = feature_set_label,
      Feature = as.character(Feature),
      Feature_Key = clean_feature_key(Feature, feature_set_label)
    ) %>%
    group_by(Task, Feature_Set, Feature_Key) %>%
    summarise(
      Feature = first(Feature),
      SVM_Max_Times_Selected = max(Times_Selected, na.rm = TRUE),
      SVM_Total_Times_Selected = sum(Times_Selected, na.rm = TRUE),
      SVM_Methods_Selected = n_distinct(Method[Times_Selected > 0]),
      SVM_Best_Method = Method[which.max(Times_Selected)][1],
      .groups = "drop"
    )
}

# ============================================================
# 4. Read Random Forest importance
# ============================================================

read_rf_features <- function(feature_set_id) {
  feature_set_label <- feature_labels[[feature_set_id]]

  path <- file.path(
    "results/tables/classification/random_forest",
    feature_set_id,
    "rf_5methods_importance_summary.csv"
  )

  df <- safe_read_csv(path)

  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  if (!"Feature" %in% colnames(df)) {
    warning("RF importance file has no Feature column: ", path)
    return(NULL)
  }

  if (!"Mean_RF_Importance" %in% colnames(df)) {
    warning("RF importance file has no Mean_RF_Importance column: ", path)
    return(NULL)
  }

  if (!"Method" %in% colnames(df)) {
    df$Method <- NA_character_
  }

  df %>%
    mutate(
      Feature_Set = feature_set_label,
      Feature = as.character(Feature),
      Feature_Key = clean_feature_key(Feature, feature_set_label)
    ) %>%
    group_by(Task, Feature_Set, Feature_Key) %>%
    summarise(
      Feature = first(Feature),
      RF_MeanDecreaseAccuracy = mean(Mean_RF_Importance, na.rm = TRUE),
      RF_MaxDecreaseAccuracy = max(Mean_RF_Importance, na.rm = TRUE),
      RF_Methods_Selected = n_distinct(Method),
      RF_Best_Method = Method[which.max(Mean_RF_Importance)][1],
      .groups = "drop"
    )
}

# ============================================================
# 5. Build feature_selection_main_summary_table
# ============================================================

en_all <- bind_rows(map(feature_sets, read_en_features))
svm_all <- bind_rows(map(feature_sets, read_svm_features))
rf_all <- bind_rows(map(feature_sets, read_rf_features))

feature_summary <- full_join(
  en_all,
  svm_all,
  by = c("Task", "Feature_Set", "Feature_Key"),
  suffix = c("_EN", "_SVM")
) %>%
  full_join(
    rf_all,
    by = c("Task", "Feature_Set", "Feature_Key")
  ) %>%
  mutate(
    Feature = coalesce(Feature_EN, Feature_SVM, Feature),
    Feature_Set = standardize_feature_set(Feature_Set),

    EN_Times_Selected = replace_na(as.numeric(EN_Times_Selected), 0),
    SVM_Total_Times_Selected = replace_na(as.numeric(SVM_Total_Times_Selected), 0),
    SVM_Max_Times_Selected = replace_na(as.numeric(SVM_Max_Times_Selected), 0),
    SVM_Methods_Selected = replace_na(as.numeric(SVM_Methods_Selected), 0),

    RF_MeanDecreaseAccuracy = replace_na(as.numeric(RF_MeanDecreaseAccuracy), 0),
    RF_MaxDecreaseAccuracy = replace_na(as.numeric(RF_MaxDecreaseAccuracy), 0),
    RF_Methods_Selected = replace_na(as.numeric(RF_Methods_Selected), 0),

    Feature_Display = make_display_name(Feature, Feature_Set),

    EN_scaled = ave(
      EN_Times_Selected,
      Task,
      Feature_Set,
      FUN = scale_0_1
    ),

    SVM_scaled = ave(
      SVM_Total_Times_Selected,
      Task,
      Feature_Set,
      FUN = scale_0_1
    ),

    RF_scaled = ave(
      pmax(RF_MeanDecreaseAccuracy, 0),
      Task,
      Feature_Set,
      FUN = scale_0_1
    ),

    Consensus_Score = EN_scaled + SVM_scaled + RF_scaled
  ) %>%
  select(
    Task,
    Feature_Set,
    Feature,
    Feature_Key,
    Feature_Display,

    EN_Times_Selected,
    EN_Mean_Coefficient,
    EN_Mean_Abs_Coefficient,

    SVM_Max_Times_Selected,
    SVM_Total_Times_Selected,
    SVM_Methods_Selected,
    SVM_Best_Method,

    RF_MeanDecreaseAccuracy,
    RF_MaxDecreaseAccuracy,
    RF_Methods_Selected,
    RF_Best_Method,

    EN_scaled,
    SVM_scaled,
    RF_scaled,
    Consensus_Score
  ) %>%
  arrange(Task, Feature_Set, desc(Consensus_Score))

write_csv(
  feature_summary,
  file.path(out_dir, "feature_selection_main_summary_table.csv")
)

cat("\nFeature-selection summary table built:\n")
print(dim(feature_summary))
print(feature_summary %>% count(Feature_Set))

# ============================================================
# 6. Read model performance weights
# ============================================================

if (file.exists(performance_file)) {

  performance_weights <- read_csv(
    performance_file,
    show_col_types = FALSE
  ) %>%
    mutate(
      Task = as.character(Task),
      Feature_Set = standardize_feature_set(Feature_Set),
      Model = as.character(Model),
      Pooled_AUC = as.numeric(Pooled_AUC)
    ) %>%
    group_by(Task, Feature_Set, Model) %>%
    summarise(
      Performance_Weight = max(Pooled_AUC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Performance_Weight = ifelse(
        is.infinite(Performance_Weight) | is.nan(Performance_Weight),
        1,
        Performance_Weight
      ),
      Performance_Weight = pmin(pmax(Performance_Weight, 0), 1)
    )

} else {

  warning("Performance file not found. All model weights set to 1.")

  performance_weights <- expand.grid(
    Task = unique(feature_summary$Task),
    Feature_Set = unique(feature_summary$Feature_Set),
    Model = c("Elastic Net", "SVM", "Random Forest"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Performance_Weight = 1)
}

write_csv(
  performance_weights,
  file.path(out_dir, "performance_weights.csv")
)

# ============================================================
# 7. Build weighted evidence table
# ============================================================

en_evidence <- feature_summary %>%
  transmute(
    Task,
    Feature_Set,
    Feature,
    Feature_Display_Clean = Feature_Display,
    Model = "Elastic Net",
    Raw_Evidence = EN_Times_Selected
  )

svm_evidence <- feature_summary %>%
  transmute(
    Task,
    Feature_Set,
    Feature,
    Feature_Display_Clean = Feature_Display,
    Model = "SVM",
    Raw_Evidence = SVM_Total_Times_Selected
  )

rf_evidence <- feature_summary %>%
  transmute(
    Task,
    Feature_Set,
    Feature,
    Feature_Display_Clean = Feature_Display,
    Model = "Random Forest",
    Raw_Evidence = pmax(RF_MeanDecreaseAccuracy, 0)
  )

weighted_evidence <- bind_rows(
  en_evidence,
  svm_evidence,
  rf_evidence
) %>%
  mutate(
    Task = as.character(Task),
    Feature_Set = standardize_feature_set(Feature_Set),
    Feature = as.character(Feature),
    Raw_Evidence = replace_na(as.numeric(Raw_Evidence), 0),
    Raw_Evidence = pmax(Raw_Evidence, 0)
  ) %>%
  group_by(Task, Feature_Set, Model) %>%
  mutate(
    Evidence_Scaled = scale_0_1(Raw_Evidence)
  ) %>%
  ungroup() %>%
  left_join(
    performance_weights,
    by = c("Task", "Feature_Set", "Model")
  ) %>%
  mutate(
    Performance_Weight = replace_na(Performance_Weight, 1),
    Weighted_Evidence = Evidence_Scaled * Performance_Weight
  )

write_csv(
  weighted_evidence,
  file.path(out_dir, "weighted_feature_evidence_rebuilt.csv")
)

# ============================================================
# 8. Global feature ranking
# ============================================================

global_feature_ranking <- weighted_evidence %>%
  group_by(Feature_Set, Feature, Feature_Display_Clean) %>%
  summarise(
    Weighted_Consensus_Score = sum(Weighted_Evidence, na.rm = TRUE),
    Mean_Weighted_Consensus_Score = mean(Weighted_Evidence, na.rm = TRUE),

    Total_Raw_Evidence = sum(Raw_Evidence, na.rm = TRUE),

    N_Tasks_With_Evidence = n_distinct(Task[Raw_Evidence > 0]),
    N_Models_With_Evidence = n_distinct(Model[Raw_Evidence > 0]),
    N_Task_Model_Combinations = sum(Raw_Evidence > 0, na.rm = TRUE),

    EN_Total_Evidence = sum(Raw_Evidence[Model == "Elastic Net"], na.rm = TRUE),
    SVM_Total_Evidence = sum(Raw_Evidence[Model == "SVM"], na.rm = TRUE),
    RF_Total_Evidence = sum(Raw_Evidence[Model == "Random Forest"], na.rm = TRUE),

    Tasks_Selected = paste(sort(unique(Task[Raw_Evidence > 0])), collapse = "; "),
    Models_Selected = paste(sort(unique(Model[Raw_Evidence > 0])), collapse = "; "),

    .groups = "drop"
  ) %>%
  arrange(desc(Weighted_Consensus_Score)) %>%
  mutate(
    Global_Rank = row_number(),
    Feature_ID = paste(Feature_Set, Feature, sep = "::")
  )

write_csv(
  global_feature_ranking,
  file.path(out_dir, "global_feature_ranking_rebuilt.csv")
)

# ============================================================
# 9. Build Top50 free
# ============================================================

top_50_free <- global_feature_ranking %>%
  filter(N_Task_Model_Combinations >= 2) %>%
  arrange(desc(Weighted_Consensus_Score)) %>%
  slice_head(n = 50)

if (nrow(top_50_free) < 50) {
  top_50_free <- global_feature_ranking %>%
    arrange(desc(Weighted_Consensus_Score)) %>%
    slice_head(n = 50)
}

# ============================================================
# 10. Build Top50 balanced
# 20 motifs + 15 regional MDS + 15 Bioanalyzer bins
# ============================================================

top_50_balanced <- global_feature_ranking %>%
  filter(N_Task_Model_Combinations >= 1) %>%
  group_by(Feature_Set) %>%
  arrange(desc(Weighted_Consensus_Score), .by_group = TRUE) %>%
  mutate(Rank_In_Set = row_number()) %>%
  ungroup() %>%
  filter(
    (Feature_Set == "End motif frequency" & Rank_In_Set <= 20) |
      (Feature_Set == "Regional MDS" & Rank_In_Set <= 15) |
      (Feature_Set == "Bioanalyzer length bins" & Rank_In_Set <= 15)
  ) %>%
  arrange(desc(Weighted_Consensus_Score))

# If one feature set has fewer available features, fill up to 50
if (nrow(top_50_balanced) < 50) {

  already_selected <- top_50_balanced %>%
    mutate(id = paste(Feature_Set, Feature, sep = "||")) %>%
    pull(id)

  fill_features <- global_feature_ranking %>%
    mutate(id = paste(Feature_Set, Feature, sep = "||")) %>%
    filter(!id %in% already_selected) %>%
    arrange(desc(Weighted_Consensus_Score)) %>%
    slice_head(n = 50 - nrow(top_50_balanced)) %>%
    select(-id)

  top_50_balanced <- bind_rows(top_50_balanced, fill_features)
}

top_50_balanced <- top_50_balanced %>%
  arrange(desc(Weighted_Consensus_Score)) %>%
  slice_head(n = 50)

# ============================================================
# 11. Save with both original and standardized filenames
# ============================================================

write_csv(
  top_50_free,
  file.path(out_dir, "top_50_free_rebuilt.csv")
)

write_csv(
  top_50_balanced,
  file.path(out_dir, "top_50_balanced_rebuilt.csv")
)

# These are the filenames used by 07_run_consensus_top50_rf.R
write_csv(
  top_50_free,
  file.path(out_dir, "04_top_50_free_features.csv")
)

write_csv(
  top_50_balanced,
  file.path(out_dir, "05_top_50_balanced_features.csv")
)

cat("\nTop 50 free composition:\n")
print(top_50_free %>% count(Feature_Set))

cat("\nTop 50 balanced composition:\n")
print(top_50_balanced %>% count(Feature_Set))

cat("\n============================================================\n")
cat("Objects rebuilt: top_50_free and top_50_balanced\n")
cat("Outputs saved in:\n")
cat(out_dir, "\n")
cat("============================================================\n")
