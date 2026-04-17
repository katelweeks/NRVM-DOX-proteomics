suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pcaMethods)
  library(tidyr)
  library(dplyr)
  library(tibble)
})

# ---------------- Configuration ----------------
set.seed(24)

condition_file <- "Data/Experimental_Annotation.csv"
protein_file   <- "Data/combined_protein_n22.tsv"
output_dir     <- "Results/"

fdr_thr   <- 0.05
lfc_thr   <- log2(1.3) # log2 fold-change threshold
use_imputation    <- FALSE
imputation_nPcs   <- 4
imputation_center <- TRUE
min_prop_non_na   <- 1.0  # keep proteins quantified in at least this proportion of samples

save_plot <- function(path, plot, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, ...)
}

# ---------------- Load metadata ----------------
meta <- read.delim(condition_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("FragPipe_ID", "Batch", "Bioreplicate", "Treatment", "Sex", "Condition")
missing_cols  <- setdiff(required_cols, colnames(meta))
if (length(missing_cols)) stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))

meta <- meta %>%
  mutate(
    Sample = paste0(FragPipe_ID, "_MaxLFQ_Intensity"),
    Treatment = factor(Treatment),
    Batch = factor(Batch),
    Bioreplicate = factor(Bioreplicate),
    Sex = factor(Sex)
  )

# ---------------- Load protein intensities ----------------
prot <- read.delim(protein_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

sample_cols <- meta$Sample
missing_samples <- setdiff(sample_cols, colnames(prot))
if (length(missing_samples)) stop("Sample columns missing from protein file: ", paste(missing_samples, collapse = ", "))

exprs <- as.matrix(prot[, sample_cols, drop = FALSE])
rownames(exprs) <- prot[["Protein"]]
suppressWarnings(mode(exprs) <- "numeric")
exprs[exprs == 0] <- NA_real_

# ---------------- Colours & Ordering ----------------
colors <- c("VehicleMale" = "#81D4FA",    # Light blue
            "DoxMale" = "#1976D2",        # Vibrant blue
            "VehicleFemale" = "#EF9A9A",  # Light red
            "DoxFemale" = "#D32F2F")      # Vibrant red

condition_order <- c("VehicleMale", "DoxMale", "VehicleFemale", "DoxFemale")

# ---------------- QC: missingness ----------------
non_missing <- colSums(!is.na(exprs))
missing_counts <- colSums(is.na(exprs))

non_missing_df <- tibble(Sample = names(non_missing), NonMissing = as.numeric(non_missing))
missing_df     <- tibble(Sample = names(missing_counts), Missing = as.numeric(missing_counts))

plot_missingness <- function(df, title) {
  df %>%
    left_join(meta, by = "Sample") %>%
    mutate(Sample_short = gsub("_MaxLFQ_Intensity", "", Sample),
           Sample_short = factor(Sample_short, levels = gtools::mixedsort(unique(Sample_short))),
           Condition = factor(Condition, levels = condition_order)) %>%
    ggplot(aes(x = Sample_short, y = Value, fill = Condition)) +
    geom_col(alpha = 0.7, width = 0.6) +
    scale_x_discrete(limits = rev) +
    scale_fill_manual(values = colors) +
    facet_grid(Condition ~ ., scales = "free_y", space = "free") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Count") +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 20),
      panel.spacing.y = unit(0.5, "lines")
    )
}

non_missing_long <- tibble(Sample = names(non_missing), Value = as.numeric(non_missing)) 
missing_long     <- tibble(Sample = names(missing_counts), Value = as.numeric(missing_counts))

save_plot(file.path(output_dir, "non_missing_values.png"),
          plot_missingness(non_missing_long, "Non-missing values per sample"),
          width = 5, height = 8, dpi = 300, bg = "transparent")

save_plot(file.path(output_dir, "missing_values.png"),
          plot_missingness(missing_long, "Missing values per sample"),
          width = 5, height = 8, dpi = 300, bg = "transparent")

# ---------------- Filter proteins ----------------
min_non_na <- ceiling(ncol(exprs) * min_prop_non_na)
exprs <- exprs[rowSums(!is.na(exprs)) >= min_non_na, , drop = FALSE]
if (!nrow(exprs)) stop("No proteins remain after missingness filtering.")

# ---------------- Impute + normalize ----------------
if (use_imputation) {
  message("--- Using PCA/SVD imputation for missing values ---")
  exprs_log <- log2(exprs + 1)
  all_na_rows <- rowSums(is.na(exprs_log)) == ncol(exprs_log)
  all_na_cols <- colSums(is.na(exprs_log)) == nrow(exprs_log)
  exprs_log <- exprs_log[!all_na_rows, !all_na_cols, drop = FALSE]
  
  pc_fit <- pcaMethods::pca(
    exprs_log,
    method = "svdImpute",
    nPcs = imputation_nPcs,
    center = imputation_center,
    scale = "none"
  )
  exprs_imputed <- pc_fit@completeObs
} else {
  message("--- Using complete-case filtering (no imputation) ---")
  exprs_imputed <- log2(exprs + 1)
  exprs_imputed <- exprs_imputed[complete.cases(exprs_imputed), , drop = FALSE]
  if (!nrow(exprs_imputed)) stop("No proteins remain after complete-case filtering.")
}

exprs_imputed <- exprs_imputed[, sample_cols, drop = FALSE]
exprs_norm <- limma::normalizeBetweenArrays(exprs_imputed, method = "quantile")
message(sprintf("Final expression matrix: %d proteins x %d samples", nrow(exprs_norm), ncol(exprs_norm)))


# ---------------- PCA and QC ----------------
qc_plots<-function(exprs_norm_df,sample_info,outfile_prefix_path){
  #exprs_norm_df = exprs_norm
  #sample_columns = sample_cols
  #sample_info = meta
  #outfile_prefix_path = paste0(output_dir,"n22_")
  exprs_pca <- exprs_norm_df
  pca <- prcomp(t(exprs_pca), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  
  pca_df<-pca$x%>%data.frame()%>%rownames_to_column(var="Sample")
  
  pca_df=left_join(pca_df,sample_info,by="Sample")
  
  batch_colors <- c("1" = "#4CAF50", # Darker Pastel Green
                    "2" = "#1976D2") # Darker Pastel Blue
  
  # remove '_MaxLFQ_Intensity' from sample labels for plotting
  pca_df$Sample <- gsub("_MaxLFQ_Intensity", "", pca_df$Sample)
  
  # PCA coloured by Condition
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition)) +
    geom_point(size = 4) +
    scale_color_manual(values = colors) +
    ggrepel::geom_text_repel(aes(label = Sample), size = 4) +
    labs(
      title = "PCA (samples)",
      x = sprintf("PC1 (%.1f%%)", var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", var_expl[2])
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
  save_plot(paste0(outfile_prefix_path,"_pca_by_condition.png"), p_pca, width = 8, height = 6, dpi = 300, bg = "transparent")
  
  # PCA coloured by Batch
  p_pca_batch <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Batch)) +
    geom_point(size = 4) +
    scale_color_manual(values = batch_colors) +
    ggrepel::geom_text_repel(aes(label = Sample), size = 4) +
    labs(
      title = "PCA (samples) coloured by Batch",
      x = sprintf("PC1 (%.1f%%)", var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", var_expl[2])
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
  save_plot(paste0(outfile_prefix_path,"_pca_by_batch.png"), p_pca_batch, width = 7, height = 6, dpi = 300, bg = "transparent")
  
  # Cumulative variance plot
  var_df <- tibble(
    PCs = seq_along(pca$sdev),
    CumulativeVariance = cumsum(pca$sdev^2) / sum(pca$sdev^2)
  )
  p_var <- ggplot(var_df, aes(x = PCs, y = CumulativeVariance)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq_along(var_df$PCs)) +
    labs(
      title = "Cumulative variance explained",
      x = "Number of PCs",
      y = "Cumulative variance explained"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
  save_plot(paste0(outfile_prefix_path, "_variance_explained.png"), p_var, width = 7, height = 6, dpi = 300, bg = "transparent")
  
  # ---------------- %CV calculation ----------------
  message("Calculating %CV values (back-transforming from log2)...")
  
  # Back-transform to linear scale
  exprs_linear <- 2 ^ exprs_norm_df
  
  # Helper function to compute %CV on linear values
  calc_cv <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2 || mean(x) <= 0) return(NA_real_)
    sd(x) / mean(x) * 100
  }
  
  # --- Overall CV per protein ---
  cv_overall <- apply(exprs_linear, 1, calc_cv)
  
  # --- CV by Condition (Vehicle, Dox, etc.) ---
  conditions <- unique(sample_info$Condition)
  
  cv_by_condition <- lapply(conditions, function(cond) {
    samples <- sample_info$Sample[sample_info$Condition == cond]
    mat_sub <- exprs_linear[, samples, drop = FALSE]
    apply(mat_sub, 1, calc_cv)
  })
  names(cv_by_condition) <- conditions
  
  # Combine into a table
  cv_df <- tibble(Protein = rownames(exprs_linear), CV_all = cv_overall)
  
  for (cond in conditions) {
    cv_df[[paste0("CV_", cond)]] <- cv_by_condition[[cond]]
  }
  
  # Save table
  write.csv(cv_df, paste0(outfile_prefix_path,"_protein_CV_by_condition.csv"), row.names = FALSE)
  
  # ---------------- Violin by sex-condition ----------------
  cv_long_sex <- cv_df %>%
    select(Protein, CV_VehicleMale, CV_DoxMale, CV_VehicleFemale, CV_DoxFemale) %>%
    pivot_longer(cols = starts_with("CV_"), names_to = "Group", values_to = "CV") %>%
    filter(!is.na(CV)) %>%
    mutate(Group = case_when(
      Group == "CV_VehicleMale" ~ "VehMale",
      Group == "CV_DoxMale"     ~ "DoxMale",
      Group == "CV_VehicleFemale" ~ "VehFemale", 
      Group == "CV_DoxFemale"   ~ "DoxFemale"
    ),
    Group = factor(Group, levels = c("VehMale", "DoxMale", "VehFemale", "DoxFemale")))
  
  p_cv_violin_sex <- ggplot(cv_long_sex, aes(x = Group, y = CV, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.75, colour = "black") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = c("VehMale" = "#A1D4FA",
                                 "DoxMale" = "#1976D2",        
                                 "VehFemale" = "#EF9A9A",  
                                 "DoxFemale" = "#D32F2F")) +   
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)) +
    labs(title = "%CV by Sex-Treatment (violin plot)",
         x = "", y = "%CV") +
    ylim(0, 100)
  
  
  save_plot(paste0(outfile_prefix_path,"_CV_violin_sex_treatment.png"), p_cv_violin_sex,
            width = 8.7, height = 6, dpi = 300, bg = "transparent")
  
  
}

qc_plots(exprs_norm_df = exprs_norm,sample_info = meta,outfile_prefix_path = paste0(output_dir,"n22"))



#Samples “Dox_5” and “Vehicle_9” had the highest number of missing values. 
#Therefore, these samples were excluded from the analysis.
exclude_samples   <- c("Dox_5_MaxLFQ_Intensity", "Vehicle_9_MaxLFQ_Intensity") # samples to exclude from analysis

if (length(exclude_samples)) {
  
  meta_n20 <- filter(meta, !Sample %in% exclude_samples)
  exprs_n20 <- as.matrix(prot[, meta_n20$Sample, drop = FALSE])
  rownames(exprs_n20) <- prot[["Protein"]]
  exprs_n20[exprs_n20 == 0] <- NA_real_
  suppressWarnings(mode(exprs_n20) <- "numeric")
  #remove NA values
  exprs_n20<-na.omit(exprs_n20)
  exprs_n20_ltr <- log2(exprs_n20 + 1)
  exprs_n20_ltr_norm <- limma::normalizeBetweenArrays(exprs_n20_ltr, method = "quantile")
  qc_plots(exprs_norm_df = exprs_n20_ltr_norm,sample_info = meta_n20,outfile_prefix_path = paste0(output_dir,"n20"))
}

# ---------------- Export normalized + imputed intensities ----------------
exprs_export <- as.data.frame(exprs_n20_ltr_norm)
exprs_export$Protein <- rownames(exprs_n20_ltr_norm)
exprs_export <- exprs_export %>% select(Protein, everything()) # Move Protein to first column

# Add gene annotations if available
if ("Gene" %in% colnames(prot)) {
  exprs_export <- exprs_export %>%
    left_join(prot %>% select(Protein, Gene) %>% distinct(), by = "Protein")
}

# Safer column naming with explicit matching
sample_data_cols <- setdiff(colnames(exprs_export), c("Protein", "Gene"))
for(i in seq_along(sample_data_cols)) {
  col_name <- sample_data_cols[i]
  matching_row <- which(meta$Sample == col_name)
  if(length(matching_row) == 1) {
    new_name <- paste0(col_name, " (", meta$Condition[matching_row], ")")
    colnames(exprs_export)[colnames(exprs_export) == col_name] <- new_name
  }
}

write.csv(exprs_export, file.path(output_dir, "n20_normalised_imputed_intensities.csv"),
          row.names = FALSE, na = "")

message("Exported normalised + imputed intensities to: ", file.path(output_dir, "n20_normalised_imputed_intensities.csv"))

# ---------------- QC: boxplots ----------------
plot_box <- function(mat, title) {
  as.data.frame(mat) %>%
    rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Sample", values_to = "Abundance") %>%
    left_join(meta, by = "Sample") %>%
    mutate(Sample_short = gsub("_MaxLFQ_Intensity", "", Sample),
           Sample_short = factor(Sample_short, levels = gtools::mixedsort(unique(Sample_short))),
           Condition = factor(Condition, levels = condition_order)) %>%
    ggplot(aes(x = Sample_short, y = Abundance, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    scale_x_discrete(limits = rev) +  # REVERSED: 1,2,3 top-to-bottom
    scale_fill_manual(values = colors) +
    facet_grid(Condition ~ ., scales = "free_y", space = "free") +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(log[2]("abundance"))) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "top",
      plot.margin = margin(10, 10, 10, 20),
      panel.spacing.y = unit(0.5, "lines")
    )
}

save_plot(file.path(output_dir, "n20_boxplot_log2_imputed.png"),
          plot_box(exprs_n20_ltr, "Protein abundance (log2)"),
          width = 6, height = 8, dpi = 300, bg = "transparent")

save_plot(file.path(output_dir, "n20_boxplot_log2_norm.png"),
          plot_box(exprs_n20_ltr_norm, "Protein abundance (log2 + quantile norm)"),
          width = 6, height = 8, dpi = 300, bg = "transparent")

# ---------------- Design 1: main effects only ----------------
expected_treat <- c("Vehicle", "Dox")
if (!all(expected_treat %in% meta_n20$Treatment)) {
  stop("Treatment column must contain: ", paste(expected_treat, collapse = ", "))
}

meta_n20$Treatment <- relevel(meta_n20$Treatment, ref = "Vehicle")
meta_n20$Sex       <- relevel(meta_n20$Sex, ref = "Female")

prot_20<-prot%>%select(-exclude_samples)


design_main <- model.matrix(~ 0 + Treatment + Sex + Batch, data = meta_n20)
colnames(design_main)[seq_along(levels(meta_n20$Treatment))] <- levels(meta_n20$Treatment)

corfit_main <- duplicateCorrelation(exprs_n20_ltr_norm,
                                    design = design_main,
                                    block  = meta_n20$Bioreplicate)

fit_main <- lmFit(exprs_n20_ltr_norm,
                  design = design_main,
                  block  = meta_n20$Bioreplicate,
                  correlation = corfit_main$consensus)

cont.matrix.main <- makeContrasts(
  DoxVsVehicle   = Dox - Vehicle,
  MaleVsFemale   = SexMale,
  Batch2VsBatch1 = Batch2,
  levels = design_main
)

fit_main2 <- contrasts.fit(fit_main, cont.matrix.main)
fit_main2 <- eBayes(fit_main2, trend = TRUE)

# ---------------- Volcano helper ----------------
make_volcano <- function(fit_obj, coef_name, annotation_df, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  tt <- topTable(fit_obj, coef = coef_name, number = Inf) %>%
    rownames_to_column("Protein") %>%
    left_join(annotation_df %>% select(Protein, Gene), by = "Protein") %>%
    mutate(
      neglog10FDR = -log10(adj.P.Val),
      category = case_when(
        adj.P.Val < fdr_cutoff & logFC >= lfc_cutoff  ~ "Up",
        adj.P.Val < fdr_cutoff & logFC <= -lfc_cutoff ~ "Down",
        TRUE ~ "NotSig"
      ),
      rank_sig = ifelse(category != "NotSig", rank(-neglog10FDR, ties.method = "min"), Inf),
      rank_fc  = ifelse(category != "NotSig", rank(-abs(logFC), ties.method = "min"), Inf),
      label_me = (rank_sig <= 25 | rank_fc <= 25) & category != "NotSig",
      point_labels = Gene
    )
  
  p <- ggplot(tt, aes(x = logFC, y = neglog10FDR, color = category)) +
    geom_point(aes(alpha = category != "NotSig"), size = 3) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 3) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(tt, label_me),
      aes(label = point_labels),
      max.overlaps = 15,
      size = 4,
      box.padding = 0.3
    ) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
    scale_colour_manual(
      values = c(
        "Up"     = "#D880D0",
        "Down"   = "#70C870",
        "NotSig" = "grey70"
      )
    ) +
    # Y-LIMIT for Male/Female ONLY
    { if (grepl("Male|Female", coef_name)) coord_cartesian(ylim = c(0, 4)) } +
    labs(
      title = coef_name,
      x = expression(log[2]("fold change")),
      y = expression(-log[10]("adjusted p-value"))
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.text  = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.position = "none"
    )
  
  list(plot = p, table = tt)
}

message("Results written to ", output_dir)

# ---------------- RESULTS + VOLCANO: main effects ----------------
for (cn in colnames(cont.matrix.main)) {
  res_main <- make_volcano(fit_main2, cn, prot_20,
                           fdr_cutoff = fdr_thr,
                           lfc_cutoff = lfc_thr)
  write.csv(res_main$table,
            file.path(output_dir, paste0("results_main_", cn, ".csv")),
            row.names = FALSE)
  save_plot(file.path(output_dir, paste0("volcano_main_", cn, ".png")),
            res_main$plot,
            width = 6.5, height = 6, dpi = 300, bg = "transparent")
}

message("Main‑effects results written to ", output_dir)

# ---------------- Design 2: Sex × Treatment interaction ----------------
design_int <- model.matrix(~ 0 + Batch + Sex * Treatment, data = meta_n20)
colnames(design_int) <- make.names(colnames(design_int))

corfit_int <- duplicateCorrelation(exprs_n20_ltr_norm,
                                   design = design_int,
                                   block  = meta_n20$Bioreplicate)

fit_int <- lmFit(exprs_n20_ltr_norm,
                 design = design_int,
                 block  = meta_n20$Bioreplicate,
                 correlation = corfit_int$consensus)

cont.matrix.int <- makeContrasts(
  DoxVsVehicle_Female = TreatmentDox,
  DoxVsVehicle_Male   = TreatmentDox + SexMale.TreatmentDox,
  SexDifference_Dox   = SexMale.TreatmentDox,
  levels = design_int
)

fit_int2 <- contrasts.fit(fit_int, cont.matrix.int)
fit_int2 <- eBayes(fit_int2, trend = TRUE)

# ---------------- RESULTS + VOLCANO: interaction model ----------------
for (cn in colnames(cont.matrix.int)) {
  res_int <- make_volcano(fit_int2, cn, prot_20,
                          fdr_cutoff = fdr_thr,
                          lfc_cutoff = lfc_thr)
  write.csv(res_int$table,
            file.path(output_dir, paste0("results_interaction_", cn, ".csv")),
            row.names = FALSE)
  save_plot(file.path(output_dir, paste0("volcano_interaction_", cn, ".png")),
            res_int$plot,
            width = 6.5, height = 6, dpi = 300, bg = "transparent")
}

message("Interaction‑model results written to ", output_dir)






