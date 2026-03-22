##############################################################################
#  HBV-IT miRNA Analysis Pipeline — STEP 1
#  GEO Dataset Download + IT-specific Differentially Expressed Gene (DEG)
#  Datasets: GSE65359 (IT/IA/IC annotated), GSE83148, GSE84044
#  Author: YoungMin (with Claude)
#  Purpose: Extract IT-specific transcriptomic signature as seed for miRNA mapping
##############################################################################

# ─── 0. Package Setup ───────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("GEOquery", "limma", "Biobase", "WGCNA", "sva")
pkgs_cran <- c("tidyverse", "data.table", "ggplot2", "pheatmap",
               "RColorBrewer", "openxlsx", "ggrepel")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(GEOquery); library(limma); library(Biobase); library(sva)
library(tidyverse); library(data.table); library(pheatmap)
library(ggplot2); library(ggrepel); library(openxlsx)

# ─── 1. Download GEO Datasets ───────────────────────────────────────────────
# NOTE: GSE65359 = KEY dataset with IT/IA/IC phase annotation
#       GSE83148  = 122 HBV infected + 6 normal (Affymetrix U133 Plus 2.0)
#       GSE84044  = 124 CHB fibrosis patients (Affymetrix U133 Plus 2.0)

download_geo <- function(gse_id, destdir = "./GEO_data") {
  message("Downloading: ", gse_id)
  gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = destdir, getGPL = TRUE)
  if (length(gse) > 1) gse <- gse[[1]]
  return(gse)
}

dir.create("./GEO_data", showWarnings = FALSE)
dir.create("./results",  showWarnings = FALSE)

gse65359 <- download_geo("GSE65359")  # IT/IA/IC annotated — PRIMARY
gse83148  <- download_geo("GSE83148") # HBV vs Normal
gse84044  <- download_geo("GSE84044") # CHB fibrosis cohort

# ─── 2. Phenotype Extraction — GSE65359 ────────────────────────────────────
extract_pheno <- function(gse, gse_id) {
  pheno <- pData(gse)
  message("\n=== ", gse_id, " Phenotype columns ===")
  print(colnames(pheno))
  message("Sample count: ", nrow(pheno))
  return(pheno)
}

pheno65359 <- extract_pheno(gse65359, "GSE65359")
pheno83148 <- extract_pheno(gse83148, "GSE83148")
pheno84044 <- extract_pheno(gse84044, "GSE84044")

# ─── 3. Disease Phase Annotation — GSE65359 ────────────────────────────────
# GSE65359: Inspect characteristics column for IT/IA/IC labels
# Typical column name: characteristics_ch1, source_name_ch1, etc.

# Check all characteristics columns
char_cols <- grep("characteristics", colnames(pheno65359), value = TRUE)
for (col in char_cols) {
  message("\n--- Column: ", col, " ---")
  print(table(pheno65359[[col]]))
}

# Manual annotation (adjust after inspecting the output above)
# EXPECTED: IT (immune tolerant), IA (immune active), IC (immune control), HC (healthy)
annotate_phase <- function(pheno) {
  pheno_ann <- pheno %>%
    mutate(disease_phase = case_when(
      grepl("immune.tolerant|IT|tolerant", characteristics_ch1, ignore.case = TRUE) ~ "IT",
      grepl("immune.active|IA|active",     characteristics_ch1, ignore.case = TRUE) ~ "IA",
      grepl("immune.control|IC|control",   characteristics_ch1, ignore.case = TRUE) ~ "IC",
      grepl("healthy|normal|HC",           characteristics_ch1, ignore.case = TRUE) ~ "HC",
      TRUE ~ "unknown"
    ))
  message("\nPhase distribution:")
  print(table(pheno_ann$disease_phase))
  return(pheno_ann)
}

# NOTE: Run annotate_phase AFTER checking exact column names above
# pheno65359_ann <- annotate_phase(pheno65359)

# ─── 4. Expression Matrix Extraction ────────────────────────────────────────
get_exprs <- function(gse, log_transform = TRUE) {
  mat <- exprs(gse)
  if (log_transform && max(mat, na.rm = TRUE) > 100) {
    message("  Log2-transforming expression matrix...")
    mat <- log2(mat + 1)
  }
  # Remove probes with >20% missing values
  na_frac <- rowMeans(is.na(mat))
  mat <- mat[na_frac < 0.2, ]
  message("  Expression matrix: ", nrow(mat), " probes × ", ncol(mat), " samples")
  return(mat)
}

exprs65359 <- get_exprs(gse65359)
exprs83148 <- get_exprs(gse83148)
exprs84044 <- get_exprs(gse84044)

# ─── 5. DEG Analysis — GSE65359: IT vs IA, IT vs IC, IT vs HC ───────────────
run_limma_deg <- function(exprs_mat, group_vec, contrast_list,
                          fdr_cut = 0.05, logfc_cut = 0.5) {
  # group_vec: factor vector of group labels
  design <- model.matrix(~ 0 + group_vec)
  colnames(design) <- levels(group_vec)
  
  fit <- lmFit(exprs_mat, design)
  
  # Build contrast matrix
  cont_str  <- paste(names(contrast_list),
                     sapply(contrast_list, function(x) paste(x[1], x[2], sep = "-")),
                     sep = "=")
  cont_mat  <- makeContrasts(contrasts = cont_str, levels = design)
  fit2      <- contrasts.fit(fit, cont_mat)
  fit2      <- eBayes(fit2)
  
  # Extract results for each contrast
  results_list <- lapply(names(contrast_list), function(cname) {
    tt <- topTable(fit2, coef = paste0(cname, "=", contrast_list[[cname]][1],
                                       "-", contrast_list[[cname]][2]),
                   n = Inf, sort.by = "P")
    tt$contrast <- cname
    tt$sig      <- tt$adj.P.Val < fdr_cut & abs(tt$logFC) > logfc_cut
    return(tt)
  })
  names(results_list) <- names(contrast_list)
  return(results_list)
}

# ─── 6. IT-Specific DEG Filtering ───────────────────────────────────────────
# Strategy: 
#   IT_core = significant in IT vs IA AND IT vs HC
#   NOT significant in IA vs IC (i.e., truly IT-specific, not just HBV-general)

filter_IT_specific <- function(deg_IT_vs_IA, deg_IT_vs_HC, deg_IA_vs_IC,
                                fdr_cut = 0.05, logfc_cut = 0.5) {
  
  # Genes significant in IT vs IA
  it_vs_ia_sig <- deg_IT_vs_IA %>%
    filter(adj.P.Val < fdr_cut, abs(logFC) > logfc_cut) %>%
    pull(ID)
  
  # Genes significant in IT vs HC  
  it_vs_hc_sig <- deg_IT_vs_HC %>%
    filter(adj.P.Val < fdr_cut, abs(logFC) > logfc_cut) %>%
    pull(ID)
  
  # Genes significant in IA vs IC (to EXCLUDE — not IT-specific)
  ia_vs_ic_sig <- deg_IA_vs_IC %>%
    filter(adj.P.Val < fdr_cut, abs(logFC) > logfc_cut) %>%
    pull(ID)
  
  # IT-core: present in both IT comparisons, absent from IA vs IC
  it_core <- intersect(it_vs_ia_sig, it_vs_hc_sig)
  it_core_specific <- setdiff(it_core, ia_vs_ic_sig)
  
  message("IT vs IA significant genes:      ", length(it_vs_ia_sig))
  message("IT vs HC significant genes:       ", length(it_vs_hc_sig))
  message("IT-core (both comparisons):       ", length(it_core))
  message("IT-specific (not in IA vs IC):    ", length(it_core_specific))
  
  # Return annotated table
  it_specific_df <- deg_IT_vs_IA %>%
    filter(ID %in% it_core_specific) %>%
    mutate(it_specific = TRUE)
  
  return(it_specific_df)
}

# ─── 7. Cross-Dataset Validation ────────────────────────────────────────────
# Use GSE83148 (HBV vs Normal) to validate IT-core genes
# Strategy: IT-core genes should show consistent directional DE in HBV vs Normal

cross_validate <- function(it_specific_genes, exprs_validation, pheno_validation,
                           normal_label = "normal", hbv_label = "HBV") {
  
  # Subset to shared genes
  shared_genes <- intersect(rownames(exprs_validation), it_specific_genes$ID)
  message("Genes shared with validation dataset: ", length(shared_genes))
  
  # Simple DE in validation
  group_val <- factor(ifelse(grepl(normal_label, pheno_validation$source_name_ch1,
                                   ignore.case = TRUE), "Normal", "HBV"))
  design_val <- model.matrix(~ group_val)
  fit_val    <- lmFit(exprs_validation[shared_genes, ], design_val)
  fit_val    <- eBayes(fit_val)
  tt_val     <- topTable(fit_val, coef = 2, n = Inf)
  
  # Directional concordance
  merged <- it_specific_genes %>%
    filter(ID %in% shared_genes) %>%
    left_join(tt_val %>% rownames_to_column("ID") %>% select(ID, logFC, adj.P.Val),
              by = "ID", suffix = c(".discovery", ".validation"))
  
  concordant <- merged %>%
    filter(sign(logFC.discovery) == sign(logFC.validation),
           adj.P.Val.validation < 0.1)
  
  message("Directionally concordant genes: ", nrow(concordant),
          " (", round(100 * nrow(concordant)/nrow(merged), 1), "%)")
  
  return(list(merged = merged, concordant = concordant))
}

# ─── 8. Visualization: Volcano Plots ────────────────────────────────────────
plot_volcano <- function(deg_df, title = "Volcano Plot",
                          fdr_cut = 0.05, logfc_cut = 0.5,
                          top_n_labels = 20) {
  deg_df <- deg_df %>%
    mutate(
      color = case_when(
        adj.P.Val < fdr_cut & logFC >  logfc_cut ~ "Up in IT",
        adj.P.Val < fdr_cut & logFC < -logfc_cut ~ "Down in IT",
        TRUE ~ "NS"
      ),
      neg_log10_p = -log10(adj.P.Val + 1e-300)
    )
  
  # Top genes to label
  top_genes <- deg_df %>%
    filter(color != "NS") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n_labels)
  
  ggplot(deg_df, aes(x = logFC, y = neg_log10_p, color = color)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_text_repel(data = top_genes, aes(label = ID),
                    size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c("Up in IT"    = "#D85A30",
                                  "Down in IT"  = "#185FA5",
                                  "NS"          = "#888780")) +
    geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    labs(title = title,
         x = "log2 Fold Change (IT vs comparison)",
         y = "-log10(FDR)",
         color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(size = 12, face = "bold"))
}

# ─── 9. Save IT-Specific Gene List ──────────────────────────────────────────
save_IT_genes <- function(it_specific_df, filename = "./results/IT_specific_genes.xlsx") {
  # This list becomes the SEED for Step 2 (miRNA target reverse mapping)
  wb <- createWorkbook()
  addWorksheet(wb, "IT_specific_DEG")
  writeData(wb, "IT_specific_DEG", it_specific_df)
  
  # Summary sheet
  addWorksheet(wb, "Summary")
  summary_df <- data.frame(
    Category = c("Total IT-specific genes", "Upregulated in IT", "Downregulated in IT",
                 "FDR threshold", "LogFC threshold"),
    Value    = c(nrow(it_specific_df),
                 sum(it_specific_df$logFC > 0),
                 sum(it_specific_df$logFC < 0),
                 "0.05", "0.5")
  )
  writeData(wb, "Summary", summary_df)
  saveWorkbook(wb, filename, overwrite = TRUE)
  message("Saved: ", filename)
  message("→ This gene list feeds directly into Step 2 (miRNA reverse mapping)")
}

# ─── 10. Master Run Function ─────────────────────────────────────────────────
# IMPORTANT: Run step-by-step; inspect output at each stage
# before proceeding, especially the pheno column inspection

run_step1 <- function() {
  message("\n========================================")
  message("  STEP 1: GEO + IT-DEG ANALYSIS")
  message("  Primary: GSE65359 (IT/IA/IC annotated)")
  message("  Validation: GSE83148, GSE84044")
  message("========================================\n")
  
  message("[1/5] Download GEO datasets...")
  # (already done above)
  
  message("\n[2/5] Inspect GSE65359 phenotype columns...")
  message(">>> ACTION REQUIRED: Review the column output and adjust annotate_phase() <<<")
  char_cols <- grep("characteristics|source|title", colnames(pheno65359), value = TRUE)
  for (col in char_cols) {
    cat("\n--- Column:", col, "---\n")
    print(table(pheno65359[[col]]))
  }
  
  message("\n[3/5] After inspecting columns, update annotate_phase() and run DEG...")
  message("  Contrasts to run: IT vs IA, IT vs HC, IA vs IC")
  
  message("\n[4/5] Filter IT-specific genes...")
  message("  IT-specific = (IT vs IA sig) AND (IT vs HC sig) NOT (IA vs IC sig)")
  
  message("\n[5/5] Cross-validate in GSE83148 and GSE84044...")
  message("  Concordance threshold: same direction + FDR < 0.1")
  
  message("\n✓ Step 1 complete. Output: ./results/IT_specific_genes.xlsx")
  message("→ Proceed to Step 2: miRNA reverse mapping using IT gene seed list")
}

run_step1()
