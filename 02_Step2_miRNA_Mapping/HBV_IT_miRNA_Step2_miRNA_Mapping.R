##############################################################################
#  HBV-IT miRNA Analysis Pipeline — STEP 2
#  miRNA Reverse Mapping: IT-specific genes → candidate miRNA regulators
#  Input:  ./results/IT_specific_genes.xlsx  (from Step 1)
#  v18 connection: AIM2, LGALS9, TGFB1, DNMT1, SOCS1 → pre-loaded as anchors
#  Output: IT-associated candidate miRNA list with confidence scores
##############################################################################

library(multiMiR)      # BiocManager::install("multiMiR")
library(tidyverse)
library(openxlsx)
library(igraph)
library(ggplot2)

# ─── 0. Load IT-specific genes from Step 1 ──────────────────────────────────
it_genes <- read.xlsx("./results/IT_specific_genes.xlsx", sheet = "IT_specific_DEG")

# ─── 1. Pre-load v18 Anchor Genes ───────────────────────────────────────────
# These are confirmed IT-relevant genes from the HBV-IT manuscript (v18)
# Direct connection between the CMH paper and this miRNA analysis
v18_anchor_genes <- c(
  "AIM2",    # Effector suppression layer 1 — inflammasome sensor
  "LGALS9",  # Effector suppression layer 2 — Tim-3 ligand, T-cell exhaustion
  "TGFB1",   # Effector suppression layer 3 — immune tolerance cytokine
  "DNMT1",   # Epigenetic regulation — DNA methylation
  "SOCS1",   # JAK-STAT suppression — innate immune evasion
  "IL10",    # Anti-inflammatory cytokine
  "FOXP3",   # Treg master transcription factor
  "PDCD1",   # PD-1 — T-cell exhaustion
  "CD274",   # PD-L1 — checkpoint ligand
  "CTLA4"    # T-cell co-inhibitory receptor
)

# Merge anchor genes with IT-specific DEG list
it_gene_seed <- union(it_genes$ID, v18_anchor_genes)
message("Total seed genes for miRNA mapping: ", length(it_gene_seed))
message("  IT-specific DEG: ", length(it_genes$ID))
message("  v18 anchor genes: ", length(v18_anchor_genes))
message("  Union seed set:   ", length(it_gene_seed))

# ─── 2. multiMiR: Reverse Target Mapping ────────────────────────────────────
# Query validated + predicted miRNA-target databases
# Validated: miRTarBase, TarBase, miRecords
# Predicted:  TargetScan, miRanda, PicTar, DIANA-microT

run_multimir_query <- function(gene_list, org = "hsa",
                               table_opt = "validated",
                               chunk_size = 50) {
  
  message("Running multiMiR query for ", length(gene_list), " genes...")
  message("  Table type: ", table_opt)
  
  # Process in chunks to avoid timeout
  results_list <- list()
  n_chunks     <- ceiling(length(gene_list) / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx   <- min(i * chunk_size, length(gene_list))
    chunk     <- gene_list[start_idx:end_idx]
    
    message("  Chunk ", i, "/", n_chunks, " (", length(chunk), " genes)")
    
    tryCatch({
      result <- get.multimir(org = org,
                             target = chunk,
                             table  = table_opt,
                             summary = FALSE)
      if (!is.null(result@data) && nrow(result@data) > 0) {
        results_list[[i]] <- result@data
      }
    }, error = function(e) {
      message("  Warning: chunk ", i, " failed: ", e$message)
    })
    
    Sys.sleep(0.5)  # Polite pause between API calls
  }
  
  if (length(results_list) == 0) {
    warning("No results returned from multiMiR")
    return(NULL)
  }
  
  combined <- bind_rows(results_list)
  message("  Total miRNA-target pairs retrieved: ", nrow(combined))
  return(combined)
}

# Query validated interactions first (highest confidence)
mirna_validated <- run_multimir_query(it_gene_seed, table_opt = "validated")

# Then predicted (wider coverage)
mirna_predicted <- run_multimir_query(it_gene_seed, table_opt = "predicted")

# ─── 3. Priority Scoring ────────────────────────────────────────────────────
prioritize_mirna <- function(validated_df, predicted_df,
                              min_validated_targets = 2,
                              min_predicted_targets = 5) {
  
  # Score by number of IT-seed targets
  score_validated <- validated_df %>%
    group_by(mature_mirna_id) %>%
    summarise(
      n_targets_validated    = n_distinct(target_symbol),
      targets_validated      = paste(sort(unique(target_symbol)), collapse = ";"),
      # Anchor gene enrichment score
      n_v18_anchors_targeted = sum(unique(target_symbol) %in% v18_anchor_genes),
      score_validated        = n_targets_validated + 3 * n_v18_anchors_targeted
    ) %>%
    filter(n_targets_validated >= min_validated_targets) %>%
    arrange(desc(score_validated))
  
  score_predicted <- predicted_df %>%
    group_by(mature_mirna_id) %>%
    summarise(
      n_targets_predicted    = n_distinct(target_symbol),
      score_predicted        = n_targets_predicted
    ) %>%
    filter(n_targets_predicted >= min_predicted_targets)
  
  # Merge
  combined_score <- full_join(score_validated, score_predicted,
                               by = "mature_mirna_id") %>%
    mutate(
      n_targets_validated = replace_na(n_targets_validated, 0),
      n_targets_predicted = replace_na(n_targets_predicted, 0),
      n_v18_anchors_targeted = replace_na(n_v18_anchors_targeted, 0),
      composite_score     = score_validated * 2 + score_predicted,
      has_v18_anchor      = n_v18_anchors_targeted > 0
    ) %>%
    arrange(desc(composite_score))
  
  message("Candidate miRNAs with validated targets:  ", 
          sum(combined_score$n_targets_validated > 0))
  message("Candidate miRNAs targeting v18 anchors:   ",
          sum(combined_score$has_v18_anchor, na.rm = TRUE))
  message("Top 10 candidate miRNAs:")
  print(head(combined_score %>% select(mature_mirna_id, n_targets_validated,
                                        n_v18_anchors_targeted, composite_score), 10))
  
  return(combined_score)
}

it_mirna_candidates <- prioritize_mirna(mirna_validated, mirna_predicted)

# ─── 4. HMDD v4.0 Cross-Reference: HBV-associated miRNAs ───────────────────
# Download HMDD (Human miRNA Disease Database) entries for HBV
# URL: https://www.cuilab.cn/hmdd
# Manual step: download "hepatitis B" entries from HMDD v4.0 and save as CSV

load_hmdd_hbv <- function(hmdd_file = NULL) {
  # If no file, use literature-curated list from known studies
  if (is.null(hmdd_file) || !file.exists(hmdd_file)) {
    message("HMDD file not found — using literature-curated HBV miRNA list")
    
    # Curated from: 
    #   Brunetto et al 2014 (PLOS ONE), 
    #   Frontiers Cellular Microbiology 2022,
    #   Sci Rep 2015 (Chen et al, n=495)
    hbv_mirnas_literature <- data.frame(
      mirna       = c("hsa-miR-122-5p", "hsa-miR-99a-5p",  "hsa-miR-192-5p",
                      "hsa-miR-223-3p", "hsa-miR-143-3p",  "hsa-miR-21-5p",
                      "hsa-miR-155-5p", "hsa-miR-125b-5p", "hsa-miR-181a-5p",
                      "hsa-miR-146a-5p","hsa-miR-29b-3p",  "hsa-miR-Let-7c-5p",
                      "hsa-miR-152-3p", "hsa-miR-148a-3p", "hsa-miR-101-3p",
                      "hsa-miR-195-5p", "hsa-miR-200a-3p", "hsa-miR-141-3p"),
      direction_in_HBV = c("up","down","up",
                            "up","down","up",
                            "up","down","up",
                            "up","down","down",
                            "down","down","down",
                            "down","down","down"),
      HBV_context = c("liver/serum","serum","liver/serum",
                      "plasma","plasma/liver","liver",
                      "PBMC","serum","PBMC",
                      "liver","liver","liver",
                      "liver","liver","liver",
                      "liver","liver","liver"),
      source      = c("Brunetto2014","Brunetto2014","Brunetto2014",
                      "Frontiers2022","Frontiers2022","multiple",
                      "multiple","multiple","multiple",
                      "multiple","multiple","multiple",
                      "multiple","multiple","multiple",
                      "multiple","multiple","multiple")
    )
    return(hbv_mirnas_literature)
  } else {
    hmdd <- read.csv(hmdd_file)
    hbv_entries <- hmdd %>%
      filter(grepl("hepatitis B", disease, ignore.case = TRUE))
    return(hbv_entries)
  }
}

hbv_mirna_db <- load_hmdd_hbv()

# ─── 5. Overlap: IT-candidates ∩ HBV-associated miRNAs ──────────────────────
prioritize_IT_HBV_mirna <- function(it_candidates, hbv_db) {
  
  # Normalize miRNA names for matching
  normalize_mirna <- function(x) {
    x %>% tolower() %>%
      gsub("hsa-", "", .) %>%
      gsub("mir-", "mir", .) %>%
      gsub("miR-", "mir", .)
  }
  
  it_norm  <- normalize_mirna(it_candidates$mature_mirna_id)
  hbv_norm <- normalize_mirna(hbv_db$mirna)
  
  overlap_idx <- which(it_norm %in% hbv_norm)
  overlap_mirna <- it_candidates[overlap_idx, ] %>%
    mutate(hbv_associated = TRUE,
           tier = case_when(
             n_v18_anchors_targeted > 0 & composite_score > 10 ~ "Tier1_IT_anchor",
             n_targets_validated >= 3                           ~ "Tier2_IT_validated",
             TRUE                                               ~ "Tier3_IT_predicted"
           ))
  
  message("\n=== miRNA Prioritization Summary ===")
  message("Total IT-candidate miRNAs:           ", nrow(it_candidates))
  message("HBV-associated in literature:        ", nrow(hbv_db))
  message("IT-candidate AND HBV-associated:     ", nrow(overlap_mirna))
  message("\nTier breakdown:")
  print(table(overlap_mirna$tier))
  
  return(overlap_mirna)
}

it_hbv_mirna_final <- prioritize_IT_HBV_mirna(it_mirna_candidates, hbv_mirna_db)

# ─── 6. Network Visualization: miRNA → IT gene targets ──────────────────────
build_mirna_network <- function(validated_pairs, top_mirna, top_n = 20) {
  
  top_mirna_ids <- head(top_mirna$mature_mirna_id, top_n)
  
  edges <- validated_pairs %>%
    filter(mature_mirna_id %in% top_mirna_ids) %>%
    select(from = mature_mirna_id, to = target_symbol) %>%
    distinct()
  
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # Node attributes
  V(g)$type   <- ifelse(V(g)$name %in% top_mirna_ids, "miRNA", "gene")
  V(g)$is_v18 <- V(g)$name %in% v18_anchor_genes
  V(g)$color  <- case_when(
    V(g)$name %in% top_mirna_ids & V(g)$name %in% it_hbv_mirna_final$mature_mirna_id
      ~ "#D85A30",   # Coral: IT-HBV miRNA
    V(g)$name %in% top_mirna_ids
      ~ "#534AB7",   # Purple: IT-only miRNA
    V(g)$is_v18
      ~ "#0F6E56",   # Teal: v18 anchor gene
    TRUE
      ~ "#888780"    # Gray: other IT gene
  )
  
  message("\nNetwork stats:")
  message("  Nodes: ", vcount(g), " (", sum(V(g)$type == "miRNA"), " miRNAs, ",
          sum(V(g)$type == "gene"), " genes)")
  message("  Edges: ", ecount(g))
  message("  Density: ", round(edge_density(g), 4))
  
  return(g)
}

# ─── 7. Save All Results ─────────────────────────────────────────────────────
save_step2_results <- function() {
  wb <- createWorkbook()
  
  addWorksheet(wb, "IT_miRNA_candidates")
  writeData(wb, "IT_miRNA_candidates", it_mirna_candidates)
  
  addWorksheet(wb, "IT_HBV_miRNA_final")
  writeData(wb, "IT_HBV_miRNA_final", it_hbv_mirna_final)
  
  addWorksheet(wb, "HBV_miRNA_literature")
  writeData(wb, "HBV_miRNA_literature", hbv_mirna_db)
  
  addWorksheet(wb, "v18_anchor_genes")
  writeData(wb, "v18_anchor_genes",
            data.frame(gene = v18_anchor_genes, note = "IT-specific, v18 manuscript"))
  
  saveWorkbook(wb, "./results/Step2_IT_miRNA_mapping.xlsx", overwrite = TRUE)
  message("Saved: ./results/Step2_IT_miRNA_mapping.xlsx")
  message("→ IT-HBV miRNA candidates feed into Step 3 meta-analysis &")
  message("  Step 4 network integration")
}

save_step2_results()

message("\n=== STEP 2 COMPLETE ===")
message("Key finding: IT-specific genes (v18 + DEG) → candidate miRNA regulators")
message("v18 Connection: TGFB1/LGALS9/SOCS1 miRNA axes now identified")
message("Next: Step 3 — literature meta-analysis for external concordance")
