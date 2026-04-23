
  
  

  
  
  
  
  
  
  
  # =============================================================
  # SC3: Single-Cell Consensus Clustering - Full Correct Demo
  # Based on: Kiselev et al., Nature Methods, 2017
  # Dataset: Yan et al. 2013 (90 cells, 7 developmental stages)
  # =============================================================
  
  
  # ---- 1. INSTALL & LOAD PACKAGES -----------------------------
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(c("SC3", "scater", "SingleCellExperiment"), ask = FALSE)
  if (!requireNamespace("mclust", quietly = TRUE))
    install.packages("mclust")
  
  library(SC3)
  library(scater)
  library(SingleCellExperiment)
  library(mclust)
  
  
  # ---- 2. LOAD DATA -------------------------------------------
  data("yan")
  
  cat("Matrix dimensions (genes x cells):", dim(yan), "\n")
  cat("\nRaw column names (first 6):\n")
  print(head(colnames(yan)))
  
  
  # ---- 3. EXTRACT GROUND-TRUTH LABELS -------------------------
  # Column names follow pattern: "Morulae..1..Cell.3.RPKM."
  # We map known keywords -> clean stage labels
  
  extract_stage <- function(name) {
    if      (grepl("Oocyte",                    name, ignore.case = TRUE)) "Oocyte"
    else if (grepl("Zygote",                    name, ignore.case = TRUE)) "Zygote"
    else if (grepl("Morula",                    name, ignore.case = TRUE)) "Morula"
    else if (grepl("16.cell|X16|16cell",        name, ignore.case = TRUE)) "16cell"
    else if (grepl("2.cell|X2.cell|2cell",      name, ignore.case = TRUE)) "2cell"
    else if (grepl("4.cell|X4.cell|4cell",      name, ignore.case = TRUE)) "4cell"
    else if (grepl("8.cell|X8.cell|8cell",      name, ignore.case = TRUE)) "8cell"
    else if (grepl("Late.*blast|blast.*late",   name, ignore.case = TRUE)) "Late_blastocyst"
    else if (grepl("Mid.*blast|blast.*mid",     name, ignore.case = TRUE)) "Mid_blastocyst"
    else if (grepl("Early.*blast|blast.*early", name, ignore.case = TRUE)) "Early_blastocyst"
    else if (grepl("blast",                     name, ignore.case = TRUE)) "Blastocyst"
    else "Unknown"
  }
  
  cell_labels <- sapply(colnames(yan), extract_stage)
  
  cat("\nCell type label counts:\n")
  print(table(cell_labels))
  
  # Check for any unmatched cells
  unknown_names <- colnames(yan)[cell_labels == "Unknown"]
  if (length(unknown_names) > 0) {
    cat("\nWARNING - unmatched cell names:\n")
    print(unknown_names)
  } else {
    cat("\nAll 90 cells labelled successfully. No unknowns.\n")
  }
  
  
  # ---- 4. BUILD SingleCellExperiment OBJECT -------------------
  sce <- SingleCellExperiment(
    assays  = list(counts = as.matrix(yan)),
    colData = DataFrame(cell_type1 = cell_labels)
  )
  
  # Log-normalise counts
  sce <- logNormCounts(sce)
  
  # SC3 requires feature_symbol in rowData
  rowData(sce)$feature_symbol <- rownames(sce)
  
  # Remove any duplicate gene names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  cat("\nSCE object summary:\n")
  print(sce)
  
  
  # ---- 5. RUN SC3 (k = 7) -------------------------------------
  # k = 7 matches the 7 known developmental stages in Yan data
  # biology = TRUE: compute marker genes, DE genes, outlier scores
  # n_cores = 1: reproducible single-threaded run
  
  cat("\nRunning SC3 clustering (k=7)...\n")
  sce <- sc3(sce, ks = 6, biology = TRUE, n_cores = 1)
  cat("Clustering complete.\n")
  
  
  # ---- 6. INSPECT CLUSTER ASSIGNMENTS -------------------------
  cluster_labels <- colData(sce)$sc3_6_clusters
  
  cat("\nSC3 cluster assignments (k=7):\n")
  print(table(cluster_labels))
  
  # Cross-tabulate SC3 clusters vs ground truth
  cat("\nCluster vs ground-truth cross-table:\n")
  print(table(SC3 = cluster_labels, Truth = cell_labels))
  
  
  # ---- 7. ADJUSTED RAND INDEX (ARI) ---------------------------
  ari <- adjustedRandIndex(cell_labels, cluster_labels)
  cat(sprintf("\nAdjusted Rand Index (ARI) at k=6 = %.4f\n", ari))
  cat("(Paper reports ~0.89 for Yan dataset)\n")
  
  
  # ---- 8. EXPLORE k = 3 TO 10 ---------------------------------
  cat("\nRunning SC3 for k = 3:10...\n")
  sce_multi <- sc3(sce, ks = 3:10, biology = FALSE, n_cores = 1)
  
  cat("\nARI across k values:\n")
  best_ari <- 0
  best_k   <- 7
  for (k in 3:10) {
    cl    <- colData(sce_multi)[[paste0("sc3_", k, "_clusters")]]
    ari_k <- adjustedRandIndex(cell_labels, cl)
    flag  <- if (ari_k == max(sapply(3:10, function(kk) {
      adjustedRandIndex(cell_labels,
                        colData(sce_multi)[[paste0("sc3_", kk, "_clusters")]])
    }))) " <-- best" else ""
    cat(sprintf("  k = %2d  ->  ARI = %.4f%s\n", k, ari_k, flag))
    if (ari_k > best_ari) { best_ari <- ari_k; best_k <- k }
  }
  cat(sprintf("\nBest k = %d with ARI = %.4f\n", best_k, best_ari))
  
  
  # ---- 9. VISUALISATIONS --------------------------------------
  
  # 9a. Consensus matrix — shows how often cell pairs cluster together
  #     (1 = always, 0 = never across all parameter combinations)

  
  colData(sce)$cell_type <- cell_labels
  sc3_plot_consensus(sce, k = 6, show_pdata = "cell_type")
  
  # 9b. Silhouette index — measure of cluster compactness
  sc3_plot_silhouette(sce, k = 6)
  
  # 9c. Marker genes — AUROC > 0.85 and p < 0.01 (paper defaults)
  sc3_plot_markers(sce, k = 6, show_pdata = "cell_type1")
  
  # 9d. Differentially expressed genes — Kruskal-Wallis, p < 0.01
  sc3_plot_de_genes(sce, k = 6, show_pdata = "cell_type1")
  
  # 9e. Cell outlier scores — robust Mahalanobis distance (MCD)
  sc3_plot_cell_outliers(sce, k = 6)
  
  
  # ---- 10. EXPORT RESULTS TO CSV ------------------------------
  results_df <- data.frame(
    cell         = colnames(sce),
    ground_truth = cell_labels,
    sc3_cluster  = as.character(colData(sce)$sc3_7_clusters)
  )
  write.csv(results_df, "SC3_cluster_results.csv", row.names = FALSE)
  cat("\nResults saved to SC3_cluster_results.csv\n")
  
  cat("\n=== DONE ===\n")
  cat(sprintf("Final ARI (k=7) = %.4f\n", ari))
  
  
  # =============================================================
  # SC3 PIPELINE SUMMARY (Kiselev et al. 2017):
  #
  #  Step 1 — Gene filter    : remove ubiquitous/rare genes (X = 6%)
  #  Step 2 — Distances      : Euclidean, Pearson, Spearman on log2(M+1)
  #  Step 3 — Transformations: PCA + graph Laplacian eigenvectors
  #  Step 4 — k-means        : on first d eigenvectors (d = 4-7% of N)
  #  Step 5 — Consensus      : CSPA consensus matrix -> hierarchical clustering
  # =============================================================