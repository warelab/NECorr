#' NECorr
#' @author Christophe Liseron-Monfils, Andrew Olson
#' @param expression file or data.table of Expression file in log2 (ratio expression) with row: gene,
#' first column: type of sample,second column: sample names
#' @param networkFile file or data.table of Molecular network file with source in the first column,
#' targets in the second column
#' with the same name than the network file and the same localization
#' @param description.file file or data.table of the genome description
#' @param condition Condition from expression to study the network co-expression correlation
#' @param meta_data file or data.table with the meta_data
#' @param permutation permutation number used for all significance calculation
#' @param method correlation method: GCC, PCC, KCC, SCC, SPP; c("GCC", "PCC", "SCC", "KCC")
#' @param visualize TRUE/FALSE whether to generate visualizations of the results
#' @param output_dir Optional directory to save results and plots if
#' @param interactive_net TRUE/FALSE whether to generate an interactive network
#' @param top_n Number of top hub genes to highlight in plots
#' @param useBestGCC TRUE/FALSE whether to use the best GCC method
#' @param asymmetricGCC TRUE/FALSE whether to use asymmetric GCC, only if useBestGCC is TRUE, see multiCorr documentation
#' @param save_results TRUE/FALSE whether to save results to CSV files
#' @param ... Additional arguments passed to correlation functions.
#' @description NECorr helps discover candidate genes that could be
#' important for specific conditions.
#' The principal inputs are the expression data and the network file.
#' The expression data should start with 3 header columns.
#' The first column describes the conditions. Each condition will be
#' treated separately for the co-expression analysis
#' The output of the program will be generated in a result folder generated
#' in the working path
#' Create the output directory if not existing; generate "./results" dir and
#' "./results/tmp"
#' C.Liseron-Monfils - Ware lab Sept2013 - CSHL
#' partly based on rsgcc package for the GCC, PCC,KCC and SCC correlation
#' Ma et al, 2012, plant Physiology
#' @return res
#' @export
NECorr <- function(networkFile, expression, description.file,
                   condition, meta_data, permutation = 1000,
                   method="GCC", visualize = TRUE,
                   output_dir = NULL, interactive_net = FALSE,
                   top_n = 20, save_results = TRUE,
                   useBestGCC = FALSE, asymmetricGCC = FALSE, ...) {

  requireNamespace("data.table", quietly = TRUE)
  requireNamespace("matrixStats", quietly = TRUE)
  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("stats", quietly = TRUE)
  n_threads <- detectCores(logical = TRUE) - 2
  if (n_threads < 1) n_threads <- 1
  # Define the steps of the process
  steps <- c(
    "Loading network",
    "Loading expression",
    "Loading description",
    "Loading metadata   ",
    "Filtering expression",
    "Calculating correlation",
    "Tissue specificity     ",
    "Integrating p-values.  ",
    "Aggregating gene metrics",
    "DE ranking              ",
    "Network topology.       ",
    "Merging parameters.     ",
    "Hub ranking             ",
    "Edge significance       ",
    "Visualization           ",
    "Saving results          "
  )

  # update_progress <- function(step) {
  #   cat(sprintf("\n[Step %d] %s", step, steps[step]))
  #   flush.console()
  # }
  total_steps <- length(steps)
  start_time <- Sys.time()
  # Progress bar function - Define the update_progress function
  update_progress <- function(step) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    percent <- step / total_steps
    eta <- (elapsed / step) * (total_steps - step)
    bar_width <- 40
    filled <- round(bar_width * percent)
    #bar <- paste0(strrep("█", filled),
    #              strrep("░", bar_width - filled))
    bar <- paste0(strrep("#", filled),
                  strrep("-", bar_width - filled))
    eta_str <- sprintf("%02d:%02d", floor(eta / 60), round(eta %% 60))
    cat(sprintf("\r[%s] %3d%% | ETA: %s | %s", bar, round(percent * 100),
                eta_str, steps[step]))
    flush.console()
    if (step == total_steps) {
      cat("\n")
    }
  }


  #### 1. Load network ####
  # --- Load network ---
  if (is.character(networkFile) && file.exists(networkFile)) {
    network.int <- fread(networkFile, header = FALSE,
                         encoding = "Latin-1",
                         colClasses = c("character", "character"),
                         nThread = n_threads, showProgress = FALSE)
    if (ncol(network.int) == 2) {
      setnames(network.int, c("source", "target"))
    } else if (ncol(network.int) >= 3) {
      # Take only columns 1 and 3
      network.int <- network.int[, .(source = V1, target = V3)]
    } else {
      stop("networkFile must have at least 2 columns.")
    }
  } else if (is.data.frame(networkFile) || is.data.table(networkFile)) {
    network.int <- safe_as_data_table(networkFile)
    # Ensure correct colnames
    if (!all(c("source", "target") %in% names(network.int))) {
      setnames(network.int, 1:2, c("source", "target"))
    }
  } else {
    stop("networkFile must be a file path or a data.frame/data.table")
  }
  # Remove NA edges
  network.int <- network.int[complete.cases(source, target)]
  update_progress(1) # Loading network
  #print(paste("Network loaded with", nrow(network.int), "edges."))
  #### 2. Load expression ####
  if (is.character(expression) && file.exists(expression)) {
    eset <- fread(expression, header = TRUE,
                  nThread = n_threads, showProgress = FALSE)
    gene_ids <- eset[[1]]
    eset[, 1 := NULL]
    eset[, gene_id := gene_ids]
    setkey(eset, gene_id)
  } else { eset <- safe_as_data_table(expression) }
  update_progress(2) # Loading expression
  #print(paste("Expression data loaded with", nrow(eset), "genes and", ncol(eset)-1, "samples."))
  #### 3. Load description ####
  if (is.character(description.file) && file.exists(description.file)) {
    description_df <- fread(description.file, header = TRUE,
                            nThread = n_threads, showProgress = FALSE)
    colnames(description_df)[1] <- "gene_id"
  } else {
    description_df <- safe_as_data_table(description.file)
    if (!"gene_id" %in% colnames(description_df)) {
      stop("Description data must have a 'gene_id' column")
    }
  }
  setkey(description_df, gene_id)
  update_progress(3) # Loading description
  #print(paste("Description data loaded with", nrow(description_df), "entries."))
  #### 4. Load metadata ####
  if (is.character(meta_data) && file.exists(meta_data)) {
    factortab <- fread(meta_data, header = TRUE,
                       nThread = n_threads, showProgress = FALSE)
    # First column is sample IDs
    sample_ids <- factortab[[1]]
    # Remove the first column from the data.table
    factortab[, 1 := NULL]
    # Rename the current first column to "condition"
    setnames(factortab, 1, "condition")
    # Add sample_id as a separate column
    factortab[, sample_id := sample_ids]
    # Set rownames for internal use
    rownames(factortab) <- sample_ids
  } else if (is.data.frame(meta_data) || is.data.table(meta_data)) {
    factortab <- meta_data
    # Rename second column to "condition"
    setnames(factortab, 1, "condition")
    # Set rownames from sample_id
    factortab[, sample_id := rownames(factortab)]
    setcolorder(factortab, c("sample_id", "condition"))
  } else {
      stop("Metadata must have at least two columns: sample_id and condition")
  }
    # Verify sample IDs match expression column names (excluding gene_id)
    expr_sample_cols <- setdiff(colnames(eset), "gene_id")
    if (!setequal(rownames(factortab), expr_sample_cols)) {
      stop("Mismatch between metadata sample IDs and expression data column names.")
  }
  update_progress(4) # Loading metadata
  #print(paste("Metadata loaded with", nrow(factortab), "samples."))
  #### 5. Check sample existence ####
  if (!all(factortab$sample_id %in% colnames(eset)[-which(names(eset) == "gene_id")])) {
    stop("Some samples in meta_data are not present in the expression data")
  }

  #### 6. Prepare data ####
  Genelist <- unique(c(network.int$source, network.int$target))
  methods <- toupper(method)
  valid_methods <- c("GCC", "PCC", "SCC", "KCC", "SPP")
  if (!all(methods %in% valid_methods)) stop("Invalid methods specified.")

  expr_cols <- setdiff(names(eset), "gene_id")
  eset[, (expr_cols) := lapply(.SD, as.numeric), .SDcols = expr_cols]
  eset[, (expr_cols) := lapply(.SD, function(x) fifelse(is.infinite(x), 1e300, x)), .SDcols = expr_cols]
  #print(paste("Expression data prepared with", nrow(eset), "genes."))
  #### 7. Handle replicates ####
  # --- Handle replicates ---
  sample.names <- unique(factortab$condition)
  idx_list <- unlist(lapply(sample.names, function(sn) {
    idx <- which(factortab$condition == sn)
    if (length(idx) == 1) rep(idx, 3) else idx
  }))

  # Create factor for DE analysis
  treatment.f <- factor(rep(sample.names, times = ifelse(table(factortab$condition) == 1, 3, table(factortab$condition))),
                        levels = sample.names)

  #### 8. Filter expression matrix ####
  eset <- eset[gene_id %in% Genelist]
  m.eset <- as.matrix(eset[, ..expr_cols])[ , idx_list, drop = FALSE]
  rownames(m.eset) <- eset$gene_id
  # Tissues vector matches m.eset columns
  tissues <- factortab$condition[idx_list]
  update_progress(5) # Filtering expression
  #print(paste("Expression data filtered to", nrow(m.eset), "genes and", ncol(m.eset), "samples."))
  #### 9. Correlation (multiCorr) ####
  int.sig.full <- multiCorr(expression = m.eset, edges = network.int,
                            pernum = permutation, methods = method, useBestGCC = useBestGCC, asymmetricGCC = asymmetricGCC)
  update_progress(6) #Calculating correlation
  # print(paste("Correlation calculated for", nrow(int.sig.full), "edges using methods:", paste(methods, collapse = ", ")))
  #### 10. Tissue specificity ####
  ts.res <- ts.IUT(name = paste(condition, permutation, sep = "_"),
                   eset = m.eset, tissues = tissues,
                   target = condition, threshold = 0.05, filter = 3, minGenes = 5)
  tsi.order <- if (!is.null(ts.res)) {
    #message("Head of tissue specificity index (TSI)")
    #print(head(ts.res$filtTSI))
    #message("Head of tissue specificity index (TSE)")
    #print(head(ts.res$filtTSE))
    colnames(ts.res$filtTSE) <- c("gene", "TSI")
    ts <- rbind(ts.res$filtTSI, ts.res$filtTSE)
    ts <- ts[!duplicated(ts$gene), ]
    ScalN(setNames(ts$TSI, ts$gene))
  } else {
    data.frame(Gene = Genelist, TSI = 0)
  }
  colnames(tsi.order)[2] <- "TSI"
  update_progress(7)
  # print(paste("Tissue specificity calculated for", nrow(tsi.order), "genes."))
  #### 11. Integrate p-values (no rowwise) ####
  N_edges <- nrow(int.sig.full)
  rank_matrix <- safe_as_data_table(lapply(int.sig.full[, ..methods], function(x) frankv(-x, ties.method = "average")))
  AvgRank <- rowMeans(as.matrix(rank_matrix), na.rm = TRUE)
  Correlation <- 1 - 2 * ((AvgRank - 1) / (N_edges - 1))

  pval_cols <- paste0(methods, "_pvalue")
  log_p <- log(as.matrix(int.sig.full[, ..pval_cols]))
  Fisher_chi2 <- -2 * rowSums(log_p, na.rm = TRUE)
  Fisher_df <- 2 * rowSums(!is.na(log_p))
  pvalue <- pchisq(Fisher_chi2, df = Fisher_df, lower.tail = FALSE)

  int.sig <- data.table(source = int.sig.full$source,
                        target = int.sig.full$target,
                        Correlation = Correlation,
                        pvalue = pvalue)
  update_progress(8)
  # print(paste("P-values integrated for", nrow(int.sig), "edges."))
  #### 12. Aggregate to gene level ####
  edge_long <- rbind(
    int.sig[, .(Gene = source, Correlation, pvalue)],
    int.sig[, .(Gene = target, Correlation, pvalue)]
  )
  agg <- edge_long[, .(
    mean_abs_corr = mean(abs(Correlation), na.rm = TRUE),
    min_pval = fishersMethod(pvalue)
  ), by = Gene]

  coexpr_df <- ScalN(setNames(agg$mean_abs_corr, agg$Gene))
  neglogp_df <- ScalN(setNames(-log10(agg$min_pval), agg$Gene))
  coexprs.pvals <- merge(coexpr_df, neglogp_df, by = "Gene", all = TRUE)
  setnames(coexprs.pvals, c("Gene", "mean_abs_corr", "coexprs.pvals"))
  update_progress(9)
  # print(paste("Gene-level metrics aggregated for", nrow(coexprs.pvals), "genes."))
  #### 13. DE ranking ####
  DE.ranks <- DE.ranking(m.eset, Genelist, treatment.f, condition, sample.names)
  if (is.null(DE.ranks)) DE.ranks <- data.frame(Gene = Genelist, DE.ranks = 0)
  update_progress(10)
  # print(paste("DE ranking calculated for", nrow(DE.ranks), "genes."))
  #### 14. Network topology with igraph ####
  g <- graph_from_data_frame(network.int, directed = FALSE)
  netstat <- data.table(
    Gene = V(g)$name,
    BetwC = betweenness(g),
    Conn = degree(g),
    ClusCoef = transitivity(g, type = "local", isolates = "zero"),
    EigenC = eigen_centrality(g)$vector,
    PageRank = page_rank(g)$vector
  )
  # Scale numeric columns between 0 and 1

  cols_to_scale <- c("BetwC", "Conn", "ClusCoef", "EigenC", "PageRank")
  netstat[, (cols_to_scale) := lapply(.SD, function(x) {
    if (all(is.na(x))) {
      return(x) # keep NA if column is all NA
    } else if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      return(rep(0, length(x))) # avoid division by zero if constant
    } else {
      return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    }
  }), .SDcols = cols_to_scale]
  update_progress(11)
  # print(paste("Network topology calculated for", nrow(netstat), "genes."))
  #### 15. Merge parameters ####
  m.param <- mergeMetricsByGene(coexprs.pvals, tsi.order, netstat, DE.ranks)
  m.param[is.na(m.param)] <- 0
  if (inherits(m.param, "data.table")) {
    m.param <- as.data.frame(m.param)
  }
  rownames(m.param) <- m.param$Gene
  m.param <- m.param[, -which(names(m.param) == "Gene"),
                     drop = FALSE]
  #
  #   m.param <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
  #                     list(coexprs.pvals, tsi.order, netstat, DE.ranks))
  #   m.param[is.na(m.param)] <- 0
  #   update_progress(12)
  #   if (inherits(m.param, "data.table")) {
  #     m.param <- as.data.frame(m.param)
  #   }
  #   if ("Gene" %in% colnames(m.param)) {
  #     rownames(m.param) <- m.param$Gene
  #     m.param$Gene <- NULL
  #   }
  #### 16. Hub ranking ####
  weights <- entropy_weight(m.param)
  hub.m.param <- sweep(as.matrix(m.param), 2, weights, `*`)
  gene.rank.h <- sort(rowSums(hub.m.param), decreasing = TRUE)
  gene_rank_df <- data.frame(Gene = names(gene.rank.h), Score = gene.rank.h)
  update_progress(13)
  # print(paste("Hub ranking calculated for", nrow(gene_rank_df), "genes."))
  #### 17. Edge significance ####
  hub.int.ranks <- hub_edge_significance(as.data.frame(network.int), as.list(gene.rank.h))
  hub.int.significant <- hub.int.ranks[hub.int.ranks$p2 < 0.05, , drop = FALSE]

  actres <- if (nrow(hub.int.significant) > 0)
    activator_significant(hub.int.significant, as.data.frame(network.int), as.data.frame(description_df)) else NULL
  update_progress(14) # Edge significance
  # print(paste("Edge significance calculated with", nrow(hub.int.significant), "significant edges."))
  #### 18. Compile results ####
  results <- list(
    necorrHub_nodes = gene_rank_df,
    necorrReg = actres,
    necorrEdges = hub.int.significant,
    netstat = as.data.frame(netstat),
    coexpres = as.data.frame(int.sig.full),
    hub.m.param = m.param
  )
  #### 19. Visualization ####
  if (visualize) {

    # Check for empty network results
    edges_empty <- is.null(results$necorrEdges) || nrow(results$necorrEdges) == 0
    hubs_empty  <- is.null(results$necorrHub_nodes) || nrow(results$necorrHub_nodes) == 0
    if (edges_empty || hubs_empty) {
      message("Visualization skipped: No network data available.\n",
              "- necorrEdges empty: ", edges_empty, "\n",
              "- necorrHub_nodes empty: ", hubs_empty)
      results$plots <- NULL   # or NA, depending on how you use it
    } else {
      # Safe call only when data exists
      results$plots <- visualize_necorr_results(
        results,
        top_n = top_n,
        interactive_net = interactive_net,
        highlight_regulators = TRUE,
        output_dir = output_dir
      )
    }
  }

  update_progress(15) # Visualization
  # print("Visualization completed.")
  #### 20. Save results ####
  if (save_results) {
    if (is.null(output_dir)) stop("Please provide an 'output_dir' to save results")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(results$necorrHub_nodes, file.path(output_dir, "necorrHub_nodes.csv"))
    if (!is.null(results$necorrReg)) fwrite(results$necorrReg, file.path(output_dir, "necorrReg.csv"))
    if (!is.null(results$necorrEdges)) fwrite(results$necorrEdges, file.path(output_dir, "necorrEdges.csv"))
    fwrite(results$coexpres, file.path(output_dir, "coexpres.csv"))
    fwrite(results$hub.m.param, file.path(output_dir, "hub.m.param.csv"))
    fwrite(results$netstat, file.path(output_dir, "netstat.csv"))
  }
  update_progress(16) # Saving results
  # print("Results saved.")
  return(results)

}

