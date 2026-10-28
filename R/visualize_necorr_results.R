#' visualize_necorr_results
#' @description Visualize NECorr Analysis Results Generates publication-ready
#' plots from NECorr output, including
#' top hub genes, significant hub interaction network, co-expression
#' significance scatter, and degree distribution.
#'
#' @param necorr_results List output from NECorr().
#' @param top_n Integer, number of top hub genes to plot.
#' @param interactive_net Logical, whether to produce an interactive
#' network plot (visNetwork).
#' @param highlight_regulators Logical, whether to highlight regulators in network plot.
#' @param max_network_edges maximum edge that is shown in the graph
#' @param output_dir Optional, directory to save plots as PDF/PNG.
#'
#' @return A list of ggplot/ggraph/visNetwork objects.
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal
#' @importFrom ggplot2 geom_point geom_histogram scale_size_continuous scale_color_manual theme_void
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom igraph "V<-"
#' @importFrom dplyr arrange filter mutate select
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' res <- NECorr(networkFile, expression, description.file, condition, metadata)
#' plots <- visualize_necorr_results(res, top_n = 20, interactive_net = FALSE)
#' plots$top_hubs_plot
#' }
visualize_necorr_results <- function(necorr_results,
                                     top_n = 20,
                                     interactive_net = FALSE,
                                     highlight_regulators = TRUE,
                                     output_dir = NULL,
                                     max_network_edges = 100) {
  requireNamespace("igraph")
  requireNamespace("scales")
  requireNamespace("dplyr")
  requireNamespace("rlang")
  requireNamespace("patchwork")
  has_vis <- requireNamespace("visNetwork", quietly = TRUE)
  requireNamespace("htmlwidgets")

  # --- 1. Top hub genes ---
  hub_df <- necorr_results$necorrHub_nodes
  hub_df$Score <- as.numeric(hub_df$Score)
  hub_df <- dplyr::arrange(hub_df, dplyr::desc(Score))
  hub_df <- dplyr::filter(hub_df, !is.na(Score))
  if (nrow(hub_df) < top_n) top_n <- nrow(hub_df)
  hub_df_top <- head(hub_df, top_n)

  p1 <- ggplot2::ggplot(hub_df_top,
                        ggplot2::aes(x = reorder(Gene, Score), y = Score)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(title = paste("Top", top_n, "Hub Genes (NECorr)"),
                  x = "Gene", y = "Hub Score") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(hub_df_top$Score) * 1.05))

  # --- 2. Filter edges early ---
  edges_df <- necorr_results$necorrEdges
  edges_df <- edges_df[order(-edges_df$edge_score), , drop = FALSE]
  edges_df <- head(edges_df, max_network_edges)

  # Build igraph from filtered edges
  nodes <- data.frame(name = unique(c(edges_df$sourceIDs, edges_df$targetIDs)))
  g <- igraph::graph_from_data_frame(edges_df, vertices = nodes, directed = FALSE)
  igraph::V(g)$hub_score <- hub_df$Score[match(igraph::V(g)$name, hub_df$Gene)]
  igraph::E(g)$weight <- edges_df$edge_score

  reg <- necorr_results$necorrReg
  if (highlight_regulators && !is.null(reg) && nrow(reg) > 0) {
    regulators <- reg$Gene[reg$FDR < 0.05]
    igraph::V(g)$is_regulator <- ifelse(igraph::V(g)$name %in% regulators,
                                        "Regulator", "Other")
  } else {
    igraph::V(g)$is_regulator <- "Other"
  }

  # --- 3. Prepare visNetwork data ---
  nodes_vis <- data.frame(
    id = igraph::V(g)$name,
    label = igraph::V(g)$name,
    value = igraph::V(g)$hub_score,  # node size
    group = igraph::V(g)$is_regulator,
    title = paste("Hub Score:", round(igraph::V(g)$hub_score, 2))
  )

  edges_vis <- data.frame(
    from = igraph::as_edgelist(g)[,1],
    to = igraph::as_edgelist(g)[,2],
    width = scales::rescale(igraph::E(g)$weight, to = c(1, 8)),  # thickness
    color = "#999999",
    title = paste("Edge score:", round(igraph::E(g)$weight, 2))
  )

  # --- 4. Network plot: interactive (visNetwork) when available, static otherwise ---
  p2 <- NULL
  if (has_vis) {
    p2 <- visNetwork::visNetwork(nodes_vis, edges_vis)
    p2 <- visNetwork::visNodes(p2, scaling = list(min = 5, max = 30),
                               color = list(background = "grey", border = "black", highlight = "yellow"))
    p2 <- visNetwork::visGroups(p2, groupname = "Regulator", color = list(background = "red", border = "black"))
    p2 <- visNetwork::visGroups(p2, groupname = "Other", color = list(background = "grey", border = "black"))
    p2 <- visNetwork::visEdges(p2, smooth = FALSE)
    p2 <- visNetwork::visOptions(p2, highlightNearest = TRUE, nodesIdSelection = TRUE)
    p2 <- visNetwork::visLayout(p2, randomSeed = 42)
    if (!interactive_net && !is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      htmlwidgets::saveWidget(p2, file = file.path(output_dir, "network_plot.html"), selfcontained = TRUE)
    }
  }
  if (!has_vis || !interactive_net) {
    p2_static <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_width = weight), colour = "grey70", alpha = 0.7) +
      ggraph::geom_node_point(ggplot2::aes(size = hub_score, colour = is_regulator)) +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
      ggraph::scale_edge_width(range = c(0.3, 2), guide = "none") +
      ggplot2::scale_colour_manual(values = c(Regulator = "firebrick", Other = "grey40"), name = NULL) +
      ggplot2::scale_size_continuous(name = "NECorr score") +
      ggplot2::theme_void()
    if (!interactive_net) p2 <- p2_static
  }

  # --- 6. Co-expression scatter plots ---
  coexp_df <- necorr_results$coexpres
  corr_types <- list(
    GCC = "GCC_pvalue",
    PCC = "PCC_pvalue",
    SCC = "SCC_pvalue",
    KCC = "KCC_pvalue"
  )
  p3_list <- list()
  for (corr in names(corr_types)) {
    corr_col <- corr
    pval_col <- corr_types[[corr]]
    if (all(c(corr_col, pval_col) %in% colnames(coexp_df))) {
      df_plot <- coexp_df[!is.na(coexp_df[[corr_col]]) & !is.na(coexp_df[[pval_col]]) &
                            coexp_df[[pval_col]] > 0, , drop = FALSE]
      if (nrow(df_plot) > 0) {
        p <- ggplot2::ggplot(df_plot,
                             ggplot2::aes(x = .data[[corr_col]],
                                          y = -log10(.data[[pval_col]]))) +
          ggplot2::geom_point(alpha = 0.4, color = "darkgreen") +
          ggplot2::labs(title = paste("Co-expression Significance -", corr),
                        x = corr_col, y = "-log10(p-value)") +
          ggplot2::theme_minimal(base_size = 14)
        p3_list[[corr]] <- p
      }
    }
  }
  p3_combined <- if (length(p3_list) > 0) {
    n_plots <- length(p3_list)
    ncol_layout <- ifelse(n_plots == 1, 1,
                          ifelse(n_plots <= 4, 2, 3))
    Reduce(`+`, p3_list) + patchwork::plot_layout(ncol = ncol_layout)
  } else NULL

  # --- 7. Degree distribution ---
  edges_df <- necorr_results$necorrEdges
  # Build igraph from filtered edges
  nodes <- data.frame(name = unique(c(edges_df$sourceIDs, edges_df$targetIDs)))
  g <- igraph::graph_from_data_frame(edges_df, vertices = nodes, directed = FALSE)
  deg_values <- igraph::degree(g)
  deg_df <- as.data.frame(table(deg_values))
  colnames(deg_df) <- c("Degree", "Count")
  deg_df$Degree <- as.numeric(as.character(deg_df$Degree))

  p4 <- ggplot2::ggplot(deg_df, ggplot2::aes(x = Degree, y = Count)) +
    ggplot2::geom_point(color = "orange", size = 3) +
    ggplot2::geom_line(color = "orange") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Degree Distribution (Log Scaled)",
                  x = "Degree (log scale)",
                  y = "Number of Nodes (log scale)") +
    ggplot2::theme_minimal(base_size = 14)

  return(list(
    top_hubs_plot = p1,
    network_plot = p2,
    coexpression_combined = p3_combined,
    degree_distribution = p4
  ))
}
