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
#' @param output_dir Optional, directory to save plots as PDF/PNG.
#'
#' @return A list of ggplot/ggraph/visNetwork objects.
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal
#' @importFrom ggplot2 geom_point geom_histogram scale_size_continuous scale_color_manual theme_void
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom igraph "V<-"
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions
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
                                     output_dir = NULL) {
  requireNamespace("ggraph")
  requireNamespace("igraph")
  requireNamespace("scales")
  requireNamespace("dplyr")
  requireNamespace("rlang")
  requireNamespace("patchwork") # For multi-panel layout

  if (interactive_net) requireNamespace("visNetwork")

  # --- 1. Top hub genes ---
  hub_df <- necorr_results$necorrHub_nodes
  hub_df$Score = as.numeric(hub_df$Score)
  hub_df <- as.data.frame(hub_df, stringsAsFactors = FALSE)
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

  # --- 2. Significant hub interaction network ---
  edges <- necorr_results$necorrEdges
  nodes <- data.frame(name = unique(c(edges$source, edges$target)))
  g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  igraph::V(g)$hub_score <- hub_df$Score[match(igraph::V(g)$name, hub_df$Gene)]

  if (highlight_regulators && "necorrReg" %in% names(necorr_results)) {
    regulators <- necorr_results$necorrReg$Gene
    igraph::V(g)$is_regulator <- ifelse(igraph::V(g)$name %in% regulators,
                                        "Regulator", "Other")
  } else {
    igraph::V(g)$is_regulator <- "Other"
  }

  if (!interactive_net) {
    reg_nodes <- igraph::V(g)$name[igraph::V(g)$is_regulator == "Regulator"]
    label_nodes <- head(reg_nodes, 20)

    p2 <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(alpha = 0.4) +
      ggraph::geom_node_point(ggplot2::aes(size = hub_score, color = is_regulator)) +
      ggplot2::scale_size_continuous(range = c(2, 10), name = "Hub Score") +
      ggplot2::scale_color_manual(values = c("Regulator" = "red", "Other" = "grey"),
                                  name = "Node Type") +
      ggraph::geom_node_text(
        ggplot2::aes(label = ifelse(name %in% label_nodes, name, "")),
        repel = TRUE, size = 3, max.overlaps = Inf
      ) +
      ggplot2::labs(title = "Significant Hub Interaction Network") +
      ggplot2::theme_void()
  } else {
    nodes_vis <- data.frame(
      id = igraph::V(g)$name,
      label = igraph::V(g)$name,
      value = igraph::V(g)$hub_score,
      group = igraph::V(g)$is_regulator,
      title = igraph::V(g)$name
    )
    edges_vis <- data.frame(from = edges$source, to = edges$target)
    p2 <- visNetwork::visNetwork(nodes_vis, edges_vis) %>%
      visNetwork::visNodes(scaling = list(min = 5, max = 30)) %>%
      visNetwork::visEdges(color = list(color = "#ccc", highlight = "red")) %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
  }

  # --- 3. Flexible Co-expression scatter plots ---
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

      # Filter valid rows only
      df_plot <- coexp_df %>%
        dplyr::filter(!is.na(.data[[corr_col]]),
                      !is.na(.data[[pval_col]]),
                      .data[[pval_col]] > 0)

      if (nrow(df_plot) > 0) {
        p <- ggplot2::ggplot(df_plot,
                             ggplot2::aes_string(x = corr_col,
                                                 y = paste0("-log10(", pval_col, ")"))) +
          ggplot2::geom_point(alpha = 0.4, color = "darkgreen") +
          ggplot2::labs(title = paste("Co-expression Significance -", corr),
                        x = corr_col, y = "-log10(p-value)") +
          ggplot2::theme_minimal(base_size = 14)
        p3_list[[corr]] <- p
      }
    }
  }

  # Combine p3_list into dynamic multipanel plot
  if (length(p3_list) > 0) {
    n_plots <- length(p3_list)
    ncol_layout <- ifelse(n_plots == 1, 1,
                          ifelse(n_plots <= 4, 2, 3))
    p3_combined <- Reduce(`+`, p3_list) + patchwork::plot_layout(ncol = ncol_layout)
  } else {
    p3_combined <- NULL
    n_plots <- 0
    ncol_layout <- 0
  }

  # --- 4. Degree distribution ---
  deg_df <- data.frame(Degree = igraph::degree(g))
  p4 <- ggplot2::ggplot(deg_df, ggplot2::aes(x = Degree)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "orange", color = "black") +
    ggplot2::labs(title = "Degree Distribution", x = "Degree", y = "Count") +
    ggplot2::theme_minimal(base_size = 14)

  # --- Print plots when function is called ---
  #print(p1)
  #if (!interactive_net) print(p2)
  #if (!is.null(p3_combined)) print(p3_combined)
  #print(p4)

  # --- Save plots ---
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(output_dir, "top_hubs_plot.pdf"), p1, width = 6, height = 4)
    if (!interactive_net) {
      ggplot2::ggsave(file.path(output_dir, "network_plot.pdf"), p2, width = 7, height = 6)
    }
    if (!is.null(p3_combined)) {
      # Auto-size: each plot ~5" wide and 4" tall
      nrow_layout <- ceiling(n_plots / ncol_layout)
      width_auto <- ncol_layout * 5
      height_auto <- nrow_layout * 4
      ggplot2::ggsave(file.path(output_dir, "coexpression_plots_combined.pdf"),
                      p3_combined, width = width_auto, height = height_auto)
    }
    if (length(p3_list) > 0) {
      for (nm in names(p3_list)) {
        ggplot2::ggsave(file.path(output_dir, paste0("coexpression_plot_", nm, ".pdf")),
                        p3_list[[nm]], width = 6, height = 4)
      }
    }
    ggplot2::ggsave(file.path(output_dir, "degree_distribution.pdf"), p4, width = 6, height = 4)
  }

  return(list(
    top_hubs_plot = p1,
    network_plot = p2,
    #coexpression_plots = p3_list,
    coexpression_combined = p3_combined,
    degree_distribution = p4
  ))
}
