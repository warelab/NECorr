# Small simulated data set with one planted, condition-specific module
sim_data <- function(seed = 1, G = 200, directed = TRUE, effect = 3, labels = paste0("c", 1:4)) {
  set.seed(seed)
  genes <- paste0("G", seq_len(G))
  g <- igraph::sample_pa(G, m = 2, directed = TRUE)
  el <- igraph::as_edgelist(g)
  net <- data.frame(source = genes[el[, 2]], target = genes[el[, 1]], stringsAsFactors = FALSE)
  if (!directed) {  # arbitrary orientation, as in an interaction list
    flip <- stats::runif(nrow(net)) < 0.5
    net[flip, ] <- net[flip, 2:1]
  }
  outdeg <- table(c(net$source, net$target))
  src_deg <- table(net$source)
  cand <- names(src_deg)[src_deg >= 5 & src_deg <= 12]
  hub <- if (directed) cand[1] else names(outdeg)[outdeg >= 5 & outdeg <= 12][1]
  module <- unique(c(hub, net$target[net$source == hub], net$source[net$target == hub]))
  cond <- rep(labels, each = 4)
  E <- matrix(stats::rnorm(G * 16, 7, 1), G, 16, dimnames = list(genes, paste0("S", 1:16)))
  prof <- stats::rnorm(16); prof[cond == labels[1]] <- prof[cond == labels[1]] + effect
  E[hub, ] <- 7 + prof + stats::rnorm(16, 0, 0.3)
  for (v in setdiff(module, hub)) E[v, ] <- 7 + 0.8 * prof + stats::rnorm(16, 0, 0.6)
  list(net = net, expr = data.frame(GeneID = genes, E, check.names = FALSE),
       desc = data.frame(gene_id = genes, Description = paste("desc", genes)),
       meta = data.frame(SampleID = colnames(E), Condition = cond, stringsAsFactors = FALSE),
       hub = hub, module = module, target = labels[1])
}
run_necorr <- function(d, ...) {
  args <- utils::modifyList(list(networkFile = d$net, expression = d$expr, description.file = d$desc,
                                 condition = d$target, meta_data = d$meta, permutation = 0, method = "PCC",
                                 visualize = FALSE, save_results = FALSE, verbose = FALSE,
                                 # the six-criteria scheme of NECorr 3.3, which most tests exercise
                                 criteria = c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree"), weights = "ahp",
                                 aggregation = "sum", direction = "both", scaling = "softcap", degree_type = "total"), list(...))
  do.call(NECorr, args)
}
