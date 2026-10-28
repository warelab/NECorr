# NECorr: Network Expression Correlation Ranking

NECorr ranks the genes of a molecular network for their importance in one
condition (a tissue, a treatment, a developmental stage). It combines
expression data with the network: genes score high when they are co-expressed
with their network partners, specific to or differentially expressed in the
condition, and well connected.

NECorr prioritises candidates for experimental follow-up. It does not infer the
network: you supply one that was measured independently of the expression data.

## Installation

Requires R >= 4.0 and a C++ compiler (Rtools on Windows, Xcode command line
tools on macOS).

```r
install.packages("BiocManager")
BiocManager::install("limma")
install.packages("remotes")
remotes::install_github("warelab/NECorr")
```

`visNetwork` is optional and only needed for the interactive network plot.

## Quick start

```r
library(NECorr)
data(necorr_example)

res <- NECorr(
  networkFile      = necorr_example$networkFile,
  expression       = necorr_example$expression,
  description.file = necorr_example$description.file,
  meta_data        = necorr_example$meta_data,
  condition        = "flower",
  method           = "GCC",
  permutation      = 1000,
  save_results     = FALSE
)

head(res$necorrHub_nodes)   # ranked genes
head(res$necorrReg)         # candidate regulators of the top genes
head(res$necorrEdges)       # top-ranked edges
res$criterion_weights       # weights used
```

Each input can be a file path, a `data.frame` or a `data.table`.

| Argument | Content |
| --- | --- |
| `networkFile` | Edges: source in column 1, target in column 2 |
| `expression` | Gene IDs in column 1, then one column per sample (log2 scale) |
| `description.file` | Gene IDs in column 1, then any annotation columns |
| `meta_data` | Sample IDs in column 1 (or as row names), condition in the next column |
| `condition` | The condition to rank genes for |

The same example is shipped as text files:
`system.file("extdata", "example_network.txt", package = "NECorr")`, and
likewise `example_expression.txt`, `example_gene_description.txt` and
`example_sample_metadata.txt`.

## How genes are scored

By default the score is the rank product of two quantities:

| Criterion | Meaning |
| --- | --- |
| `DE` (up-regulation) | Moderated t-statistic of the condition against the mean of the other conditions (limma), set to zero when negative |
| `Degree` | Number of targets of the gene in a directed network; number of partners in an undirected network |

Each is converted to a rank divided by the number of genes, and the score is
the geometric mean of the two. There is no weight to choose. In a directed
network the up-regulated regulators with many targets come first.

This default was chosen on human cancer cell lines, where genome-wide CRISPR
screens give a complete ground truth. Counting only up-regulation
(`direction = "up"`) and combining it with the number of targets placed more
condition-specific dependencies at the top of the ranking than expression
alone, network propagation or any of the other settings below. Counting
down-regulation as well (`direction = "both"`) lowered the number of hits. The
same rule was then applied to rice tissues and Arabidopsis root cell types.

Four further criteria can be added with `criteria =`:

| Criterion | Meaning |
| --- | --- |
| `TS` | Selective expression in the condition: the moderated t-statistic that 80% of the pairwise comparisons with the other conditions reach |
| `Coexpr` | Number of network partners that are strongly co-expressed with the gene, compared with random partners drawn as in a degree-preserving randomisation of the network |
| `CondCoexpr` | As `Coexpr`, for the part of the correlation that comes from the samples of the condition |
| `Bridging` | Betweenness centrality beyond what the gene's degree predicts |

`direction = "both"` counts down-regulation and exclusion from the condition as
well. The six-criteria weighted sum of version 3.3 is obtained with:

```r
NECorr(..., criteria = c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree"),
       weights = "ahp", aggregation = "sum", direction = "both", scaling = "softcap",
       degree_type = "total")
```

`degree_type = "total"` counts regulators and targets together, as in earlier
versions:

```r
NECorr(..., directed = TRUE, degree_type = "total")
```

A ranking that uses only the wiring and the expression, and does not increase
with degree:

```r
c3 <- c("Coexpr", "TS", "CondCoexpr")
NECorr(..., criteria = c3, aggregation = "product", scaling = "softcap",
       weights = necorr_ahp_matrix(c3, preset = "wiring"))
```

### Weights and aggregation

| `weights =` | Method |
| --- | --- |
| `"equal"` (default) | Equal weights |
| `"ahp"` | Analytic Hierarchy Process, see `necorr_ahp_matrix()`. The consistency ratio is in `res$settings`. |
| `"critic"` | CRITIC: data-driven, down-weights criteria that duplicate each other |
| `"entropy"` | Entropy weight method |
| named numeric vector | Your own weights |
| pairwise matrix | Your own AHP judgements: edit `necorr_ahp_matrix()` and pass it |

`aggregation = "product"` (default) multiplies the weighted criteria, `"sum"`
adds them, `"waspas"` averages the two. `scaling = "rank"` (default) uses ranks,
`"softcap"` a bounded 0 to 1 scale.

## Network types

Set `directed = TRUE` for regulatory networks (transcription factor to target:
DAP-seq, ChIP-seq, yeast one-hybrid) and `directed = FALSE` for interaction
networks (yeast two-hybrid, AP-MS). With `directed = NULL` NECorr guesses and
tells you what it chose.

`res$necorrReg` tests whether a gene's partners are enriched in top-scoring
genes. In a directed network the partners are its targets and the genes are
reported as candidate regulators. In an undirected network direction is
unknown and they are reported as "associated".

Do not use a network built from the same expression data (co-expression
networks, or networks inferred with GENIE3, ARACNe or similar): the
co-expression criterion would re-measure how the edges were made.

## Significance

`res$necorrEdges` holds the top `edge_top` fraction of edges by score. The
`percentile` column is a rank, not a p-value.

With `n_null > 0`, expression profiles are shuffled among the network genes and
each gene gets an empirical p-value against genes of similar degree
(`p_empirical`, `FDR`). Use `n_null = 200` or more if you need FDR values.
Few genes pass an FDR threshold on large networks; the p-values are best used
to see where the ranking stops being distinguishable from chance.

## Functions

`NECorr()`, `multiCorr()`, `gene_coexpression()`, `coexpression_background()`,
`condition_specificity()`, `ts.IUT()`, `DE.ranking()`,
`rank_regulators()`, `hub_edge_percentile()`, `detect_direction()`,
`necorr_ahp_matrix()`, `ahp_weights()`, `critic_weight()`,
`sanitize_expression()`, `visualize_necorr_results()`. Each has a help page,
for example `?NECorr`. A tutorial is in `vignette("NECorr_tutorial")`.

## Citation

Liseron-Monfils C, Olson A, Ware D (2018). NECorr, a tool to rank gene
importance in biological processes using molecular networks and transcriptome
data. bioRxiv 326868. <https://doi.org/10.1101/326868>

The preprint describes version 1. Version 3 keeps the aim, ranking the genes
of a network for one condition, with a default validated against experimental
data and corrected statistics.

## License

MIT, see [LICENSE](./LICENSE).
