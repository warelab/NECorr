if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "AnnotationDbi",
  "GO.db",
  "preprocessCore",
  "Biobase",
  "limma",
  "impute"
), update = FALSE, ask = FALSE)
