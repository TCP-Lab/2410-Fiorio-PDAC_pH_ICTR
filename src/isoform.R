

# --- Package stuff ------------------------------------------------------------

library(tidyverse)
requireNamespace("assertthat")
requireNamespace("DeNumerator")
requireNamespace("DESeq2")

# --- Load and clean the input data --------------------------------------------

sanitize_SampleNames <- function(x) {
  x |> str_remove_all("PANC1-?") |>
    str_remove_all("(expected_count|TPM)(\\..*)?") |>
    str_remove_all("_R_")
}

# Metadata path is default for this project
load_data <- function(data_path,
                      metadata_path = "./data/in/sample_metadata.tsv")
{
  raw_data <- read_tsv(data_path)
  
  # Good for both genes and isoforms
  main_id <- colnames(raw_data)[1]
  gene_metadata <- raw_data[c(main_id, "SYMBOL", "GENENAME", "GENETYPE")]
  
  data <- raw_data %>% select(-all_of(c("SYMBOL", "GENENAME", "GENETYPE"))) %>%
    column_to_rownames(main_id)
  colnames(data) |> sanitize_SampleNames() -> colnames(data)
  
  sample_metadata <- read_tsv(
    metadata_path,
    skip = 1, comment="#"
  )
  
  sample_metadata$treatment <- factor(
    ifelse(grepl("Control", sample_metadata$sample_alias), "Ctrl",
           ifelse(grepl("Acute", sample_metadata$sample_alias), "Acute", "Selected")),
    levels = c("Ctrl", "Acute", "Selected")
  )
  sample_metadata <- sample_metadata %>% column_to_rownames("file_id")
  
  # Filter out batched samples
  old_cols <- ncol(data)
  data <- data %>% select(! c("4-days-pH-6_6-7_5", "pH-selected-17"))
  sample_metadata <- sample_metadata %>% filter(row.names(sample_metadata) %in% colnames(data))
  cat(paste0("Deleted ", old_cols - ncol(data), " samples!\n"))
  
  list(
    counts = data,
    metadata = list(
      gene = gene_metadata,
      sample = sample_metadata
    )
  )
}

# # TPMs
# gene_data <- load_data(data_path = "./data/in/pH_CountMatrix_genes_TPM.tsv")
# isoform_data <- load_data(data_path = "./data/in/pH_CountMatrix_isoforms_TPM.tsv")

# Expected Counts
gene_data <- load_data(data_path = "./data/in/pH_CountMatrix_genes_expected_count.tsv")
isoform_data <- load_data(data_path = "./data/in/pH_CountMatrix_isoforms_expected_count.tsv")

# --- Count check --------------------------------------------------------------

# Gene of Interest
# Fibroblast Growth Factor Receptor 2
# HGNC: 3689 NCBI Gene: 2263 Ensembl: ENSG00000066468 OMIM: 176943 UniProtKB/Swiss-Prot: P21802

# Retrieve gene info and isoform counts
cat("\n", rep("-", 90), "\n", sep = "")
gene_data$metadata$gene[gene_data$metadata$gene$ENSEMBL == "ENSG00000066468",] |> print()
cat(rep("-", 90), "\n", sep = "")

all_isoforms <- read.csv("./data/in/all_FGFR2_transcripts.txt",
                         header= FALSE) |> unlist()
isoform_counts <- isoform_data$counts[all_isoforms,]
# # Show all isoform counts
# isoform_counts |> print()
cat("\n")
data.frame(gene_level = as.numeric(gene_data$counts["ENSG00000066468",]),
           isoform_sum = colSums(isoform_counts)) |> print()
cat("\n")

# --- Isoform stats ------------------------------------------------------------

# Get stats for the transcripts of interest within a given experimental group
mean_level <- function(isoform_dataset, transcripts, group) {
  isoform_dataset$metadata$sample |> filter(treatment == group) |>
    rownames() -> samples
  isoform_dataset$counts[transcripts, samples] -> counts
  data.frame(mean = rowMeans(counts), SD = apply(counts, 1, sd)) |>
      round(digits = 4) -> df
  
  cat("\n", ">>> ", group, "\n", sep = "")
  return(df)
}

# Isoforms of Interest
ioi <- c("ENST00000358487", "ENST00000457416")
mean_level(isoform_data, ioi, "Ctrl")
mean_level(isoform_data, ioi, "Acute")
mean_level(isoform_data, ioi, "Selected")

cat("\n")



