

# --- Package stuff ------------------------------------------------------------

library(tidyverse)
requireNamespace("assertthat")
requireNamespace("DeNumerator")
requireNamespace("DESeq2")

# --- Load and clean the input data --------------------------------------------

sanitize_names <- function(x) {
  x %>% str_replace_all(" |-", "_")
}

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
  
  main_id <- colnames(raw_data)[1]
  gene_metadata <- raw_data[c(main_id, "SYMBOL", "GENENAME", "GENETYPE")]
  
  data <- raw_data %>% select(-all_of(c("SYMBOL", "GENENAME", "GENETYPE"))) %>%
    column_to_rownames(main_id)
  colnames(data) |> sanitize_SampleNames() -> colnames(data)
  
  sample_metadata <- read_tsv(
    metadata_path,
    skip = 1, comment="#"
  )
  
  # IMPORTANT
  # ---------
  # By default, R will choose a reference level for factors based on
  # alphabetical order. To set the reference group, you have to explicitly set
  # the factor levels, keeping in mind that the first one will be taken as
  # the reference (denominator in the Fold Change computation), by using
  # factor(data, levels = c("untreated","treated"))
  # or
  # relevel(factor, ref = "untreated")
  sample_metadata$treatment <- factor(
    ifelse(grepl("Control", sample_metadata$sample_alias), "Ctrl",
           ifelse(grepl("Acute", sample_metadata$sample_alias), "Acute", "Selected")),
    levels = c("Ctrl", "Acute", "Selected")
  )
  sample_metadata <- sample_metadata %>% column_to_rownames("file_id")
  
  # # Filter out batched samples (Replicate 1)
  # old_cols <- ncol(data)
  # data <- data %>% select(! c("4-days-pH-6_6-7_5", "pH-selected-17"))
  # sample_metadata <- sample_metadata %>% filter(row.names(sample_metadata) %in% colnames(data))
  # cat(paste0("Deleted ", old_cols - ncol(data), " samples!\n"))
  
  list(
    counts = data,
    metadata = list(
      gene = gene_metadata,
      sample = sample_metadata
    )
  )
}

gene_data <- load_data(data_path = "./data/in/pH_CountMatrix_genes_TPM.tsv")
isoform_data <- load_data(data_path = "./data/in/pH_CountMatrix_isoforms_TPM.tsv")

# Check the colnames/rownames order
data$counts <- data$counts[, rownames(data$metadata$sample)]
assertthat::are_equal(colnames(data$counts), rownames(data$metadata$sample))

# Fibroblast Growth Factor Receptor 2
# HGNC: 3689 NCBI Gene: 2263 Ensembl: ENSG00000066468 OMIM: 176943 UniProtKB/Swiss-Prot: P21802


gene_data$metadata$gene[gene_data$metadata$gene$ENSEMBL == "ENSG00000066468",]
gene_data$counts["ENSG00000066468",]




all_transcripts <- read.csv("./data/in/all_FGFR2_transcripts.txt",
                            header= FALSE) |> unlist()



isoform_counts <- isoform_data$counts[rownames(isoform_data$counts) %in% all_transcripts,]



data.frame(gene_level = as.numeric(gene_data$counts["ENSG00000066468",]),
           isoform_sum = colSums(isoform_counts))


mean_level <- function(isoform_dataset, transcript, group) {
  isoform_dataset$metadata$sample |> filter(treatment == group) |>
    rownames() -> samples
  isoform_dataset$counts[transcript, samples] |> rowMeans()
  }

mean_level(isoform_data, "ENST00000358487", "Ctrl")
mean_level(isoform_data, "ENST00000457416", "Ctrl")

mean_level(isoform_data, "ENST00000358487", "Acute")
mean_level(isoform_data, "ENST00000457416", "Acute")

mean_level(isoform_data, "ENST00000358487", "Selected")
mean_level(isoform_data, "ENST00000457416", "Selected")





