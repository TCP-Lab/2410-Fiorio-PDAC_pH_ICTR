

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

# Paths are default for this project
load_data <- function(
        data_path = "./data/in/pH_CountMatrix_genes_expected_count.tsv",
        metadata_path = "./data/in/sample_metadata.tsv"
    ){
    raw_data <- read_tsv(data_path)
    gene_metadata <- raw_data[c("ENSEMBL", "SYMBOL", "GENENAME", "GENETYPE")]
    data <- raw_data %>% select(-all_of(c("SYMBOL", "GENENAME", "GENETYPE"))) %>%
        column_to_rownames("ENSEMBL")
    colnames(data) |> sanitize_SampleNames() -> colnames(data)
    
    # Pre-filtering
    # Ref: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#Pre-filtering
    # As per Ref. above, pre-filtering has no (or minimal) effect on DESeq but
    # it makes the object smaller and the computation (slightly) faster.
    # Remove zero-counts only 
    #data <- data[rowSums(data) > 0,]
    # keep only rows that have a count of at least 10, in at least 3 samples.
    smallest_group_size <- 3
    data <- data[rowSums(data >= 10) >= smallest_group_size,]
    
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

data <- load_data()

# Check the colnames/rownames order since DeSeq2 is dumb
data$counts <- data$counts[, rownames(data$metadata$sample)]
assertthat::are_equal(colnames(data$counts), rownames(data$metadata$sample))

# --- Run DESeq ----------------------------------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(data$counts),
                                      colData = data$metadata$sample,
                                      design = ~treatment)
dds_res <- DESeq2::DESeq(dds)

# --- Routine plots ------------------------------------------------------------

alpha_level <- 0.05

if (! dir.exists("data/out")) {
    dir.create("data/out")
}
pdf(file = file.path("data", "out", "MAplot.pdf"),
    width = 10, height = 6)
DESeq2::plotMA(dds_res, main = "MA Plot", alpha = alpha_level)
dev.off()

# --- Enumerate the results ----------------------------------------------------

enum <- DeNumerator::denumerate(dds_res, alpha = alpha_level, fold_change = 0)
pdf(file = file.path("data", "out", "Denumerator_Plot.pdf"),
    width = 10, height = 6)
DeNumerator::plot_enumeration_frame(enum)
dev.off()

# --- Save gene lists as annotated CSV -----------------------------------------

betas <- DESeq2::resultsNames(dds_res)[-1]
for (beta in betas) {
    DESeq2::results(dds_res, name = beta) |>
        as.data.frame() |> filter(padj < alpha_level) -> DEGs
    rich_DEGs <- merge(data$metadata$gene, DEGs,
                       by.x = "ENSEMBL", by.y = 0, all.y = TRUE)
    write.csv(rich_DEGs, file = file.path("data", "out",
                                          paste0("DEGs", beta, ".csv")))
}




