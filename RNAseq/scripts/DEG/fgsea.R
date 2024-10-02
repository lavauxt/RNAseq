# Workflow for fgsea analysis using DE results as input #

## Authoring : T LAVAUX
## Licence MIT

# Sources for reference
# https://stephenturner.github.io/deseq-to-fgsea/
# https://github.com/hamidghaedi/Enrichment-Analysis
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

# Guide for gmt sets
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
# https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/

message('Start fgsea script')

# Load library
suppressPackageStartupMessages({
  library(tools)
  library(tidyverse)
  library(htmlwidgets)
  library(ggplot2)
  library(optparse)
  library(fgsea)
  library(DT)
  library(enrichplot)
  library(clusterProfiler)
})

######Parsing input options and setting defaults########
# Load DE analysis from DESeq2
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = "input_file.tsv", 
              help = "Input file path for DE results", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "output_folder", 
              help = "Output folder for plots", metavar = "character"),
  make_option(c("-p", "--padj"), type = "numeric", default = 0.05, 
              help = "Adjusted p-value cutoff for filtering", metavar = "numeric"),
  make_option(c("-g", "--gmt"), type = "character", default = "hallmark.gmt", 
              help = "GMT file for pathways", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

output_folder <- opt$output
input_file <- opt$input
gmt_file_name  <- opt$gmt
padj_cutoff <- opt$padj

dir.create(output_folder, showWarnings = FALSE)

message('Extract gmt file name')
gmt_name <- file_path_sans_ext(basename(gmt_file_name))

message('Load DE results')
res <- read_tsv(paste(input_file))

message('Load the pathways for GSEA')
pathways.GSEA <- read.gmt(paste(gmt_file_name))

message('Load the pathways for fgsea')
pathways.fgsea <- gmtPathways(paste(gmt_file_name))

message('Convert the res data frame to a named vector for GSEA')
gene_list <- res$log2FoldChange
names(gene_list) <- res$gene

message('Remove NA values and sort the gene list in decreasing order')
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

message('Perform GSEA analysis')
gsea_results <- GSEA(gene_list, TERM2GENE = pathways.GSEA, pvalueCutoff = padj_cutoff)

message('Prepare data for fgsea')
res2 <- res %>%
  dplyr::select(gene, log2FoldChange) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(log2FoldChange = mean(log2FoldChange))
ranks <- deframe(res2)

message('Perform fgsea analysis')
fgseaRes <- fgseaMultilevel(pathways = pathways.fgsea, stats = ranks)

message('Save in an HTML table')
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
datatable_object <- fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()
output_file <- paste0("GSEA_Table_with_", gmt_name, ".html")
saveWidget(datatable_object, file = paste(output_folder, output_file, sep = "/"))

message('Plot the fgsea result into barplot')
# Filter significant pathways
fgseaResTidy_filtered <- fgseaResTidy %>% filter(padj < padj_cutoff)
# Replace underscores with spaces in pathway names
fgseaResTidy_filtered$pathway <- gsub("_", " ", fgseaResTidy_filtered$pathway)
# Wrap long pathway names to fit them into the plot
fgseaResTidy_filtered$pathway <- stringr::str_wrap(fgseaResTidy_filtered$pathway, width = 60)
# Plotting
plot <- ggplot(fgseaResTidy_filtered, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0), width = 0.6) +  # Adjust bar width if needed
  coord_flip() +  
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 8),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  # Set fill colors for positive and negative NES (less pale colors)
  scale_fill_manual(values = c("TRUE" = "#1E90FF", "FALSE" = "#FF6347")) +  # Dodger Blue and Tomato
  guides(fill = "none") # Remove the fill legend
# Save the plot as a PDF with adjusted dimensions (increase height for readability)
output_file <- file.path(output_folder, paste0("Barplot_GSEA_Table_with_", gmt_name, ".pdf"))
num_pathways <- nrow(fgseaResTidy_filtered)
width_adjustment <- max(8, num_pathways * 0.25)  
ggsave(output_file, plot = plot, width = width_adjustment, height = 8, device = "pdf")

message('Iterate over each pathway and save GSEA plots with gseaplot2')
gene_set_ids <- gsea_results@result$ID
for (i in gene_set_ids) {
    # Check for valid enrichment score in gseaplot2
    if (!is.null(gsea_results@result[i, ]) && !any(is.na(gsea_results@result[i, ]$ES))) {
      p <- tryCatch(
        {
          pathway_name <- gsub("_", " ", gsea_results@result[gsea_results@result$ID == i, "Description"])
          gseaplot2(gsea_results, geneSetID = i, title = as.character(pathway_name))
        },
        error = function(e) {
          message(paste("Error generating gseaplot2 for pathway:", i, "Skipping..."))
          return(NULL)
        }
      )
    }
    
  # Save plots if valid
  if (!is.null(p)) {
    # Save as PDF
    pdf_file <- paste0(output_folder, "/GSEA_plot_", gmt_name, "_", i, ".pdf")
    ggsave(filename = pdf_file, plot = p, device = "pdf", width = 8, height = 6)

    # Save as TIFF with compression and specified DPI
    tiff_file <- paste0(output_folder, "/GSEA_plot_", gmt_name, "_", i, ".tiff")
    ggsave(filename = tiff_file, plot = p, device = "tiff", width = 8, height = 6, 
          compression = "lzw", dpi = 600)
  } else {
    message(paste("Skipping pathway:", i, "due to invalid results or plotting error."))
  }
}

message('All done')
### Output the versions of all tools used in the analysis
sessionInfo()