# Workflow for fgsea analysis using DE results as input

## Authoring : T LAVAUX
## Licence MIT

# Sources for reference
# https://stephenturner.github.io/deseq-to-fgsea/
# https://github.com/hamidghaedi/Enrichment-Analysis

# Guide for gmt sets
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
# https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/

suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(dplyr)
  library(tools)
  library(htmlwidgets)
  library(maditr)
  library(ggplot2)
  library(fgsea)
  library(DT)
  library(ComplexHeatmap)
  library(optparse)
  library(patchwork) 
})

###### Parsing input options and setting defaults ########
# Load DE analyze from DESeq2
option_list <- list(
  make_option('--input', default='data', help='DE analysis file', dest='datafile'),
  make_option('--output', default='results', help='Folder where to save the results', dest='resultfolder'),
  make_option('--pathway', default='', help='Pathways to use, gmt file', dest='pathway')
)
opt <- parse_args(OptionParser(option_list=option_list))

FolderOutput <- opt$resultfolder
DataInput <- opt$datafile
Pathways <- opt$pathway

dir.create(FolderOutput)

# Load gmt
gmt_file_name <- file_path_sans_ext(basename(Pathways))

# Load data
res <- read_tsv(paste(DataInput))

# Load the pathways into a named list
pathways.hallmark <- gmtPathways(paste(Pathways))

# We use log fold change cf https://support.bioconductor.org/p/129277/
res2 <- res %>%
  dplyr::select(gene, log2FoldChange) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(log2FoldChange=mean(log2FoldChange))

ranks <- deframe(res2)
head(ranks, 20)

# Perform fgsea analysis
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a html table
datatable_object <- fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

# Save in html
output_file <- paste("GSEA_Table_with", gmt_file_name, ".html")
saveWidget(datatable_object, file = paste(FolderOutput, output_file))

# Filtered results by padj value
fgseaResTidy_filtered <- fgseaResTidy %>% filter(padj < 0.05)

# Create the ggplot with the desired formatting
plot <- ggplot(fgseaResTidy_filtered, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +          # Fill by padj < 0.05
  coord_flip() +                               # Flip coordinates for horizontal bars
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = paste("Results for GMT file: ", gmt_file_name)) +
  theme_minimal() +
  
  # Bold the text elements
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm")  # Optional margin adjustment
  ) +
  
  # Remove the color legend
  guides(fill = "none")

# Save the plot as a PDF with adjusted dimensions (thinner plot)
output_file <- file.path(FolderOutput, paste0("GSEA_Table_with_", gmt_file_name, ".pdf"))
ggsave(output_file, plot = plot, width = 6, height = 8, device = "pdf")

# Notify user of success
message("Plot saved to: ", output_file)

# Table by gene
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "gene") %>%
  unnest(cols = c(gene)) %>%
  inner_join(res, by="gene")

# Filter pathways based on padj < 0.001
filtered_pathways <- fgseaResTidy %>%
  filter(padj < 0.001) %>%
  pull(pathway)

# Subset the pathways.hallmark to include only the filtered pathways
filtered_pathways_hallmark <- pathways.hallmark[filtered_pathways]

# Plot function for enrichment plots with ranked list metrics
plot_function <- function(gene_set_name) {
  gene_set <- filtered_pathways_hallmark[[gene_set_name]]
  
  # Create the enrichment plot
  enrichment_plot <- plotEnrichment(pathway = gene_set, stats = ranks) +
    ggtitle(paste("Enrichment Plot for: ", gene_set_name)) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    ) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Ranked Genes", y = "Enrichment Score")
  
  # Create the ranked list metric plot
  metric_plot <- ggplot() +
    geom_segment(aes(x = which(ranks %in% gene_set), xend = which(ranks %in% gene_set), 
                     y = 0, yend = 1), color = "black") +   # Vertical lines for ranked genes
    theme_void() +
    labs(x = "Ranked Genes", y = "") +
    theme(axis.title.x = element_text(size = 12, face = "bold"),
          axis.ticks.x = element_blank())  # No x-axis ticks for metrics
  
  # Combine the enrichment plot and ranked list metric plot
  combined_plot <- enrichment_plot / metric_plot   # Stack plots vertically
  
  # Save the combined plot
  output_file <- paste0(FolderOutput, gsub(" ", "_", gene_set_name), ".pdf")
  ggsave(output_file, plot = combined_plot, width = 8, height = 6)
}

# Apply the function to each gene set in pathways.hallmark
lapply(names(filtered_pathways_hallmark), plot_function)

output_file <- paste("GSEA_Table_with", gmt_file_name, ".pdf")
pdf(paste(FolderOutput, output_file))
plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.001]], ranks, fgseaRes, 
               gseaParam=0.5)
dev.off()

### Output the versions of all tools used in the analysis
sessionInfo()
