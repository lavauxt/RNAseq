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
  library(tidyverse)
  library(reshape2)
  library(dplyr)
  library(tools)
  library(htmlwidgets)
  library(maditr)
  library(ggplot2)
  library(optparse)

  library(fgsea)
  library(DT)
  library(ComplexHeatmap)
  library(DOSE)
  library(clusterProfiler)
  library(enrichplot)
})

######Parsing input options and setting defaults########
# Load DE analysis from DESeq2
option_list <- list(
  make_option('--input', default = 'data', help = 'DE analysis file', dest = 'datafile'),
  make_option('--output', default = 'results', help = 'Folder where to save the results', dest = 'resultfolder'),
  make_option('--pathway', default = '', help = 'Pathways to use, gmt file', dest = 'pathway')
)
opt <- parse_args(OptionParser(option_list = option_list))

FolderOutput <- opt$resultfolder
DataInput <- opt$datafile
Pathways <- opt$pathway

dir.create(FolderOutput, showWarnings = FALSE)

message('Load gmt')
gmt_file_name <- file_path_sans_ext(basename(Pathways))

message('Load DE results')
res <- read_tsv(paste(DataInput))

message('Load the pathways for GSEA')
pathways.GSEA <- read.gmt(paste(Pathways))

message('Load the pathways for fgsea')
pathways.fgsea <- gmtPathways(paste(Pathways))

message('Convert the res data frame to a named vector for GSEA')
gene_list <- res$log2FoldChange
names(gene_list) <- res$gene

message('Remove NA values and sort the gene list in decreasing order')
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

message('Perform GSEA analysis')
gsea_results <- GSEA(gene_list, TERM2GENE = pathways.GSEA, pvalueCutoff = 0.05)

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

message('Show in an HTML table')
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
datatable_object <- fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()
output_file <- paste0("GSEA_Table_with_", gmt_file_name, ".html")
saveWidget(datatable_object, file = paste(FolderOutput, output_file, sep = "/"))

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
    pdf_file <- paste0(FolderOutput, "/GSEA_plot_", gmt_file_name, "_", i, ".pdf")
    ggsave(filename = pdf_file, plot = p, width = 8, height = 6)
    
    # Save as TIFF
    tiff_file <- paste0(FolderOutput, "/GSEA_plot_", gmt_file_name, "_", i, ".tiff")
    ggsave(filename = tiff_file, plot = p, width = 8, height = 6, device = "tiff")
  } else {
    message(paste("Skipping pathway:", i, "due to invalid results or plotting error."))
  }
}

message('All done')
### Output the versions of all tools used in the analysis
sessionInfo()