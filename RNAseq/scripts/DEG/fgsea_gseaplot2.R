# Workflow for fgsea analysis using DE results as input #

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
  library(enrichplot)
})

######Parsing input options and setting defaults########
# Load DE analyse from DESeq2
option_list<-list(
	make_option('--input',default='data', help='DE analysis file',dest='datafile'),
	make_option('--output',default='results', help='Folder were to save the results',dest='resultfolder'),
	make_option('--pathway',default='', help='Pathways to use, gmt file',dest='pathway')
)
opt<-parse_args(OptionParser(option_list=option_list))

FolderOutput=opt$resultfolder
DataInput=opt$datafile
Pathways=opt$pathway

dir.create(FolderOutput)

# Load gmt
gmt_file_name <- file_path_sans_ext(basename(Pathways))

# Load datas
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
res2

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
fgseaResTidy_filtered <- fgseaResTidy %>% filter(padj < 0.001)

ggplot(fgseaResTidy_filtered, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= paste("Results for GMT file: ", gmt_file_name)) + 
  theme_minimal()

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

# Plot function
plot_function <- function(gene_set_name) {
  gene_set <- filtered_pathways_hallmark[[gene_set_name]]
  
  # Create the enrichment plot using gseaplot2
  plot <- gseaplot2(x = fgseaRes, geneSetID = gene_set_name, title = paste("Enrichment Plot for: ", gene_set_name)) +
         geom_col(aes(fill = padj < 0.05)) +  # Color by p-value
         coord_flip()                         # Flip x and y axes

  # Customize the plot further (optional)
  # plot <- plot +
  #   scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey"))  # Customize color palette
  #   theme_minimal()  # Apply a minimal theme

  output_file <- paste0(FolderOutput, gsub(" ", "_", gene_set_name), ".pdf")
  
  # Save the plot
  ggsave(output_file, plot=plot, width=8, height=6)
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
