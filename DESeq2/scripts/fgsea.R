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
  
  # Create the enrichment plot
  plot <- plotEnrichment(pathway = gene_set, stats = ranks) +
    ggtitle(paste("Enrichment Plot for: ", gene_set_name))
  
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

stop('End script')

#________ Heatmap Plot_____________#
# still not working as intended
# pathways with significant enrichment score
sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
sig.gen <- unique(na.omit(gene.in.pathway$gene[gene.in.pathway$pathway %in% sig.path]))

### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
h.dat <- dcast(gene.in.pathway[, c(1,2)], gene~pathway)
rownames(h.dat) <- h.dat$gene
#h.dat <- h.dat[, -1]

h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]

h.dat[!is.na(h.dat)] = as.numeric(1)
h.dat[is.na(h.dat)] = as.numeric(0)

row_names <- rownames(h.dat)

# keep those genes with 3  or more occurences
h.dat_numeric <- data.frame(lapply(h.dat, as.numeric))
row_sums <- rowSums(h.dat_numeric)
result <- data.frame(RowName = row_names, RowSum = row_sums)
#table(data.frame(rowSums(h.dat)))

# 1       2    3    4    5    6 
# 1604  282   65   11    1    1 
#h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]

filtered_genes <- result$RowName[result$RowSum >= 3]
h.dat <- h.dat[filtered_genes, ]

topTable <- res[res$gene %in% rownames(h.dat), ]
rownames(topTable) <- topTable$gene

# match the order of rownames in toptable with that of h.dat
topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
rownames(topTableAligned) <- topTableAligned$gene
topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
rownames(topTableAligned) <- topTableAligned$gene
all(rownames(topTableAligned) == rownames(h.dat))

# colour bar for -log10(adjusted p-value) for sig.genes
dfMinusLog10FDRGenes <- data.frame(-log10(
 topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'padj']))
dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0

# colour bar for fold changes for sigGenes
dfFoldChangeGenes <- data.frame(
 topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'log2FoldChange'])

# merge both
dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')

dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                        ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
colours <- list(
 'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
haGenes <- rowAnnotation(
 df = dfGeneAnno,
 col = colours,
 width = unit(1,'cm'),
 annotation_name_side = 'top')

# Now a separate colour bar for the GSEA enrichment padj. This will 
# also contain the enriched term names via annot_text()
# colour bar for enrichment score from fgsea results
dfEnrichment <- fgseaRes[, c("pathway", "NES")]
dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
dd <- dfEnrichment$pathway
dfEnrichment <- dfEnrichment[, -1]
if (length(dd) == nrow(dfEnrichment)) {
  dfEnrichment_new <- data.frame(dfEnrichment, row.names = dd)
} else {
  stop("Length of 'dd' does not match the number of rows in 'dfEnrichment'")
}

colnames(dfEnrichment_new) <- 'Normalized\n Enrichment score'

haTerms <- HeatmapAnnotation(
 df = dfEnrichment_new,
 Term = anno_text(
 colnames(h.dat),
 rot = 45,
 just = 'right',
 gp = gpar(fontsize = 10)),
 height = unit.c(unit(8, 'cm')),
 annotation_name_side = 'left')

# now generate the heatmap
hmapGSEA <- Heatmap(h.dat,
                   name = 'GSEA hallmark pathways enrichment',
                   split = dfGeneAnno[,2],
                   col = c('0' = 'white', '1' = 'forestgreen'),
                   rect_gp = gpar(col = 'grey85'),
                   cluster_rows = TRUE,
                   show_row_dend = TRUE,
                   row_title = 'Top Genes',
                   row_title_side = 'left',
                   row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                   row_title_rot = 90,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                   row_names_side = 'left',
                   row_dend_width = unit(35, 'mm'),
                   cluster_columns = TRUE,
                   show_column_dend = TRUE,
                   column_title = 'Enriched terms',
                   column_title_side = 'top',
                   column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                   column_title_rot = 0,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   clustering_distance_columns = 'euclidean',
                   clustering_method_columns = 'ward.D2',
                   clustering_distance_rows = 'euclidean',
                   clustering_method_rows = 'ward.D2',
                   bottom_annotation = haTerms)

tiff("name.tiff", units="in", width=13, height=33, res=400)
draw(hmapGSEA + haGenes,
    heatmap_legend_side = 'right',
    annotation_legend_side = 'right')
dev.off()

### Output the versions of all tools used in the analysis
sessionInfo()