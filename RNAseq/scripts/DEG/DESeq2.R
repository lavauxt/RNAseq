# Workflow for DE analysis using DESeq2 #
# Import quantification via any tools (salmon/alevin is default tools)

## Original author : L. RIGOLOT
## Automomatisation and maintenance : T LAVAUX

# 15/09/2022 : add argparse options and refactor some code - TL
# 01/12/2022 : add KEGG patwhay graphs - TL
# 08/12/2022 : add Ensembl update package (v107), improve files naming & graphs (GO) - TL
# 19/01/2022 : add covar option & update SPIA PlotP to be more robust ; refactor some code - TL
# 02/02/2023 : add interaction option for the statistic model, sample table filter samples list from files ; some variabilisation & options added - TL

# Sources for reference
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html
# http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html

# Multifactor Designs in DESeq2 check https://www.youtube.com/watch?v=X6p3E-QTcUc

# In DTE, differential expression between conditions is assessed at the individual transcript level, while in DTU the relative expression of the isoforms of a gene are compared between conditions; i.e. a DTU analysis aims at discovering differences in the proportions of the expressed isoforms of a gene

# DTU/DTE http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html

# IsoformSwitchAnalyzeR https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html

# https://www.bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html

# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(optparse)
library(stringr)
library(tximport)
library(tidyverse)
library(gdata)
})

######Parsing input options and setting defaults########
option_list <- list(
  make_option(c("--input"), type="character", default='data', help='Folder where the data are stored [default %default]', dest='FolderOutput'),
  make_option(c("--output"), type="character", default='results', help='Folder where to save the results [default %default]', dest='DataInput'),
  make_option(c("--count"), type="character", default='salmon', help='Type of RNA count [default %default]', dest='CountType'),
  make_option(c("--sampletable"), type="character", default='./sample_table.csv', help='Path of the sample table [default %default]', dest='SampleTable'),
  make_option(c("--compare"), type="character", default='condition', help='Comparison to be analysed [default %default]', dest='Condition2Compare'),
  make_option(c("--level2compare"), type="character", default='level_to_compare', help='Level to compare with [default %default]', dest='level_to_compare'),
  make_option(c("--baselevel"), type="character", default='base_level', help='Base level [default %default]', dest='base_level'),
  make_option(c("--annotdb"), type="character", default='./database/annotations_ahb.csv', help='Annotation database [default %default]', dest='annotdbfilepath'),
  make_option(c("--gene"), type="character", default='CD180', help='Gene name for plot construction [default %default]', dest='geneplot'),
  make_option(c("--topplot"), type="integer", default=30, help='Number of top genes for top plots [default %default]', dest='topgene'),
  make_option(c("--topheatmap"), type="integer", default=50, help='Number of top genes for interactive heatmap [default %default]', dest='topheatgene'),
  make_option(c("--covar"), type="character", default=NULL, help='Covariable for the statistic model (ex batch, cell) [default %default]', dest='covariable'),
  make_option(c("--inter"), type="character", default=NULL, help='Interaction for the statistic model, separated by a : (ex condition:cell) [default %default]', dest='interaction'),
  make_option(c("--shrink"), type="character", default='ashr', help='Shrinking algorithm for DE analysis: normal, ashr (default), or apeglm [default %default]', dest='shrink_method'),
  make_option(c("--ensembl"), type="character", default='./database/EnsDb.Hsapiens.v107.tar.gz', help='Ensembl package file in .tar.gz format [default %default]', dest='ensembl_package'),
)

# Create the OptionParser and parse the arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Directly assign the options to variables
FolderOutput = opt$FolderOutput
DataInput = opt$DataInput
CountType = opt$CountType
SampleTable = opt$SampleTable
Condition2Compare = opt$Condition2Compare
level_to_compare = opt$level_to_compare
base_level = opt$base_level
annotdbfilepath = opt$annotdbfilepath
geneplot = opt$geneplot
topgene = opt$topgene
topheatgene = opt$topheatgene
covariable = opt$covariable
interaction = opt$interaction
shrink_method = opt$shrink_method
ensembl_package = opt$ensembl_package

### STEP 1 ### Setup to import datas into proper files
# Create directories
dir.create(FolderOutput)
subdirs <- c("Counts", "FoldChanges", "Plots", "KEGG", "GSEA", "SPIA", "QC", "Databases", "Log")
sapply(subdirs, function(subdir) dir.create(file.path(FolderOutput, subdir)))

# Extracting samples from sample_table
sample_table <- read.csv(paste(SampleTable), sep=";")
write.table((sample_table), paste(FolderOutput ,"/Log/Copy of the sample table.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
samples_to_analyse <- data.frame(sample_table$sample)
write.table((samples_to_analyse), paste(FolderOutput ,"/Log/List of samples from the sample table.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

# Extracting samples name from files
files_name <- list.files(path = DataInput, full.names = FALSE)
# Read samples name from thoses files name 
samples_files <- gsub("\\..*","",files_name) # sample name is the first part of the name before the first dot
# Extracting full path of the sample files count
fullpathfiles <- list.files(path = DataInput, full.names = TRUE)
# Creating table for tximport : name of the sample and full path of the file
names(fullpathfiles) <- samples_files
write.table((fullpathfiles), paste(FolderOutput ,"/Log/Tximport_Unfiltered list.txt", sep = ""), sep="\t", quote=FALSE, row.names = TRUE, col.names = FALSE)

# Filtering the list with the sample list from the sample table
tximport_file_list <- fullpathfiles[ samples_files %in% sample_table$sample ]
write.table((tximport_file_list), paste(FolderOutput ,"/Log/Tximport_Final list.txt", sep = ""), sep="\t", quote=FALSE, row.names = TRUE, col.names = FALSE)

## Create a dataframe with transcript ID & gene ID with ENSEMBL ref
install.packages(paste(ensembl_package), type = "source", repos = NULL)
package_name <- gsub(".*/(.*?)\\.tar\\.gz", "\\1", ensembl_package)
suppressPackageStartupMessages(library(paste(package_name),character.only=TRUE))
mv(from = package_name, to = "edb") 
# Reformat
k <- keys(edb, keytype = "TXNAME")
tx2gene <- select(edb, k, "GENENAME", "TXNAME")
tx2gene<-tx2gene[,1:2]
write.table(tx2gene, paste(FolderOutput ,"/Databases/Tx2gene_",package_name,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

## Tximport
# Use option countsFromAbundance="lengthScaledTPM" to obtain "raw" counts values (non-normalized) from the TPM (recommended for DE analysis)
# Option countsFromAbundance="scaledTPM" to take the TPM scaled up to library size as "raw" counts
# Default option : takes TPM as scaled values (Abundance) and NumReads as "raw" counts (Counts)
# type = c("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie")
# countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM")
tx2gene <- read.table(paste(FolderOutput ,"/Databases/Tx2gene_",package_name,".txt", sep = ""), header=T, stringsAsFactors = FALSE, sep="\t")
txi <- tximport(tximport_file_list, type = paste(CountType), tx2gene=tx2gene, ignoreTxVersion = TRUE, countsFromAbundance="lengthScaledTPM")

## Export counts data
data <- txi$counts %>% round(digits=0)
data <- cbind(rownames(data),data)
colnames(data)[1] <- "Gene"
write.table(data, paste(FolderOutput ,"/Counts/Tximport_Counts_per_gene_",level_to_compare,"_vs_",base_level,".txt", sep = ""), quote=FALSE, row.names = FALSE, col.names = TRUE)

## Sampletable/metadata check
# Import from the working directory a file matching the samples to the corresponding groups that are investigated
meta <- data.frame(sample_table, row.names = colnames(txi$counts))

# Check that the row names of the metadata are the same as the column names of the counts data 
# and that the column names of the counts data are the same as the sample names of metadata
# if not there's incoherence between sample table & files, files are treated in alphabetic order, so there's possibly some mismatches

print(colnames(txi$counts))
print(rownames(meta))

sampletest <- all(sort(colnames(txi$counts)) == sort(rownames(meta)))
if (!sampletest) {
  print("Check your sample's files count (missing file or wrong name compared to the sample table)")
  stop()
}

sampletabletest <- all(sort(colnames(txi$counts)) == sort(meta$sample))
if (!sampletabletest) {
  print("Check the sample's name in the sample table (missing file or wrong name compared to the sample's file count)")
  stop()
}

# Library
suppressPackageStartupMessages({
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(heatmaply)
library(EnhancedVolcano)
library(genefilter)
library(ReportingTools)
library(regionReport)
library(pcaExplorer)
library(apeglm)
})

### Step 2 ### Create the dds object
# design is the name(s) of the column(s) of the sample_table with the info on the groups/conditions you want to compare
# Need to figure out how to allow complex design like ~ cell+condition+cell:condition

if (is.null(covariable)) {
  print("No covarialbe")
	DESeq2Model <-  paste("~ ",Condition2Compare)
}else{
	DESeq2Model <-  paste("~ ",Condition2Compare," + ",covariable)
}

if (is.null(interaction)) {
	print("No effect")
}else{
	DESeq2Model <-  paste("~ ",Condition2Compare," + ",covariable," + ",interaction)
}

print("Modele is")
print(DESeq2Model)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = as.formula(DESeq2Model))

# Specify the base level
dds$condition <- relevel( dds$condition, paste(base_level))

## Pre-filtering the dataset (dds object) before DE analysis
## Delete row if 0 counts for all samples
keep <- rowSums(counts(dds)) >= 1 # Keep if sum of row values >= 1 
dds <- dds[keep,]

## Other filters (select only genes with significant amount of counts in a significant number of samples)
## Ex : Delete row if less than 3 samples have 3 counts or less for this gene
# Keep if >= 3 counts in >= 3 samples ; can be a variable
#keep <- rowSums(counts(dds) >= 3) >= 3
#dds <- dds[keep,]

### Step 3 ### Exploratory data analysis (PCA & hierarchical clustering) - QC control
# Identifying outliers and sources of variation in the data
# Hierarchical clustering places similar samples together, represented by the 
# tree structure. High correlations across the board (> 0.999) suggest no outlying sample(s).
# The samples are supposed to be clustering together by sample group

# Transform counts for data visualization (rlog transformation)
rld <- rlog(dds, blind=TRUE)

# PCA uses the top 500 most variable genes to determine the similarity of the samples
pdf(file = paste(FolderOutput ,"PCA_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 9, height = 7)

# Generate PCA plot without text labels first
pca_plot <- pcaplot(rld, intgroup = Condition2Compare, text_labels = FALSE, ellipse = TRUE, ntop = 500, title = "PCA plot using the top 500 variable genes")
# Add custom labels using ggrepel to avoid overlaps
pca_plot <- pca_plot + 
  geom_point(size = 3) + 
  geom_text_repel(aes(label = rownames(colData(rld))), size = 3, max.overlaps = Inf)
# Ensure the legend is included
pca_plot <- pca_plot + theme(legend.position = "right")
print(pca_plot)
dev.off()

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap for unsupervised hierarchical clustering (sample to sample distance)
pdf(file = paste(FolderOutput ,"HeatMap_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
pheatmap(rld_cor, annotation = meta)
dev.off()

###########################################################
### Step 4 ### Differential expression analysis with DESeq2
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis
# For example, include other sources of variation (if the info is in the sample_table) using
# "dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ covariable + condition)"
# Run DESeq2 differential expression analysis 
dds <- DESeq(dds)
resultsNames(dds)  # Check the default group comparison
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- cbind(rownames(normalized_counts),normalized_counts) # Add a first column with gene names
colnames(normalized_counts)[1]<-"gene"
write.table(normalized_counts, paste(FolderOutput ,"/Counts/DESeq_Normalized_Counts_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=F, col.names=TRUE, row.names=FALSE)

### Step 5 ### Check the fit of the dispersion estimates (Quality check)
# The dispersion estimates reflect the variance in gene expression for a given
# mean value. With few replicates per group, the estimates of variation for each gene
# are often unreliable. DESeq2 shares information across genes to generate more accurate
# estimates of variation using a "shrinkage" method (DESeq2 assumes that genes with
# similar expression levels should have similar dispersion).
# This is a good plot to examine to ensure the data is a good fit for the DESeq2 model. 
# Black dots = dispersion estimates for each gene, Blue dots = shrunken dispersion values.
# We expect the data to generally scatter around the curve, with the dispersion decreasing 
# with increasing mean expression levels.
# With more replicates per condition, less shrinkage is applied to the dispersion estimates.

# Plot dispersion estimates
pdf(file = paste(FolderOutput ,"/QC/DESeq2_Dispersion_Estimates_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
plotDispEsts(dds)
dev.off()

### Step 6 ### Create contrasts to perform Wald testing on the shrunken log2 foldchanges
# DESeq2 uses Wald test to identify genes that are differentially expressed 
# between two sample classes. Results for different comparisons can be extracted 
# depending on how many factors are used in the design formula and how many factor 
# levels are present.
# By default, DESeq compares the last condition against the first condition.
# Example of typical use : contrast <- c("condition", "level_to_compare", "base_level")
# Here, for the factor "condition", we want to compare "level_to_compare" versus "base_level"
# The fold change will be calculated with "level_to_compare" as numerator and "base_level"
# as denominator. Fold-change = level_to_compare/base_level.
# By default, DESeq2 allows for the shrinkage of the LFC estimates toward zero 
# when the information for a gene is low (low counts or high dispersion values). 
# LFC shrinkage uses information from all genes to generate more accurate estimates.
# The coef will be dependent on what your contrast was and should be identical to what 
# is stored in resultsNames()

# Use "resultsNames(dds)" to see the values that you can provide to the "coef" argument
# comparison <- c(Condition2Compare, level_to_compare, base_level) 

comparison_raw <- resultsNames(dds)[-1]

write.table(comparison_raw, paste(FolderOutput ,"/Log/Possible_Models.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
# Specify contrast for comparison of interest
# Output results of Wald test for contrast of interest
# ex : condition_S_LPS_vs_NS cell_LLC_vs_LB conditionS_LPS.cellLLC or comparison <- list(c("condition_S_LPS_vs_NS","cell_LLC_vs_LB"))
if (!is.null(interaction)) {
  selected_pattern <- paste0(Condition2Compare,level_to_compare)
  # Find indices where pattern matches
  indices <- grep(selected_pattern, comparison_raw, ignore.case = TRUE, value = FALSE)
  # Subset based on indices
  selected_values <- comparison_raw[indices]
  # Create a list
  comparison <- list(selected_values = selected_values)
  print('Values in model with interaction')
  print(selected_pattern)
  print(indices)
  print(selected_values)
} else {
  comparison <- list(c(comparison_raw))
   print('Values in model without interaction')
  print(comparison)
}

write.table(comparison, paste(FolderOutput ,"/Log/Final_Model.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

# Compute res with contrast selected
res <- results(dds, contrast = comparison , alpha = 0.05)
# Save unshrunken results
res_tableOE_unshrunken <- res 
# Shrink the log2 fold changes (LFC) to be more accurate
# "apeglm" is the adaptive Student's t prior shrinkage estimator from the 'apeglm' package
# apeglm requires use of coef
# For type="apeglm": Specifying apeglm passes along DESeq2 MLE log2 fold changes and standard errors to the apeglm function in the apeglm package, and re-estimates posterior LFCs for the coefficient specified by coef.
# "ashr" is the adaptive shrinkage estimator from the 'ashr' package, using a fitted mixture of normals prior - see the Stephens (2016) reference below for citation
# For type="ashr": Specifying ashr passes along DESeq2 MLE log2 fold changes and standard errors to the ash function in the ashr package, with arguments mixcompdist="normal" and method="shrink".
# "normal" is the 2014 DESeq2 shrinkage estimator using a Normal prior, this algorithm is out of date compared to the 2 below
res <- lfcShrink(dds, coef = 2, type = paste(shrink_method)) # recommandation for max coeff is 4
# Save shrunken results, res is now shrunken
res_tableOE_shrunken <- res 
# Compare unshrunken vs shrunken results with MA plot
# MA plot shows the mean of the normalized counts versus the log2foldchanges for all genes tested. The genes that are significantly DE are colored in blue
pdf(file = paste(FolderOutput ,"/QC/MAplot_unshrunken_foldchanges_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
plotMA(res_tableOE_unshrunken, ylim=c(-2,2)) # MA plot using unshrunken fold changes ; scale is +/- 4
abline(h=c(-1,1), col="dodgerblue", lwd=2) # add a blue line for +1/-1 fold change
dev.off()
pdf(file = paste(FolderOutput ,"/QC/MAplot_shrunken_foldchanges_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
plotMA(res_tableOE_shrunken, ylim=c(-2,2)) # MA plot using shrunken fold changes ; scale is +/- 4
abline(h=c(-1,1), col="dodgerblue", lwd=2) # add a blue line for +1/-1 fold change 
dev.off()

### Step 7 ### Output significant results
# This function reports the number of genes up- and down-regulated at the selected threshold (padj/alpha), the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Set thresholds (padj cutoff = FDR)
padj.cutoff.01 <- 0.01
padj.cutoff.05 <- 0.05
# Subset the significant results
sig_res_01 <- dplyr::filter(res_tbl, padj < padj.cutoff.01)
sig_res_05 <- dplyr::filter(res_tbl, padj < padj.cutoff.05)
sig_res_01_LFC1 <- dplyr::filter(sig_res_01, log2FoldChange > 1)
sig_res_01_LFC0 <- dplyr::filter(sig_res_01, log2FoldChange > 0)
sig_res_01_LFCminus1 <- dplyr::filter(sig_res_01, log2FoldChange < -1)
sig_res_01_LFCminus0 <- dplyr::filter(sig_res_01, log2FoldChange < 0)
sig_res_05_LFC1 <- dplyr::filter(sig_res_05, log2FoldChange > 1)
sig_res_05_LFC0 <- dplyr::filter(sig_res_05, log2FoldChange > 0)
sig_res_05_LFCminus1 <- dplyr::filter(sig_res_05, log2FoldChange < -1)
sig_res_05_LFCminus0 <- dplyr::filter(sig_res_05, log2FoldChange < 0)

# Save results
write.table(res_tbl, paste(FolderOutput ,"/FoldChanges/DEgenes_raw_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
# alpha 0,01
write.table(sig_res_01, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres01_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_01_LFC1, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres01_LFC1_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_01_LFCminus1, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres01_LFC-1_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
# alpha 0,05
write.table(sig_res_05, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres05_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_05_LFC1, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres05_LFC1_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_05_LFCminus1, paste(FolderOutput ,"/FoldChanges/DEgenes_sigres05_LFC-1_",level_to_compare,"_vs_",base_level,".txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

### Step 8 ### Visualize results: volcano plots, heatmaps, normalized counts plots of top genes...
## 1 # Create reports from DESeq2 results with RegionReport = HTML and PDF Report with visualization plots
regionreport <- DESeq2Report(dds, project = "DE report with RegionReport", res=res, intgroup = Condition2Compare, outdir = FolderOutput, output = paste("RegionReport_",level_to_compare,"_vs_",base_level, sep = ""))

## if res=NULL, then results will be used on dds with default parameters (last condition vs first condition). 

## 2 # Plotting significant DE genes
## 2.1 # Plot expression for a single gene of interest (e.g. CD180)
# Use DESeq2 plotCounts function to plot expression for a single gene and save as PDF
# The specified gene needs to match with the original input to DESeq2 (GeneID, EnsemblID...)
#pdf(file = paste(FolderOutput ,"/Plots/", geneplot, "_Expression.pdf", sep = ""))
#plotCounts(dds, gene=geneplot, intgroup=Condition2Compare)
#dev.off()

# Use ggplot2 to plot the normalized counts using the samplenames (rownames(d) as labels)
# Save plotcounts to a data frame object
#d <- plotCounts(dds, gene=geneplot, intgroup=Condition2Compare, returnData=TRUE)

# Create graph
#pdf(file = paste(FolderOutput ,"/Plots/", geneplot, "_Expression_Labeled.pdf", sep = ""))
#ggplot(d, aes(x = .data[[Condition2Compare]], y = count, color = .data[[Condition2Compare]])) + geom_point(position=position_jitter(w = 0.1,h = 0)) + geom_text_repel(aes(label = rownames(d))) + ggtitle(geneplot) + theme(plot.title = element_text(hjust = 0.5))
#dev.off()

## 2.2 # Plot expression for multiple genes of interest (e.g. top 20 DE genes)
## Normalized counts for all significant DE genes
normalized_counts <- as_tibble(normalized_counts)
AllDEgenes<-sig_res_01$gene # p < 0.01
AllDEgenes_norm <- normalized_counts %>%
  dplyr::filter(gene %in% AllDEgenes)
write.table(AllDEgenes_norm, paste(FolderOutput ,"/Counts/AllDEgenes_Normcounts_",level_to_compare,"_vs_",base_level,".txt", sep = "") , sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

## Order results by padj values 
topx_sigOE_genes <- sig_res_01 %>% 
  arrange(padj) %>% 	# Arrange rows by padj values
  pull(gene) %>% 		# Extract character vector of ordered genes
  head(n=topgene)		# Extract the first x top genes

## Normalized counts for top x significant genes
topx_sigOE_norm <- normalized_counts %>%
  dplyr::filter(gene %in% topx_sigOE_genes)

# Gathering the columns to have normalized counts to a single column , added after last column [2:numberofsamples+1]
x <- nrow(sample_table) +1

gathered_topx_sigOE <- topx_sigOE_norm %>%
  gather(colnames(topx_sigOE_norm)[2:x], key = "sample", value = "normalized_counts")

topx_sigOE_final <- gathered_topx_sigOE %>% dplyr::select(gene, sample) # keep gene and sample columns

topx_sigOE_final <- topx_sigOE_final %>% cbind(as.numeric(gathered_topx_sigOE$normalized_counts)) # join normalized_counts column as numeric format
colnames(topx_sigOE_final)[3]<-"normalized_counts"

## Merge metadata information with the normalized counts data (merge or inner_join)
# inner_join() will merge 2 data frames with respect to the "sample" column, 
# i.e. a column with the same column name in both data frames.
topx_sigOE_final <- inner_join(meta, topx_sigOE_final)

## Plot using ggplot2 top x genes
pdf(file= paste(FolderOutput ,"Top",topgene,"_DE_Genes_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
ggplot(topx_sigOE_final) +
  geom_point(aes(x = gene, y = normalized_counts, color = .data[[Condition2Compare]])) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Log Normalized Counts") +
  ggtitle(paste("Top Significant ",topgene," DE Genes ",level_to_compare,"_vs_",base_level)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## 3 # Gene clustering (Heatmap) using Heatmaply
## Use transformed counts by rlog transformation ("rld") and select top x DE genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), topheatgene)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

## Create interactive heatmap using heatmaply and save as html file 
#heatmaply(mat, file = paste(FolderOutput ,"heatmaply_plot_defaultcolors.html", sep = ""))
## Personalize colors (e.g. blue for low and red for high expression) save and open
heatmaply(mat, file =  paste(FolderOutput ,"HeatMap_Top",topheatgene,"_genes_",level_to_compare,"_vs_",base_level,".html", sep = ""), scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = c(-4, 4), ))

## 4 # Volcano plot (EnhancedVolcano)
pdf(file= paste(FolderOutput ,"/Plots/DE_Volcanoplot_All_Genes_",level_to_compare,"_vs_",base_level,".pdf", sep = ""))
EnhancedVolcano(sig_res_01,
                lab = sig_res_01$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DE_Volcanoplot',
                pCutoff = 0.01,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 4.0,
                legendLabSize = 9,
                col=c('grey','grey','grey', 'red2'),
                colAlpha = 1)
dev.off()

# TODO list of gene should be an arguement
#pdf(file= paste(FolderOutput ,"/Plots/DE_Volcanoplot_Selection_of_Genes.pdf", sep = ""))
#EnhancedVolcano(sig_res_01,
#                lab = sig_res_01$gene,
#                x = 'log2FoldChange',
#                y = 'pvalue',
#                title = 'DE_Volcanoplot_Selection_of_Genes',
#                selectLab = c('TLR1','TLR2','TLR3','TLR4', 'TLR5', 'TLR6', 'TLR7', 'TLR8', 'TLR9', 'TLR10', 'CD180', 'LY86', 'LY96', 'NFKB1', 'NFKBIA', 'NFKBID'),
#                pCutoff = 0.01,
#                FCcutoff = 1.0,
#                cutoffLineType = 'twodash',
#                cutoffLineWidth = 0.8,
#                pointSize = 1.0,
#                labSize = 3.0,
#                labCol = 'black',
#                labFace = 'bold',
#                col=c('grey', 'grey', 'grey', 'red2'),
#                colAlpha = 0.75,
#                legendLabels=c('NS','Log FC','p-value',
#                               'p-value & Log FC'),
#                legendPosition = 'right',
#                legendLabSize = 10,
#                legendIconSize = 5.0,
#                drawConnectors = TRUE,
#                widthConnectors = 1.0,
#                colConnectors = 'black')
#dev.off()


### Step 9 ### Functional analysis: GO over-representation, GSEA with clusterProfiler
## 1 # Gene Ontology (GO) over-representation analysis with clusterProfiler

# Load libraries
suppressPackageStartupMessages({
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
})
## Ensembl and EntrezID are needed for this functional analysis
# Create a gene-level dataframe 
annotations_ahb <- genes(edb, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype)
# keep the first identifier for these multiple mapping cases
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)
# Return only the non-duplicated genes using indices
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]
# Save annotations_ahb
write.csv(as.data.frame(annotations_ahb), file= paste(FolderOutput ,"/Databases/Annotations_ahb_v107.csv", sep = ""))
## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_tbl, annotations_ahb, by=c("gene"="gene_name"))
## Create background dataset for hypergeometric testing = all genes tested
allOE_genes <- as.character(res_ids$gene)
# sig_res_01$gene is shrunken
sigOE_genes <- as.character(sig_res_01$gene)
#sigOE_genes <- as.character(sig_res_05$gene)

## Run GO enrichment analysis (significant OE genes within universal OE genes)
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = F)

## Output results from GO analysis to a table
cluster_summary <- as.data.frame(ego)
write.csv(as.data.frame(cluster_summary), file= paste(FolderOutput ,"/GSEA/GO_ClusterProfiler_",level_to_compare,"_vs_",base_level,".csv", sep = ""))


## Run GO enrichment analysis for 2 fold changes only (not very accurate statistically speaking)
#sigOE_genes_2fold <- as.character(sig_res_01_LFC1$gene)
#sigOE_genes_2fold <- as.character(sig_res_01_LFCminus1$gene)

# Run GO enrichment analysis (significant OE genes 2 fold within universal OE genes)
#ego2fold <- enrichGO(gene = sigOE_genes_2fold, 
#                universe = allOE_genes,
#                keyType = "SYMBOL",
#                OrgDb = org.Hs.eg.db, 
#                ont = "BP", 
#                pAdjustMethod = "BH", 
#                qvalueCutoff = 0.05, 
#                readable = F)

## Output results from GO analysis to a table
#cluster_summary_2fold <- as.data.frame(ego2fold)
#write.csv(as.data.frame(cluster_summary_2fold), file= paste(FolderOutput ,"/GSEA/GO_ClusterProfiler_2Fold_",level_to_compare,"_vs_",base_level,".csv", sep = ""))


## Visualize clusterProfiler results
pdf(file = paste(FolderOutput ,"/GSEA/GO_ClusterProfiler_Dotplot_Top",topgene,level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
dotplot(ego, showCategory=topgene, font.size = 8)
dev.off()

# Enrichmap #
# Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)
# Enrichmap clusters the x most significant (by padj) GO terms to visualize relationships between terms
pdf(file = paste(FolderOutput ,"/GSEA/GO_Cluster_Profiler_Enrichmap_Top_",topgene,"_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
emapplot(ego, showCategory = topgene, cex.params = list(category_label = 0.5, line = 0.25))
dev.off()

# Category netplot #
# To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
sigOE <- dplyr::filter(res_ids, padj < 0.01)
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
# Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
pdf(file = paste(FolderOutput ,"/GSEA/GO_Clusterprofiler_Category_Netplot_All_Genes_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
cnetplot(ego, categorySize="pvalue", showCategory = 10, color.params = list(foldChange = OE_foldchanges), cex.params = list(category_label = 0.5, gene_label = 0.5), vertex.label.font=0.05, max.overlaps = 25)
dev.off()

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)

pdf(file = paste(FolderOutput ,"/GSEA/GO_Clusterprofiler_Category_Netplot_Filter_by_2_FoldChanges_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
cnetplot(ego, categorySize="pvalue", showCategory = 10, color.params = list(foldChange = OE_foldchanges), cex.params = list(category_label = 0.5, gene_label = 0.5), vertex.label.font=0.05, max.overlaps = 25)
dev.off()


## GSEA analysis with clusterProfiler
# Functional class scoring (FCS) tools, such as GSEA, most often use the gene-level statistics or log2 fold changes for all genes from the differential expression results
# then look to see whether gene sets for particular biological pathways are enriched among the large positive or negative fold changes.
# https://www.gsea-msigdb.org/gsea/index.jsp
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA
# https://github.com/hbctraining/DGE_workshop_salmon_online/blob/master/lessons/10_FA_over-representation_analysis.md

# Preparation for GSEA
# To use the KEGG gene sets, we need to acquire Entrez IDs and remove NA values and duplicates
# Remove any NA values
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
# Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]

# Write a table with all datas
write.csv(as.data.frame(res_entrez), file= paste(FolderOutput ,"/GSEA/Genes_DE_for_GSEA_",level_to_compare,"_vs_",base_level,".csv", sep = ""))
# Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
# Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid
# Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

# First, we will set the seed so that we all obtain the same result:
set.seed(123456)
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # hsa for Homo Sapiens
                    minGSSize = 5, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer genes default 5
                    pvalueCutoff = 0.05, # padj cutoff value (generally 0.05)
                    verbose = FALSE)

write.csv(as.data.frame(gseaKEGG), file= paste(FolderOutput ,"/KEGG/KEGG_Results_",level_to_compare,"_vs_",base_level,".csv", sep = ""))

# Print pathways
gseaKEGG_results <- gseaKEGG@result

## Output images for all significant KEGG pathways
# Define a safe version of get_kegg_plots using purrr::safely
safe_get_kegg_plots <- purrr::safely(function(x) {
   pathview(gene.data = foldchanges, 
            pathway.id = gseaKEGG_results$ID[x], 
            species = "hsa",
            limit = list(gene = 2, cpd = 1),
            kegg.dir = paste(FolderOutput, "/KEGG/", sep = ""))
})

# Filter out specific KEGG pathways
#gseaKEGG_results <- gseaKEGG_results %>%
#  dplyr::filter(!rownames(gseaKEGG_results) %in% c(
#    'hsa04215', 'hsa05206', 'hsa01240', 'hsa01200', 
#    'hsa01230', 'hsa01212', 'hsa01210', 'hsa01232', 
#    'hsa01250', 'hsa01040'))

# Set working directory
setwd(paste(FolderOutput, "/KEGG/", sep = ""))

# Map over the pathway IDs safely
results <- purrr::map(1:length(gseaKEGG_results$ID), safe_get_kegg_plots)

# Optionally, you can inspect errors in results by checking:
errors <- purrr::map(results, "error")  # To see which pathways had errors

# clusterProfiler::GSEA

# The SPIA (Signaling Pathway Impact Analysis) tool can be used to integrate the lists of differentially expressed genes, their fold changes, and pathway
# topology to identify affected pathways.
suppressPackageStartupMessages(library(SPIA))

# Significant genes is a vector of fold changes where the names are ENTREZ gene IDs. The background set is a vector of all the genes represented on the platform.
background_entrez <- res_entrez$entrezid
sig_res_entrez <- res_entrez[which(res_entrez$padj < 0.05), ] # using p value of 0.05
write.csv(as.data.frame(sig_res_entrez), file= paste(FolderOutput ,"/SPIA/SPIA_Significant_Gene_",level_to_compare,"_vs_",base_level,".csv", sep = ""))
sig_entrez <- sig_res_entrez$log2FoldChange
names(sig_entrez) <- sig_res_entrez$entrezid
spia_result <- spia(de=sig_entrez, all=background_entrez, organism="hsa")

spia_result$pG=combfunc(spia_result$pNDE,spia_result$pPERT,combine="norminv")
spia_result$pGFdr=p.adjust(spia_result$pG,"fdr")
spia_result$pGFWER=p.adjust(spia_result$pG,"bonferroni")
write.csv(as.data.frame(spia_result), file= paste(FolderOutput ,"/SPIA/Spia_Result_",level_to_compare,"_vs_",base_level,".csv", sep = ""))

## Table Legend ##
# pSize: number of genes on the pathway
# NDE: number of DE genes per pathway
# tA: observed total perturbation accumulation in the pathway
# pNDE: probability to observe at least NDE genes on the pathway using a hypergeometric model (similar to ORA)
# pPERT: probability to observe a total accumulation more extreme than tA only by chance
# pG: the p-value obtained by combining pNDE and pPERT
# pGFdr and pGFWER are the False Discovery Rate and Bonferroni adjusted global p-values, respectively
# Status: gives the direction in which the pathway is perturbed (activated or inhibited)
# KEGGLINK : web link to the KEGG website that displays the pathway image with the differentially expressed genes highlighted in red

# Plot SPIA results
# In this plot, each pathway is a point and the coordinates are the log of pNDE (using a hypergeometric model) and the p-value from perturbations, pPERT. The
# oblique lines in the plot show the significance regions based on the combined evidence
source('SPIA_plot_fork.R')
pdf(file = paste(FolderOutput ,"/SPIA/SPIA_Plot_Results_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
plotP_fork(spia_result) # we can specify other threshold value
dev.off()

### Output the versions of all tools used in the DE analysis
sessionInfo()