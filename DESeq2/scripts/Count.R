# Workflow for DE analysis using DESeq2 #
# Import quantification via any tools (salmon/alevin is default tools)

## Original author : L. RIGOLOT
## Automomatisation and maintenance : T LAVAUX

# 15/09/2022 : add argparse options and refactor some code - TL
# 01/12/2022 : add KEGG patwhay graphs - TL
# 08/12/2022 : add Ensembl update package (v107), improve files naming & graphs (GO) - TL
# 19/01/2022 : add covar option & update SPIA PlotP to be more robust ; refactor some code - TL
# 02/02/2023 : add effect option for the statistic model, sample table filter samples list from files ; some variabilisation & options added - TL

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
suppressWarnings(library(ggplot2))
suppressWarnings(library(optparse))
suppressWarnings(library(stringr))
suppressWarnings(library(tximport))
suppressWarnings(library(tidyverse))
suppressWarnings(library(gdata))

######Parsing input options and setting defaults########
option_list<-list(
	make_option('--input',default='data', help='Folder were the data are stored',dest='datafolder'),
	make_option('--output',default='results', help='Folder were to save the results',dest='resultfolder'),
	make_option('--count',default='salmon', help='Type of RNA count',dest='count'),
	make_option('--sampletable',default='./sample_table.csv', help='Path of the sample table',dest='samtable'),
	make_option('--compare',default='condition', help='Comparison to be analysed',dest='cond'),
	make_option('--level2compare',default='level_to_compare', help='Level to compare with',dest='level2com'),
	make_option('--baselevel',default='base_level', help='Base level',dest='baselvl'),
	make_option('--annotdb',default='./database/annotations_ahb.csv', help='Annotation database',dest='annotdb'),
	make_option('--gene',default='CD180', help='Gene name for plot construction',dest='gene'),
	make_option('--topplot',default=30, help='Number of top gene taken into account for top plots',dest='topgene'),
	make_option('--topheatmap',default=50, help='Number of top genes taken into account for interactive heatmap',dest='topheat'),
	make_option('--covar',default='', help='Covariable to add to the statistic model',dest='covar'),
	make_option('--eff',default='', help='Effect to add to the statistic model, separated by a :',dest='eff'),
	make_option('--shrink',default='ashr', help='Shrinking algorithm for DE analysis, choose between normal, ashr (default) or apeglm',dest='shrink'),
	make_option('--ensembl',default='./database/EnsDb.Hsapiens.v107.tar.gz', help='Ensembl package file in .tar.gz format',dest='ensembl')
)
opt<-parse_args(OptionParser(option_list=option_list))

FolderOutput=opt$resultfolder
DataInput=opt$datafolder
CountType=opt$count
SampleTable=opt$samtable
Condition2Compare=opt$cond
level_to_compare=opt$level2com
base_level=opt$baselvl
annotdbfilepath=opt$annotdb
geneplot=opt$gene
topgene=opt$topgene
topheatgene=opt$topheat
covariable=opt$covar
effect=opt$eff
shrink_method=opt$shrink
ensembl_package=opt$ensembl

### STEP 1 ### Setup to import datas into proper files
# Create output dir
dir.create(FolderOutput)
dir.create(paste(FolderOutput ,"/Counts", sep = ""))
dir.create(paste(FolderOutput ,"/FoldChanges", sep = ""))
dir.create(paste(FolderOutput ,"/Plots", sep = ""))
dir.create(paste(FolderOutput ,"/KEGG", sep = ""))
dir.create(paste(FolderOutput ,"/GSEA", sep = ""))
dir.create(paste(FolderOutput ,"/SPIA", sep = ""))
dir.create(paste(FolderOutput ,"/QC", sep = ""))
dir.create(paste(FolderOutput ,"/Databases", sep = ""))
dir.create(paste(FolderOutput ,"/Log", sep = ""))

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
suppressMessages(library(paste(package_name),character.only=TRUE))
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

#sampletest <- all(colnames(txi$counts) == rownames(meta))
#if (!sampletest) {
#	print("Check your samples files (missing file or wrong name compared to the sample table)")
#	stop()
#}

#sampletabletest <- all(colnames(txi$counts) == meta$sample)
#if (!sampletabletest) {
#	print("Check the sample list in the sample table (missing file or wrong name compared to the sample table)")
#	stop()
#}

# Library
suppressWarnings(library(DESeq2))
suppressWarnings(library(ggrepel))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(pheatmap))
suppressWarnings(library(DEGreport))
suppressWarnings(library(heatmaply))
suppressWarnings(library(EnhancedVolcano))
suppressWarnings(library(genefilter))
suppressWarnings(library(ReportingTools))
suppressWarnings(library(regionReport))
suppressWarnings(library(pcaExplorer))
suppressWarnings(library(apeglm))

### Step 2 ### Create the dds object
# design is the name(s) of the column(s) of the sample_table with the info on the groups/conditions you want to compare
# Need to figure out how to allow complex design like ~ cell+condition+cell:condition

if (covariable == ''){
	DESeq2Model <-  paste("~ ",Condition2Compare)
}else{
	DESeq2Model <-  paste("~ ",Condition2Compare," + ",covariable)
}

if (effect == ''){
	print("No effect")
}else{
	DESeq2Model <-  paste("~ ",Condition2Compare," + ",covariable," + ",effect)
}

print("Modele is")
print(DESeq2Model)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = as.formula(DESeq2Model))

# Specify the base level
# TODO condition should be a variable
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
pcaplot(rld, intgroup= Condition2Compare, text_labels = TRUE, ellipse = TRUE, ntop = 500 , title = "PCA plot using the top 500 variable genes")
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
comparison <- list(c(comparison_raw))
write.table(comparison, paste(FolderOutput ,"/Log/Possible_Models.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
# Specify contrast for comparison of interest
# Output results of Wald test for contrast of interest
# ex : condition_S_LPS_vs_NS cell_LLC_vs_LB conditionS_LPS.cellLLC or comparison <- list(c("condition_S_LPS_vs_NS","cell_LLC_vs_LB"))
# TODO rewrite the comparison to include constrast depending on the variables we specify

write.table(comparison, paste(FolderOutput ,"/Log/Final_Model.txt", sep = ""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
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
# Save shrunken results
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
# This function reports the number of genes up- and down-regulated at the selected
# threshold (padj/alpha), the number of genes that were tested (genes with non-zero 
# total read count), and the number of genes not included in multiple test correction 
# due to a low mean count

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
## 1 # Create reports from DESeq2 results with special R packages
## RegionReport = HTML and PDF Report with visualization plots
regionreport <- DESeq2Report(dds, project = "DE report with RegionReport", res=NULL, intgroup = Condition2Compare, outdir = FolderOutput, output = paste("RegionReport_",level_to_compare,"_vs_",base_level, sep = ""))

## if res=NULL, then results will be used on dds with default parameters (last 
## condition vs first condition). Use the "res" object created at step 6 to generate 
## a RegionReport for the contrast of interest = replace res by is("res object", "class")

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

## Order results by padj values #res_tbl
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
suppressWarnings(library(DOSE))
suppressWarnings(library(pathview))
suppressWarnings(library(clusterProfiler))
suppressWarnings(library(org.Hs.eg.db))
suppressWarnings(library(ggnewscale))

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
sigOE_genes_2fold <- as.character(sig_res_01_LFC1$gene)
sigOE_genes_2fold <- as.character(sig_res_01_LFCminus1$gene)

# Run GO enrichment analysis (significant OE genes 2 fold within universal OE genes)
ego2fold <- enrichGO(gene = sigOE_genes_2fold, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = F)

## Output results from GO analysis to a table
cluster_summary_2fold <- as.data.frame(ego2fold)
write.csv(as.data.frame(cluster_summary_2fold), file= paste(FolderOutput ,"/GSEA/GO_ClusterProfiler_2Fold_",level_to_compare,"_vs_",base_level,".csv", sep = ""))


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
get_kegg_plots <- function(x) {
   pathview(gene.data = foldchanges, 
            pathway.id = gseaKEGG_results$ID[x], 
            species = "hsa",
            limit = list(gene = 2, cpd = 1),
			kegg.dir = paste(FolderOutput ,"/KEGG/", sep = ""))
}

# need to remove some pathview because it crash the analysis
gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa04215')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa05206')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01240')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01200')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01230')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01212')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01210')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01232')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01250')

gseaKEGG_results <- gseaKEGG_results %>%
  dplyr::filter(rownames(gseaKEGG_results) != 'hsa01040')


setwd(paste(FolderOutput ,"/KEGG/", sep = ""))
try(purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots))

# clusterProfiler::GSEA

# The SPIA (Signaling Pathway Impact Analysis) tool can be used to integrate the lists of differentially expressed genes, their fold changes, and pathway
# topology to identify affected pathways.
suppressWarnings(library(SPIA))

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

# Fork to SPIA Plot function
plotP_fork<-function(x,threshold=0.05){

if(class(x)!="data.frame" | dim(x)[1]<1 | !all(c("ID","pNDE","pPERT","pG","pGFdr","pGFWER")%in%names(x)))
{
 stop("plotP can be applied only to a dataframe produced by spia function!!!") 
}

if(threshold<x[1,"pGFdr"]){
print(paste("The threshold value was corrected to be equal to ",x[1,"pGFdr"]))
threshold <- x[1,"pGFdr"]
}

pb<-x[,"pPERT"]
ph<-x[,"pNDE"]

#determine what combine method was used to convert ph and pb into pG
combinemethod=ifelse(sum(combfunc(pb,ph,"fisher")==x$pG)>sum(combfunc(pb,ph,"norminv")==x$pG),"fisher","norminv")

okx<-(ph<1e-6)
oky<-(pb<1e-6)
ph[ph<1e-6]<-1e-6
pb[pb<1e-6]<-1e-6

plot(-log(ph),-log(pb),xlim=c(0,max(c(-log(ph),-log(pb))+1,na.rm=TRUE)),
 ylim=c(0,max(c(-log(ph),-log(pb)+1),na.rm=TRUE)),pch=19,main="SPIA two-way evidence plot",cex=1.5,
 xlab="-log(P NDE)",ylab="-log(P PERT)")
tr<-threshold/dim(na.omit(x))[1]

#abline(v=-log(tr),lwd=1,col="red",lty=2)
#abline(h=-log(tr),lwd=1,col="red",lty=2)

if(combinemethod=="fisher"){
points(c(0,-log(getP2(tr,"fisher")^2)),c(-log(getP2(tr,"fisher")^2),0),col="red",lwd=2,cex=0.7,type="l")
}else{
somep1=exp(seq(from=min(log(ph)),to=max(log(ph)),length=200))
somep2=pnorm(qnorm(tr)*sqrt(2)-qnorm(somep1))
points(-log(somep1),-log(somep2),col="red",lwd=2,cex=0.7,type="l") 
}

trold=tr
tr<-max(x[,"pG"][x[,"pGFdr"]<=threshold])
if(tr<=trold){tr=trold*1.03}

if(combinemethod=="fisher"){
points(c(0,-log(getP2(tr,"fisher")^2)),c(-log(getP2(tr,"fisher")^2),0),col="blue",lwd=2,cex=0.7,type="l")
 }else{
somep1=exp(seq(from=min(log(ph)),to=max(log(ph)),length=200))
somep2=pnorm(qnorm(tr)*sqrt(2)-qnorm(somep1))
points(-log(somep1),-log(somep2),col="blue",lwd=2,cex=0.7,type="l") 
}

#abline(v=-log(tr),lwd=1,col="blue",lty=2)
#abline(h=-log(tr),lwd=1,col="blue",lty=2)

oks<-x[,"pGFWER"]<=threshold
oks2<-x[,"pGFdr"]<=threshold

points(-log(ph)[oks2],-log(pb)[oks2],pch=19,col="blue",cex=1.5)
points(-log(ph)[oks],-log(pb)[oks],pch=19,col="red",cex=1.5)

if(sum(oks)>0){
 text(-log(ph)[oks]+0.70,-log(pb)[oks],labels=as.vector(x$ID)[oks],cex=0.65)
}
if(sum(oks2)>0){
 text(-log(ph)[oks2]+0.70,-log(pb)[oks2],labels=as.vector(x$ID)[oks2],cex=0.65)
}

testx <- sum(okx)
if (!(is.na(testx)) && (sum(okx)>0)) {
	points(-log(ph)[okx]-0.12,-log(pb)[okx],pch="|",col="black",cex=1.5)
}

testy <- sum(oky)
if (!(is.na(testy)) && (sum(oky)>0)) {
	points(-log(ph)[oky],-log(pb)[oky]-0.12,pch="_",col="black",cex=1.5)
}

}

# Plot SPIA results
# In this plot, each pathway is a point and the coordinates are the log of pNDE (using a hypergeometric model) and the p-value from perturbations, pPERT. The
# oblique lines in the plot show the significance regions based on the combined evidence
pdf(file = paste(FolderOutput ,"/SPIA/SPIA_Plot_Results_",level_to_compare,"_vs_",base_level,".pdf", sep = ""), width = 12, height = 14)
plotP_fork(spia_result) # we can specify other threshold value
dev.off()


### Output the versions of all tools used in the DE analysis
#sessionInfo()
