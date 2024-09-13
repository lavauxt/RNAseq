# TODO
QFCZ5_WT_filtered<- SCTransform(QFCZ5_WT_filtered, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = TRUE, return.only.var.genes = FALSE, seed.use = 1448145)

QFCZ5_WT_filtered <- RunPCA(QFCZ5_WT_filtered, features = VariableFeatures(object = QFCZ5_WT_filtered))

ElbowPlot(QFCZ5_WT_filtered, ndims = 50)

QFCZ5_WT_filtered <- FindNeighbors(QFCZ5_WT_filtered, dims = 1:15)

QFCZ5_WT_filtered <- FindClusters(QFCZ5_WT_filtered, resolution = 0.5)

QFCZ5_WT_filtered <- RunUMAP(QFCZ5_WT_filtered, dims = 1:15, seed.use = 1323)

DimPlot(QFCZ5_WT_filtered, reduction = "umap", label = TRUE, label.box = TRUE)+ NoLegend()

FeaturePlot(QFCZ5_WT_filtered, features = "Gapdh", order = TRUE)

FeaturePlot(QFCZ5_WT_filtered, features = "Il2ra", order = TRUE)

Idents(QFCZ5_WT_filtered)<- QFCZ5_WT_filtered$seurat_clusters

cluster_markers<-FindAllMarkers(object = QFCZ5_WT_filtered)
View(cluster_markers)

DimPlot(object = QFCZ5_WT_filtered, reduction = "umap", label = TRUE, label.box = TRUE, group.by = "seurat_clusters") + NoLegend()

FeaturePlot(object = QFCZ5_WT_filtered, reduction = "umap", features = "Hba-a1")

FeaturePlot(object = QFCZ5_WT_filtered, reduction = "umap", features = "Sry")

# Plot filtered
plot1_filtered <- FeatureScatter(QFCZ5_WT_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_filtered <- FeatureScatter(QFCZ5_WT_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_filtered + plot2_filtered

# Testing some graphs
metadata <- QFCZ5_WT@meta.data
View(metadata)

> metadata %>% 
+     ggplot(aes(x=orig.ident, fill=orig.ident)) + 
+     geom_bar() +
+     theme_classic() +
+     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
+     theme(plot.title = element_text(hjust=0.5, face="bold")) +
+     ggtitle("NCells")

metadata %>% 
+     ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
+     geom_density(alpha = 0.2) + 
+     scale_x_log10() + 
+     theme_classic() +
+     ylab("Cell density") +
+     geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~orig.ident)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)

# Normalize
QFCZ5_WT_filtered <- NormalizeData(QFCZ5_WT_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Cell cycle
QFCZ5_WT_phase <- CellCycleScoring(QFCZ5_WT_filtered, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(QFCZ5_WT_phase@meta.data)

# Identify the most variable genes
QFCZ5_WT_phase <- FindVariableFeatures(QFCZ5_WT_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
QFCZ5_WT_phase <- ScaleData(QFCZ5_WT_phase)

# Perform PCA
QFCZ5_WT_phase <- RunPCA(QFCZ5_WT_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(QFCZ5_WT_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")


QFCZ5_WT_filtered <- NormalizeData(QFCZ5_WT_filtered, verbose = TRUE)
QFCZ5_WT_filtered <- CellCycleScoring(QFCZ5_WT_filtered, g2m.features=g2m_genes, s.features=s_genes)
QFCZ5_WT_filtered <- SCTransform(QFCZ5_WT_filtered, vars.to.regress = c("percent.mt"))

# Select the most variable features to use for integration
QFCZ5_WT_integ_features <- SelectIntegrationFeatures(object.list = QFCZ5_WT_filtered, 
                                            nfeatures = 3000)

# Prepare the SCT list object for integration
QFCZ5_WT_integ_features <- PrepSCTIntegration(object.list = QFCZ5_WT_filtered, 
                                   anchor.features = QFCZ5_WT_integ_features)

# Find significant genes
QFCZ5_WT_filtered <- FindVariableFeatures(QFCZ5_WT_filtered, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(QFCZ5_WT_filtered), 10)

plot1 <- VariableFeaturePlot(QFCZ5_WT_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# PCA
all.genes <- rownames(QFCZ5_WT_filtered)
QFCZ5_WT_filtered <- ScaleData(QFCZ5_WT_filtered, features = all.genes)
QFCZ5_WT_filtered <- RunPCA(QFCZ5_WT_filtered, features = VariableFeatures(object = QFCZ5_WT_filtered))

# Exploring PCA
VizDimLoadings(QFCZ5_WT_filtered, dims = 1:2, reduction = "pca")
DimPlot(QFCZ5_WT_filtered, reduction = "pca") + NoLegend()
DimHeatmap(QFCZ5_WT_filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(QFCZ5_WT_filtered, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(QFCZ5_WT_filtered)

# Find clusters
QFCZ5_WT_filtered <- FindNeighbors(QFCZ5_WT_filtered, dims = 1:15)
QFCZ5_WT_filtered <- FindClusters(QFCZ5_WT_filtered, resolution = 0.5)

# UMAP
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
QFCZ5_WT_filtered <- RunUMAP(QFCZ5_WT_filtered, dims = 1:15)
DimPlot(QFCZ5_WT_filtered, reduction = "umap")


# UMAP 10
QFCZ5_WT_filtered <- FindNeighbors(QFCZ5_WT_filtered, dims = 1:10)
QFCZ5_WT_filtered <- FindClusters(QFCZ5_WT_filtered, resolution = 0.5)

# UMAP
QFCZ5_WT_filtered <- RunUMAP(QFCZ5_WT_filtered, dims = 1:10)
DimPlot(QFCZ5_WT_filtered, reduction = "umap")

#nb of gene per clusters
table(QFCZ5_WT_filtered@meta.data$seurat_clusters)

saveRDS(QFCZ5_WT_filtered, file = "../output/QFCZ5_WT_filtered.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
QFCZ.markers <- FindAllMarkers(QFCZ5_WT_filtered, only.pos = TRUE)
QFCZ.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
QFCZ.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(QFCZ5_WT_filtered, features = top10$gene) + NoLegend()




# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)




# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/demo/singleCell_singleR_part2/scripts")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
DimPlot(pbmc.seurat.filtered, reduction = "umap")
View(pbmc.seurat.filtered@meta.data)


# run SingleR with multiple reference datasets (default mode) ---------

# for pbmc data, we will use two datasets
hpca <- celldex::HumanPrimaryCellAtlasData()
dice <- celldex::DatabaseImmuneCellExpressionData()

# ...1. Strategy 1: Using reference-specific labels ----------
hpca$label.main
dice$label.main

# adding ref info to labels
hpca$label.main <- paste0('HPCA.', hpca$label.main)
dice$label.main <- paste0('DICE.', dice$label.main)

# create a combined ref based on shared genes
shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])
combined
combined$label.main

# run singleR using combined ref
# savings counts into a separate object
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels)

pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), 'labels']
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)

# ...2. Strategy 2: Comparing scores across references ----------

hpca$label.main
dice$label.main
hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
dice$label.main <- gsub('DICE\\.','', dice$label.main)

com.res2 <- SingleR(test = pbmc_counts, 
        ref = list(HPCA = hpca, DICE = dice),
        labels = list(hpca$label.main, dice$label.main))

# Check the final label from the combined assignment.
table(com.res2$labels)

# which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.', com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))

# get de. genes from each individual references
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes

# Combined diagnostics
plotScoreHeatmap(com.res2)


# ...3. Strategy 3: Using Harmonized Labels ----------

hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same sets of genes:
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms:
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)), 10)

# using label.ont instead on label.main while running SingleR

com.res3 <- SingleR(test = pbmc_counts,
        ref = list(HPCA = hpca.ont, DICE = dice.ont),
        labels = list(hpca.ont$label.ont, dice.ont$label.ont))


table(com.res3$labels)



# How to map cell ontology terms? ----------------

colData(hpca.ont)
colData(dice.ont)

hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.fle, header = F)