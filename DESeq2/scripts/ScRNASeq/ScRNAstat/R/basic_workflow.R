#' Basic Workflow
#'
#' the workflow from Seurat, including: `NormalizeData`,`FindVariableFeatures`,`ScaleData`,
#' `RunPCA`,`RunTSNE`,`RunUMAP`,`FindNeighbors`,`FindClusters(sce, resolution = seq(0.1,1,by=0.1))`
#' we use `clustree` to check the different resolution for `FindClusters`.
#'
#' @param sce  An object of class Seurat
#' @param dir the path for saving the figures by `DotPlot` with known famous markers.
#'
#' @return list(p1,p2,p3,sce), the last one in the new sce with PCA,tSNE,UMAP information.
#' @import ggplot2
#' @import Seurat
#' @import clustree
#' @import patchwork
#' @importFrom Seurat NormalizeData GetAssay FindVariableFeatures ScaleData RunPCA VariableFeatures RunTSNE RunUMAP FindNeighbors FindClusters DimPlot
#' @importFrom clustree clustree
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' \dontrun{
#' basic_workflow(AJ064_small_sce,dir=tempdir())
#' }
basic_workflow <- function(sce,dir='.'){
  sce
  sce <- NormalizeData(sce,
                       normalization.method =  "LogNormalize",
                       scale.factor = 10000)
  GetAssay(sce,assay = "RNA")
  sce <- FindVariableFeatures(sce,
                              selection.method = "vst",
                              nfeatures = 2000)
  sce <- ScaleData(sce)
  sce <- RunPCA(object = sce, npcs = 10,
                pc.genes = VariableFeatures(sce))

  sce <- RunTSNE(object = sce,
                 dims = 1:10,
                 do.fast = TRUE)
  sce = RunUMAP(sce, dims = 1:10)
  sce <- FindNeighbors(sce, dims = 1:10)
  sce <- FindClusters(sce, resolution = seq(0.1,1,by=0.1))
  p1=clustree(sce)
  
  colnames(sce@meta.data)

  cg=colnames(sce@meta.data)[grepl('RNA_snn', colnames(sce@meta.data))]
  pl=lapply(cg, function(i){
    DimPlot(sce,reduction = "umap",
            group.by = i,
            label=T)
  })
 
  pdf(file.path(dir,
                paste0( 'basic_workflow_clustree.pdf')))
  print(p1)
  dev.off()
  
  pdf(file.path(dir,
                paste0( 'basic_workflow_check_resolution.pdf')))
  lapply(pl, print)
  dev.off()
   
  
  return(sce)

}
