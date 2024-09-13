#' Basic Quality Control
#'
#' add `percent_mito`,`percent_ribo`,`percent_hb` to the Seurat class.
#' And draw `VlnPlot` for these `qc` values.
#'
#' @param sce An object of class Seurat
#' @param org human or mouse, default: human
#' @param group default:`orig.ident`,you can change it to `seurat_clusters` or `celltype`
#' @param dir the path for saving the figures by `DotPlot` with known famous markers.
#'
#' @return list(p1,p2,p3,sce), the last one in the new `sce`.
#' @import ggplot2
#' @import Seurat
#' @importFrom Seurat VlnPlot NoLegend FeatureScatter
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' basic_qc(AJ064_small_sce,dir= tempdir())
#'
basic_qc <- function(sce,org='human',group='orig.ident',dir='.'){
 if(org=='human'){
   sce=PercentageFeatureSet(sce, "^MT-", col.name = "percent_mito")
   sce=PercentageFeatureSet(sce, "^RP[SL]", col.name = "percent_ribo")
   sce=PercentageFeatureSet(sce, "^HB[^(P)]", col.name = "percent_hb")
 }else if(org=='mouse'){
   sce=PercentageFeatureSet(sce, "^Mt-", col.name = "percent_mito")
   sce=PercentageFeatureSet(sce, "^Rp[sl]", col.name = "percent_ribo")
   sce=PercentageFeatureSet(sce, "^Hb[^(p)]", col.name = "percent_hb")
 }else{
   stop('So far, we only accept human and mouse')
 }

  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(sce, group.by =  group ,
             features = feats,
             pt.size = 0.01, ncol = 2) +
    NoLegend()
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(sce, group.by =  group ,
             features = feats, pt.size = 0.01,
             ncol = 3, same.y.lims=T) +
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2

  p3=FeatureScatter(sce, "nCount_RNA", "nFeature_RNA",
                    group.by = group , pt.size = 0.5)
  
  pdf(file.path(dir,
                paste0( 'basic_qc_based_on_',group,'.pdf')))
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  return( sce )

}
