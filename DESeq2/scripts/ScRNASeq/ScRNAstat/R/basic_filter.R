#' basic_filter
#'
#' filter the genes which show expression less than 3 cells.
#' filter the cells which percent_mito < 25 & percent_ribo > 3 & percent_hb < 10
#' filter the cells which  nFeature_RNA > 300 & nFeature_RNA < 8000
#'
#' @param sce  An object of class Seurat
#'
#' @return sce.all.filt  An object of class Seurat
#' @import Seurat
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
#' basic_filter(AJ064_small_sce)
#'
basic_filter <-  function(sce){
  if (!'percent_mito' %in% colnames(sce@meta.data)){
    sce=basic_qc(sce,dir = tempdir())
  }
  selected_f <- rownames(sce)[Matrix::rowSums(sce@assays$RNA@counts > 0 ) > 3]
  sce.all.filt <- subset(sce, features = selected_f)
  sce.all.filt
  sce.all.filt <- subset(sce.all.filt,
                         percent_mito < 25 & percent_ribo > 3 & percent_hb < 10)

  sce.all.filt <- subset(sce.all.filt,
                         nFeature_RNA > 300 & nFeature_RNA < 8000  )
  sce.all.filt

  return(sce.all.filt)
}
