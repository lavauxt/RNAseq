#' Basic Find Markers
#'
#' To find de `novo` markers by `FindAllMarkers` from Seurat with default setting.
#'
#' @param sce An object of class Seurat
#' @param group default:seurat_clusters, you can change it to  celltype
#' @param dir path for saving results
#'
#' @return sce.markers a data.frame of markers.
#'
#' @import Seurat
#' @importFrom ggplot2 ggsave coord_flip
#' @importFrom utils write.csv
#' @importFrom dplyr group_by top_n
#' @export
#'
#' @examples
#' \donttest{
#' basic_find_markers(AJ064_small_last_sce,dir=tempdir())
#' }
basic_find_markers <-  function(sce,group='seurat_clusters',dir='.'){

  Seurat::Idents(sce)= sce@meta.data[,group]
  message( table(Seurat::Idents(sce)) )
  n=length(table(Seurat::Idents(sce)))
  # library(future)
  # plan("multiprocess", workers = 8)
  sce.markers <- Seurat::FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
                                thresh.use = 0.25)
  utils::write.csv(sce.markers,
            file=file.path(dir,paste0(group,'_sce.markers.csv')))

  top10 <- sce.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)
  Seurat::DoHeatmap(sce,top10$gene,size=3)
  ggplot2::ggsave(filename=file.path(dir,paste0(group,'_sce.markers_check_top10_heatmap.pdf')),
          height = 15,width = 18)

  p <- Seurat::DotPlot(sce, features = unique(top10$gene),
               assay='RNA'  )  + coord_flip()


  ggplot2::ggsave(file.path(dir,paste0(group,'_DotPlot_check_top10_markers_by_clusters.pdf')),
         height = 18)

  top3 <- sce.markers %>% dplyr::group_by(cluster) %>% top_n(3, avg_log2FC)
  DoHeatmap(sce,top3$gene,size=3)
  ggplot2::ggsave(file.path(dir,paste0(group,'_DoHeatmap_check_top3_markers_by_clusters.pdf')))


  p <- DotPlot(sce, features = unique(top3$gene),
               assay='RNA'  )  + coord_flip()


  ggsave(file.path(dir,paste0(group,'_DotPlot_check_top3_markers_by_clusters.pdf')),height = 8)
  save(sce.markers,file = file.path(dir,paste0(group,'_sce.markers.Rdata')))
  return(sce.markers)

}
