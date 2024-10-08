% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_find_markers.R
\name{basic_find_markers}
\alias{basic_find_markers}
\title{Basic Find Markers}
\usage{
basic_find_markers(sce, group = "seurat_clusters", dir = ".")
}
\arguments{
\item{sce}{An object of class Seurat}

\item{group}{default:seurat_clusters, you can change it to  celltype}

\item{dir}{path for saving results}
}
\value{
sce.markers a data.frame of markers.
}
\description{
To find de `novo` markers by `FindAllMarkers` from Seurat with default setting.
}
\examples{
\donttest{
basic_find_markers(AJ064_small_last_sce,dir=tempdir())
}
}
