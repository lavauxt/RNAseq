#' Basic Markers
#'
#' @param sce An object of class Seurat
#' @param org human or mouse, default: human
#' @param group default:`orig.ident`, you can change it to `seurat_clusters` or `celltype`
#' @param dir the path for saving the figures by `DotPlot` with known famous markers.
#'
#' @return a list of figures by `DotPlot`
#' @import ggplot2
#' @import Seurat
#' @import clustree
#' @importFrom stringr str_to_title
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' \donttest{
#' basic_markers(AJ064_small_last_sce,dir=tempdir())
#' }
basic_markers <-  function(sce,org='human',group='orig.ident',dir='.'){

  # T Cells (CD3D, CD3E, CD8A),
  # B cells (CD19, CD79A, MS4A1 [CD20]),
  # Plasma cells (IGHG1, MZB1, SDC1, CD79A),
  # Monocytes and macrophages (CD68, CD163, CD14),
  # NK Cells (FGFBP2, FCG3RA, CX3CR1),
  # Photoreceptor cells (RCVRN),
  # Fibroblasts (FGF7, MME),
  # Endothelial cells (PECAM1, VWF).
  # epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
  #   immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM),
  # stromal (CD10+,MME,fibo or CD31+,PECAM1,endo)

  all_markers =c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                'CD19', 'CD79A', 'MS4A1' ,
                'IGHG1', 'MZB1', 'SDC1',
                'CD68', 'CD163', 'CD14',
                'TPSAB1' , 'TPSB2',  # mast cells,
                'RCVRN','FPR1' , 'ITGAM' ,
                'C1QA',  'C1QB',  # mac
                'S100A9', 'S100A8', 'MMP19',# monocyte
                'LAMP3', 'IDO1','IDO2',## DC3
                'CD1E','CD1C', # DC2
                'KLRB1','NCR1', # NK
                'FGF7','MME', 'ACTA2',
                'PECAM1', 'VWF',
                'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )



  # TH17 cells (C-12, CCR6+ and RORC)
  # type 1 helper T cells (TH1; C-17, CXCR3+ and IFNG and TBX21)
  # heterogeneous continuum of intermediate TH1/TH17 states (C-13, C-16 and C-19)
  # with varying degrees of CXCR3, CCR6, CCR5 and CD161 surface protein expression and RORC and TBX21 expression
  # CD161+ subset of type 2 helper T (TH2) cells (C-14),
  # described as pathogenic, with higher expression of allergy-associated HPGDS and IL17RB
  # naive (LEF1, SELL, TCF7),
  # effector (IFNG),
  # cytotoxicity (GZMB, PRF1),
  # early and general exhaustion (PDCD1, CTLA4, ENTPD1 ) .
  # antigen presentation (CD74, HLA-DRB1/5, HLA-DQA2)
  # FCGR3A (FCGR3A+ Monocyte), KLRF1 (NK), FCER1A (DC), and PF4 (MP/platelets).

  Tcells_markers = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                     'LEF1', 'SELL' , 'TCF7', # naive
                     'FOXP3',
                     'CCR6', 'RORC' ,  # TH17 cells
                     'TBX21', 'CXCR3', 'IFNG', # type 1 helper T cells
                     'CCR3','CCR4',
                     'PDCD1',  'CTLA4','ENTPD1', # early and general exhaustion
                     'GZMB', 'GZMK','PRF1', # cytotoxicity
                     'IFNG',   'CCL3' ,'CXCR6' , 'ITGA1',  # effector
                     'NKG7','KLRF1','MKI67','PF4','FCER1A','FCGR3A')

  # mast cells, TPSAB1 and TPSB2
  # B cell,  CD79A  and MS4A1 (CD20)
  # naive B cells, such as MS4A1 (CD20), CD19, CD22, TCL1A, and CD83,
  # plasma B cells, such as CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4
  Bcells_markers = c('CD3D','MS4A1','CD79A',
                     'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
                     'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
                     'TPSAB1' , 'TPSB2',  # mast cells,
                     'PTPRC' )


  Myeloid_markers =c('CD68', 'CD163', 'CD14',  'CD86','C1QA',  'C1QB',  # mac
                    'S100A9', 'S100A8', 'MMP19',# monocyte
                    'LAMP3', 'IDO1','IDO2',## DC3
                    'MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7',
                    'CD1E','CD1C', # DC2
                    'XCR1','CLEC9A','FCER1A',# DC1
                    'GZMB','TCF4','IRF7')

  # epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
  # - alveolar type I cell (AT1; AGER+)
  # - alveolar type II cell (AT2; SFTPA1)
  # - secretory club cell (Club; SCGB1A1+)
  # - basal airway epithelial cells (Basal; KRT17+)
  # - ciliated airway epithelial cells (Ciliated; TPPP3+)

  epi_markers = c(  'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' ,
                       'AGER','SFTPA1','SCGB1A1','KRT17','TPPP3',
                       'KRT4','KRT14','KRT8','KRT18',
                       'CD3D','PTPRC' )


  stromal_markers = c('TEK',"PTPRC","EPCAM","PDPN","PECAM1",'PDGFRB',
                     'CSPG4','GJB2', 'RGS5','ITGA7',
                     'ACTA2','RBP1','CD36', 'ADGRE5','COL11A1','FGF7', 'MME')



  # basal (e.g., KRT5, ACTA2, MYLK, SNAI2),
  # luminal progenitor (TNFRSF11A (RANK), KIT),
  # mature luminal cells (ESR1, PGR, FOXA1)

  brca_markers = c(  'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' ,
                       'KRT5', 'ACTA2', 'MYLK', 'SNAI2', # basal
                       'RANK', 'KIT', # luminal progenitor
                       'ESR1', 'PGR', 'FOXA1',# mature luminal cells
                       'KRT4','KRT14','KRT8','KRT18',
                       'CD3D','PTPRC' )
  if(org=='human'){
  }else if(org=='mouse'){
    all_markers=str_to_title(all_markers)
    Tcells_markers=str_to_title(Tcells_markers)
    Bcells_markers=str_to_title(Bcells_markers)
    Myeloid_markers=str_to_title(Myeloid_markers)
    epi_markers=str_to_title(epi_markers)
    stromal_markers=str_to_title(stromal_markers)
    brca_markers=str_to_title(brca_markers)
  } else {
    stop('So far, we only accept human and mouse')
  }
  genes_to_check=list(
    all_markers=unique(all_markers),
    Tcells_markers=unique(Tcells_markers),
    Bcells_markers=unique(Bcells_markers),
    Myeloid_markers=unique(Myeloid_markers),
    epi_markers=unique(epi_markers),
    stromal_markers=unique(stromal_markers),
    brca_markers=unique(brca_markers)
  )
  # dpl: dotplot, list
  dpl <- lapply(genes_to_check, function(cg){
    DotPlot(sce, assay = "RNA",
            features = cg,
            group.by =  group  

      ) + coord_flip()
  })
  
  pdf(file.path(dir,
                paste0( 'markers_based_on_',group,'.pdf')))
  lapply(1:length(genes_to_check), function(i){
    p=names(genes_to_check)[i] 
    print(dpl[[i]]+ggtitle(p))
  })
  dev.off()
  return(dpl)

}


