# how to use this function based on Seurat CellCycleScoring()

# source("/path/to/your/SexScore.R")

# Small list of sex related gene
# male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1")
# female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Utx", "Rps4x", "Jarid1c", "Prdx4")

# Extended list of sex related gene
# male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1", 
#                   "Smcy", "Eif1ay", "Ssty1", "Ssty2", "Rbmy1a1", "Sox30", "Sox9", "Amh", "Insl3", 
#                   "Hsd17b3", "Cyp17a1")

# female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Rps4x", "Jarid1c", "Prdx4", "Med14",
#                     "Zfx", "Ocrl", "Ddx39b", "Gpm6b", "Fthl17", "Esr1", "Foxl2", "Cyp19a1", "Wnt4", 
#                     "Amhr2", "Eda2r")

# male.features <- c(
#  "Abca13", "Adam22", "Arhgap31-as1", "Bcorp1", "Bpi", "Camp", "Cdy4p", "Ceacam19", 
#  "Ceacam6", "Ceacam8", "Celsr2", "Crisp3", "Ctd-2308n23.2", "Ddx3y", "Defa3", "Defa4", 
#  "Eif1ay", "Fbn1", "Gstt1", "Inhba", "Kalp", "Kdm5d", "Kif9-as1", "Lcn2", "Linc00278", 
#  "Ltf", "Mmp8", "Mpo", "Prky", "Rn7skp282", "Rp11-424g14.1", "Rp6-99m1.2", "Rps4y1", 
#  "Slc1a2", "Tmsb4y", "Ttty14", "Ttty15", "Txlng2p", "Uty", "Usp9y", "Xist", "Ythdf3-as1", 
#  "Zfy", "Zfy-as1", "Znf839p1"
#)

# female.features <- c(
#  "Ac084018.1", "Ighg3", "Ighv1-46", "Igkv1-27", "Igkv3-15", "Igkv3d-11", "Iglv1-40", 
#  "Iglv1-47", "Iglv3-1", "Iglv6-57", "Ift81", "Kdm5c", "Prkx", "Tsix", "Xist"
#)

#object <- CellCycleScoring(
#  object = object,
#  s.features = female.features,  # Female-specific genes
#  g2m.features = male.features,  # Male-specific genes
#  ctrl = NULL,  # Keep default control behavior
#  set.ident = TRUE  # Set cell identities to the identified sex
#)

# Then we can regress the Seurat object
# object$Sex.Difference <- object$Female.Score - object$Male.Score 
# object <- ScaleData( object, vars.to.regress = "Sex.Difference", features = rownames(object) 


SexScoring <- function(
  object,
  female.features,
  male.features,
  ctrl = NULL,
  set.ident = FALSE,
  ...
) {

  name <- 'Sex.Score'
  features <- list('Female.Score' = female.features, 'Male.Score' = male.features)
  
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  
  object.sex <- AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl,
    ...
  )
  
  sex.columns <- grep(pattern = name, x = colnames(x = object.sex[[]]), value = TRUE)
  sex.scores <- object.sex[[sex.columns]]
  
  rm(object.sex)
  CheckGC()
  
  # Assign cells to 'Female', 'Male', or 'Undetermined' based on scores
  assignments <- apply(
    X = sex.scores,
    MARGIN = 1,
    FUN = function(scores, female = 'Female', male = 'Male', null = 'Undetermined') {
      if (all(scores < 0)) {
        return(null)  # Assign to 'Undetermined' if both scores are low (optional)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undetermined')  # Optional: 'Undetermined' for ties
        } else {
          return(c(female, male)[which(x = scores == max(scores))])
        }
      }
    }
  )
  
  # Combine scores and assignments
  sex.scores <- merge(x = sex.scores, y = data.frame(assignments), by = 0)
  colnames(x = sex.scores) <- c('rownames', 'Female.Score', 'Male.Score', 'Sex')
  rownames(x = sex.scores) <- sex.scores$rownames
  sex.scores <- sex.scores[, c('Female.Score', 'Male.Score', 'Sex')]
  
  object[[colnames(x = sex.scores)]] <- sex.scores
  
  # Optionally, set cell identity to 'Sex'
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Sex'
  }
  
  return(object)
}