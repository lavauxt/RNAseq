# how to use this function based on Seurat CellCycleScoring()

# source("/path/to/your/SexScore.R")

# Seuratobject_sex <- SexScoring(
# object = object, # after normalization
# female.features = female.features, 
# male.features = male.features, 
# ctrl = NULL, 
# set.ident = TRUE 
# )

# Use for regression
# Seuratobject_sex$Sex.Difference <-Seuratobject_se$Female.Score - Seuratobject_se$Male.Score 
# Seuratobject_sex_final <- ScaleData(Seuratobject_sex, vars.to.regress = "Sex.Difference", features = rownames(Seuratobject_se))

# display PCA
# PCA_test_final <- RunPCA(Seuratobject_sex_final, features = c(male.features, female.features))
# DimPlot(PCA_test_final)


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