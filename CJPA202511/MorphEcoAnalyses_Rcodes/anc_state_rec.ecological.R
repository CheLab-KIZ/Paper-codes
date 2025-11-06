library(ape)
library(phytools)
library(geiger)

tree <- read.tree("G:\\longxi\\longxibioinfo\\distrib\\diplo_timetree.nwk")
phenotype_data <- read.csv("G:\\longxi\\longxibioinfo\\distrib\\eco_recon\\species_environment_factors_mean.csv")

phenotype_data <- phenotype_data[match(tree$tip.label, phenotype_data$species), ]

env_factors <- c("bio_1", "bio_12", "bio_4","bio_15", "mean_NPP", "Elevation")

tree_ultra <- chronos(tree)

for (factor in env_factors) {
  cat("\nProcessing environmental factor:", factor, "\n")
  
  factor_vector <- phenotype_data[[factor]]
  names(factor_vector) <- tree$tip.label
  
  bm_model <- fitContinuous(tree_ultra, factor_vector, model = "BM")
  ou_model <- fitContinuous(tree_ultra, factor_vector, model = "OU")
  white_model <- fitContinuous(tree_ultra, factor_vector, model = "white")
  
  cat("BM Model Result:\n")
  print(bm_model)
  
  cat("OU Model Result:\n")
  print(ou_model)
  
  cat("White Model Result:\n")
  print(white_model)
  
  bm_aicc <- bm_model$opt$aicc
  ou_aicc <- ou_model$opt$aicc
  white_aicc <- white_model$opt$aicc
  
  cat("BM Model AICc:", bm_aicc, "\n")
  cat("OU Model AICc:", ou_aicc, "\n")
  cat("White Model AICc:", white_aicc, "\n")
  
  best_aicc <- min(bm_aicc, ou_aicc, white_aicc)
  
  if (best_aicc == bm_aicc) {
    best_model <- "BM"
  } else if (best_aicc == ou_aicc) {
    best_model <- "OU"
  } else {
    best_model <- "White"
  }
  
  cat("Best model is:", best_model, "AICc value:", best_aicc, "\n")
  
  tryCatch({
    anc_state <- fastAnc(tree_ultra, factor_vector)
    
    cont_map <- contMap(tree_ultra, factor_vector, plot = FALSE, res = 200)
    
    pdf(file = paste0("Ancestor_", factor, "_Tree.pdf"), width = 12, height = 8)
    
    plotTree(tree_ultra, type = "phylogram", fsize = 0.8, main = paste("Ancestor", factor, "Estimates (ML)"))
    
    plot(cont_map, add = TRUE)
    
    nodelabels(round(anc_state, 2), bg = "lightblue", cex = 0.7)
    
    dev.off()
    
    cat("Saved ancestor trait tree for environmental factor", factor, "as: Ancestor_", factor, "_Tree.pdf\n")
  }, error = function(e) {
    cat("Ancestor trait reconstruction failed:", e$message, "\n")
  })
}