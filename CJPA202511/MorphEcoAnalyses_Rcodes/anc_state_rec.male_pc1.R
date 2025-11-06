library(ape)
library(phytools)
library(geiger)

tree <- read.tree("G:/longxi/longxibioinfo/PhenotypeStatistics/PCA_anc/diplo_timetree.nwk")

phenotype_data <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCA_anc/PCA_Male_MeanPC1.csv")

phenotype_data$species <- trimws(phenotype_data$species)
tree$tip.label <- trimws(tree$tip.label)

common_species <- intersect(tree$tip.label, phenotype_data$species)

tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
phenotype_data <- phenotype_data[phenotype_data$species %in% common_species, ]
phenotype_data <- phenotype_data[match(tree$tip.label, phenotype_data$species), ]

svl_vector <- phenotype_data$mean_PC1
names(svl_vector) <- phenotype_data$species

if (any(is.na(svl_vector))) {
  print("Missing species PC1 data:")
  print(phenotype_data[is.na(svl_vector), ])
  stop("Check data file, fill in missing PC1 values.")
}

tree_ultra <- chronos(tree)

bm_model <- fitContinuous(tree_ultra, svl_vector, model = "BM")
ou_model <- try(fitContinuous(tree_ultra, svl_vector, model = "OU"), silent = TRUE)
white_model <- fitContinuous(tree_ultra, svl_vector, model = "white")

cat("BM Model:\n")
print(bm_model)

cat("\nWhite Model:\n")
print(white_model)

if (!inherits(ou_model, "try-error")) {
  cat("\nOU Model:\n")
  print(ou_model)
  
  bm_aicc <- bm_model$opt$aicc
  ou_aicc <- ou_model$opt$aicc
  white_aicc <- white_model$opt$aicc
  
  best_aicc <- min(bm_aicc, ou_aicc, white_aicc)
  
  if (best_aicc == bm_aicc) {
    best_model <- "BM"
  } else if (best_aicc == ou_aicc) {
    best_model <- "OU"
  } else {
    best_model <- "White"
  }
  
  cat("\nBest model:", best_model, ", AICc =", best_aicc, "\n")
  
} else {
  cat("\nOU model fitting failed, only comparing BM and White models\n")
  
  bm_aicc <- bm_model$opt$aicc
  white_aicc <- white_model$opt$aicc
  
  best_model <- ifelse(bm_aicc < white_aicc, "BM", "White")
  best_aicc <- min(bm_aicc, white_aicc)
  
  cat("\nBest model:", best_model, ", AICc =", best_aicc, "\n")
}

anc_state <- fastAnc(tree_ultra, svl_vector)

cont_map <- contMap(tree_ultra, svl_vector, plot = FALSE)

custom_colors <- colorRampPalette(c(
  "#C00000",
  "#ED6547",
  "#FFD966",
  "#767171",
  "#767171",
  "#767171",
  "#95B3D4",
  "#00B0F0",
  "#24A1C2"
))(100)

cont_map <- setMap(cont_map, custom_colors)

plot(cont_map,
     fsize = 0.9,
     outline = TRUE,
     legend = 0.4 * max(svl_vector),
     lwd = 5,
     main = "Ancestral Trait Evolution Map (Custom Colors)")
     
nodelabels(text = round(anc_state, 2), 
           frame = "none",
           cex = 0.8,
           adj = c(-0.3, 0.2))