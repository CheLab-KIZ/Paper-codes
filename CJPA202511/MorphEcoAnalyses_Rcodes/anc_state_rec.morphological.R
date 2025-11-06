library(ape)
library(phytools)
library(geiger)

tree <- read.tree("G:/longxi/longxibioinfo/PhenotypeStatistics/PCA_rec/diplo_timetree.nwk")
phenotype_data <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCA_rec/input.csv")

phenotype_data$species <- trimws(phenotype_data$species)
tree$tip.label <- trimws(tree$tip.label)

traits <- grep("_mean$", names(phenotype_data), value = TRUE)
traits <- traits[traits != "SVL_mean"]

custom_colors <- colorRampPalette(c(
  "#C00000", "#ED6547", "#FFD966", "#767171",
  "#767171", "#767171", "#95B3D4", "#00B0F0", "#24A1C2"
))(100)

plot_ancestral_for_sex <- function(sex_data, sex_label) {
  common_species <- intersect(tree$tip.label, sex_data$species)
  tree_sex <- drop.tip(tree, setdiff(tree$tip.label, common_species))
  sex_data <- sex_data[sex_data$species %in% common_species, ]
  sex_data <- sex_data[match(tree_sex$tip.label, sex_data$species), ]
  tree_ultra <- chronos(tree_sex)
  
  par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
  
  for (trait in traits) {
    trait_vector <- sex_data[[trait]]
    names(trait_vector) <- sex_data$species
    
    if (any(is.na(trait_vector))) {
      plot.new()
      text(0.5, 0.5, trait, cex = 2, font = 2)
      next
    }
    
    anc_state <- fastAnc(tree_ultra, trait_vector)
    cont_map <- contMap(tree_ultra, trait_vector, plot = FALSE)
    cont_map <- setMap(cont_map, custom_colors)
    
    plot(cont_map,
         fsize = 0.8,
         outline = TRUE,
         legend = 0.4 * max(trait_vector),
         lwd = 4)
    
    usr <- par("usr")
    text(mean(usr[1:2]), mean(usr[3:4]), 
         labels = gsub("_mean", "", trait), 
         cex = 3, font = 2, col = rgb(0, 0, 0, 0.3))
    
    nodelabels(text = round(anc_state, 2),
               frame = "none", cex = 0.6, adj = c(-0.3, 0.2))
  }
}

pdf("female_traits_centerlabel.pdf", width = 12, height = 10)
female_data <- subset(phenotype_data, SEX == "F")
plot_ancestral_for_sex(female_data, "Female")
dev.off()

pdf("male_traits_centerlabel.pdf", width = 12, height = 10)
male_data <- subset(phenotype_data, SEX == "M")
plot_ancestral_for_sex(male_data, "Male")
dev.off()

cat("Generated female_traits_centerlabel.pdf and male_traits_centerlabel.pdf\n")