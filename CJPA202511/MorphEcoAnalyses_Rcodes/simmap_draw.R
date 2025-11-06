setwd("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis")
library(phytools)

tree <- read.tree("G:/longxi/longxibiaoxingyanhuan/zhuxianxingzhuangchongjian/diplo_timetree.nwk")

factors <- read.delim("SVL_PCA_ENV.txt", header = TRUE, row.names = 1, check.names = FALSE)

tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(factors)))

habit <- as.factor(setNames(factors$habitat, rownames(factors)))
levels(habit) <- c("No-Radiation", "Radiation")
cols <- setNames(c("#24A1C2", "#BE213B"), levels(habit))

habit_trees <- make.simmap(tree, habit, nsim = 100)

env_vars <- c("NPP", "bio_1", "bio_4", "bio_12", "bio_15")

for (var in env_vars) {
  cat("Processing:", var, "\n")
  
  var_data <- setNames(factors[[var]], rownames(factors))
  
  anc_var <- fastAnc(tree, var_data)
  var_all <- c(anc_var, var_data)
  
  wilcox_result <- wilcox.test(factors[[var]] ~ factors$habitat, exact = FALSE)
  z_val <- round(qnorm(wilcox_result$statistic / (length(var_data) * (length(var_data) + 1) / 2)), 3)
  p_val <- signif(wilcox_result$p.value, 3)
  
  phenogram(habit_trees[[1]], var_all,
            colors = cols,
            lwd = 3,
            spread.labels = TRUE,
            spread.cost = c(1, 0),
            fsize = 1.1,
            ftype = "off",
            ylab = var,
            main = paste0(var, " Traitgram by Habitat"))
  
  legend("topleft", legend = levels(habit),
         col = cols, lwd = 3, title = "Habitat", bty = "n")
  
  text(x = 5, y = max(var_all, na.rm = TRUE) * 0.95,
       labels = paste0("Wilcoxon test:\nZ = ", z_val, ", p = ", p_val),
       adj = 0, cex = 1.1)
  
  box(lwd = 2)
}