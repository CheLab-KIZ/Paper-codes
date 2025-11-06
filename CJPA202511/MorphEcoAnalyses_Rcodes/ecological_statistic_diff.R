library(car)
library(ggplot2)
library(ggpubr)

out_dir <- "G:/longxi/longxibioinfo/distrib/eco_recon/ecorec-20250809"

df <- read.csv(file.path(out_dir, "Species_EnvironmentalFactorsMean.csv"),
               header = TRUE, fileEncoding = "UTF-8-BOM")

env_vars <- c("mean_NPP", "bio_1", "bio_12", "bio_15", "bio_4", "Elevation")

results <- data.frame()
norm_results <- data.frame()
var_results <- data.frame()

plot_list <- list()

for (var in env_vars) {
  norm_test <- by(df[[var]], df$Clade, shapiro.test)
  norm_df <- data.frame(
    Variable = var,
    Clade = names(norm_test),
    p_value = sapply(norm_test, function(x) x$p.value)
  )
  norm_results <- rbind(norm_results, norm_df)
  
  var_test <- leveneTest(df[[var]] ~ as.factor(df$Clade))
  var_results <- rbind(var_results,
                       data.frame(Variable = var,
                                  Test = "Levene's",
                                  p_value = var_test$`Pr(>F)`[1]))
  
  pvals_norm <- norm_df$p_value
  if (all(pvals_norm > 0.05)) {
    test <- t.test(df[[var]] ~ df$Clade)
    method <- "t-test"
  } else {
    test <- wilcox.test(df[[var]] ~ df$Clade)
    method <- "Wilcoxon"
  }
  
  results <- rbind(results, data.frame(
    Variable = var,
    p_value = round(test$p.value, 4),
    Method = method,
    Significance = ifelse(test$p.value < 0.05, "YES", "NO")
  ))
  
  p <- ggplot(df, aes(x = Clade, y = .data[[var]], fill = Clade)) +
    geom_boxplot(width = 0.5, alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
    stat_compare_means(method = ifelse(method == "t-test", "t.test", "wilcox.test")) +
    scale_fill_manual(values = c("Other clade" = "#24A1C2", "Radiation" = "#BD213B")) +
    labs(title = paste(var, "by Clade"), y = var) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(out_dir, paste0(var, "_by_clade.png")),
         plot = p, width = 5, height = 4)
  
  plot_list[[var]] <- p
}

combined_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 3)

print(combined_plot)

ggsave(file.path(out_dir, "combined_plot.png"),
       plot = combined_plot, width = 10, height = 8)

write.csv(results, file.path(out_dir, "clade_comparison_results.csv"), row.names = FALSE)
write.csv(norm_results, file.path(out_dir, "normality_test_results.csv"), row.names = FALSE)
write.csv(var_results, file.path(out_dir, "variance_test_results.csv"), row.names = FALSE)