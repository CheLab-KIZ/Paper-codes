library(FactoMineR)
library(ggplot2)
library(dplyr)

data <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/all_data_imputed.csv")

output_dir <- "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/"

run_pca_analysis <- function(data, label, output_dir) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    data_pca <- data[, c("SVL", "TL", "TRL", "HW", "HL", "SEL", "HH", "FLL", "HLL", "T4L")]
    
    data_no_na <- na.omit(data)
    data_pca_no_na <- na.omit(data_pca)
    stopifnot(nrow(data_no_na) == nrow(data_pca_no_na))
    
    pca_res <- PCA(data_pca_no_na, scale.unit = TRUE, ncp = 5, graph = FALSE)
    
    pc_percent <- round(pca_res$eig[, 2], 1)
    
    write.csv(as.data.frame(pca_res$eig), 
              file.path(output_dir, paste0("PCA_", label, "_ExplainedVariance.csv")), row.names = FALSE)
    write.csv(as.data.frame(pca_res$var$coord), 
              file.path(output_dir, paste0("PCA_", label, "_Loadings.csv")), row.names = TRUE)
    
    ind_scores <- as.data.frame(pca_res$ind$coord)
    ind_scores <- cbind(data_no_na[, c("species", "zhixi", "SEX")], ind_scores)
    write.csv(ind_scores, 
              file.path(output_dir, paste0("PCA_", label, "_Scores.csv")), row.names = FALSE)
    
    mean_pc1 <- ind_scores %>%
      group_by(species) %>%
      summarise(mean_PC1 = mean(Dim.1, na.rm = TRUE))
    write.csv(mean_pc1, 
              file.path(output_dir, paste0("PCA_", label, "_MeanPC1_by_species.csv")), row.names = FALSE)
    
    my_colors <- c("Radiation" = "#BD213B", "no-Radiation" = "#24A1C2")
    
    p <- ggplot(ind_scores, aes(x = Dim.1, y = Dim.2, color = zhixi)) +
      geom_point(size = 3, alpha = 0.85) +
      scale_color_manual(values = my_colors) +
      theme_minimal() +
      labs(
        title = paste("PCA (", label, ")"),
        x = paste0("PC1 (", pc_percent[1], "%)"),
        y = paste0("PC2 (", pc_percent[2], "%)"),
        color = "Ecomorph"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2)
      )
    
    ggsave(file.path(output_dir, paste0("PCA_", label, "_PC1_PC2_plot.pdf")),
           plot = p, width = 6, height = 7.5)
    
    print(p)
}

run_pca_analysis(data, "Combined", output_dir)
run_pca_analysis(subset(data, SEX == "F"), "Female", output_dir)
run_pca_analysis(subset(data, SEX == "M"), "Male", output_dir)

library(dplyr)
library(ggplot2)
library(ggpubr)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PCA_Combined_Scores.csv")

male_scores <- scores %>% filter(SEX == "M")

if (length(unique(male_scores$zhixi)) == 2) {
  t_test_male <- t.test(Dim.1 ~ zhixi, data = male_scores)
  print(t_test_male)
  
  capture.output(t_test_male, file = "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PC1_t_test_male_by_zhixi.txt")
}

p <- ggplot(male_scores, aes(x = zhixi, y = Dim.1, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(method = "t.test", label = "p.format", size = 14,
                     label.y = max(male_scores$Dim.1) * 1.05,
                     label.x = 1.5) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "Male Lizard PC1 Score Comparison by Ecomorph",
    x = "Ecomorph (zhixi)",
    y = "PC1 Score (Body Size Index)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

print(p)

ggsave("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_PC1_male_by_zhixi.pdf", 
       plot = p, width = 7, height = 9)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PCA_Combined_Scores.csv")

male_scores <- scores %>% filter(SEX == "F")

group_var <- "zhixi"
value_var <- "Dim.1"

if (length(unique(male_scores[[group_var]])) == 2) {
  
  shapiro_results <- male_scores %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      p_value_shapiro = shapiro.test(get(value_var))$p.value,
      .groups = "drop"
    )
  
  levene_result <- leveneTest(as.formula(paste(value_var, "~", group_var)), data = male_scores)
  levene_p <- levene_result$`Pr(>F)`[1]
  
  normal_all <- all(shapiro_results$p_value_shapiro > 0.05)
  variance_equal <- (levene_p > 0.05)
  
  if (normal_all & variance_equal) {
    test_method <- "t-test"
    test_result <- t.test(as.formula(paste(value_var, "~", group_var)), data = male_scores)
  } else {
    test_method <- "Wilcoxon rank-sum test"
    test_result <- wilcox.test(as.formula(paste(value_var, "~", group_var)), data = male_scores)
  }
  
  final_results <- shapiro_results %>%
    rename(Group = !!sym(group_var),
           Shapiro_p = p_value_shapiro) %>%
    mutate(
      Levene_p = ifelse(row_number() == 1, levene_p, NA),
      Test_Method = ifelse(row_number() == 1, test_method, NA),
      Test_Statistic = ifelse(row_number() == 1, test_result$statistic, NA),
      Test_p = ifelse(row_number() == 1, test_result$p.value, NA)
    )
  
  write.csv(final_results,
            "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PC1_female_test_results.csv",
            row.names = FALSE)
}

p <- ggplot(male_scores, aes(x = zhixi, y = Dim.1, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(method = ifelse(test_method == "t-test", "t.test", "wilcox.test"),
                     label = "p.format", size = 14,
                     label.y = max(male_scores$Dim.1) * 1.05,
                     label.x = 1.5) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "Female Lizard PC1 Score Comparison by Ecomorph",
    x = "Ecomorph (zhixi)",
    y = "PC1 Score (Body Size Index)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

print(p)

ggsave("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_PC1_male_by_zhixi.pdf", 
       plot = p, width = 7, height = 9)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PCA_Combined_Scores.csv")

female_scores <- scores %>% filter(SEX == "F")

female_scores$zhixi <- as.factor(female_scores$zhixi)

test_results <- list()

shapiro_results <- female_scores %>%
  group_by(zhixi) %>%
  summarise(
    shapiro_p = shapiro.test(Dim.1)$p.value,
    .groups = "drop"
  )

levene_p <- leveneTest(Dim.1 ~ zhixi, data = female_scores)$`Pr(>F)`[1]

if (all(shapiro_results$shapiro_p > 0.05) & levene_p > 0.05) {
  test_type <- "t-test"
  stat_result <- t.test(Dim.1 ~ zhixi, data = female_scores)
} else {
  test_type <- "Wilcoxon rank-sum"
  stat_result <- wilcox.test(Dim.1 ~ zhixi, data = female_scores)
}

final_results <- data.frame(
  Group = shapiro_results$zhixi,
  Shapiro_p = round(shapiro_results$shapiro_p, 5),
  Levene_p = round(levene_p, 5),
  Test_type = test_type,
  Statistic = round(stat_result$statistic, 5),
  P_value = round(stat_result$p.value, 5)
)

write.csv(
  final_results,
  "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/PC1_female_test_results.csv",
  row.names = FALSE
)

p <- ggplot(female_scores, aes(x = zhixi, y = Dim.1, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(
    method = ifelse(test_type == "t-test", "t.test", "wilcox.test"),
    label = "p.format", size = 14,
    label.y = max(female_scores$Dim.1) * 1.05,
    label.x = 1.5
  ) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "Female Lizard PC1 Score Comparison by Ecomorph",
    x = "Ecomorph (zhixi)",
    y = "PC1 Score (Body Size Index)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

print(p)

ggsave(
  "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_PC1_female_by_zhixi.pdf",
  plot = p, width = 7, height = 9
)