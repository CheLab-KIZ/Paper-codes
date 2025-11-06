library(dplyr)
library(ggplot2)
library(ggpubr)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/all_data_imputed.csv")

male_scores <- scores %>% filter(SEX == "F")

if (length(unique(male_scores$zhixi)) == 2) {
  t_test_male <- t.test(SVL ~ zhixi, data = male_scores)
  print(t_test_male)
  
  capture.output(t_test_male, file = "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/SVL_t_test_male_by_zhixi.txt")
}

p <- ggplot(male_scores, aes(x = zhixi, y = SVL, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(method = "t.test", label = "p.format", size = 14,
                     label.y = max(male_scores$SVL, na.rm = TRUE) * 1.05,
                     label.x = 1.5) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "SVL Comparison of Male Lizards by Ecotype",
    x = "Ecotype (zhixi)",
    y = "Snout-Vent Length SVL (mm)"
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

ggsave("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_SVL_male_by_zhixi.pdf", 
       plot = p, width = 7, height = 9)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)
library(broom)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/all_data_imputed.csv")

male_scores <- scores %>% filter(SEX == "M")

male_scores$zhixi <- factor(male_scores$zhixi)

shapiro_results <- male_scores %>%
  group_by(zhixi) %>%
  summarise(
    shapiro_p = shapiro.test(SVL)$p.value,
    .groups = "drop"
  )

levene_res <- leveneTest(SVL ~ zhixi, data = male_scores)
levene_p <- levene_res[1, "Pr(>F)"]

normal_all <- all(shapiro_results$shapiro_p > 0.05)
equal_var <- levene_p > 0.05

if (normal_all) {
  if (equal_var) {
    test_res <- t.test(SVL ~ zhixi, data = male_scores, var.equal = TRUE)
    test_method <- "t-test (equal var)"
  } else {
    test_res <- t.test(SVL ~ zhixi, data = male_scores, var.equal = FALSE)
    test_method <- "t-test (unequal var)"
  }
} else {
  test_res <- wilcox.test(SVL ~ zhixi, data = male_scores)
  test_method <- "Wilcoxon rank-sum test"
}

final_results <- data.frame(
  zhixi_group = paste(shapiro_results$zhixi, collapse = " vs "),
  Shapiro_p_group1 = shapiro_results$shapiro_p[1],
  Shapiro_p_group2 = shapiro_results$shapiro_p[2],
  Levene_p = levene_p,
  Test_method = test_method,
  Test_statistic = unname(test_res$statistic),
  Test_p = test_res$p.value
)

write.csv(final_results,
          "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/SVL_Male_Test_Results.csv",
          row.names = FALSE)

p <- ggplot(male_scores, aes(x = zhixi, y = SVL, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(method = ifelse(normal_all,
                                     ifelse(equal_var, "t.test", "t.test"),
                                     "wilcox.test"),
                     label = "p.format", size = 14,
                     label.y = max(male_scores$SVL, na.rm = TRUE) * 1.05,
                     label.x = 1.5) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2",
                               "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "SVL Comparison of Male Lizards by Ecotype",
    x = "Ecotype (zhixi)",
    y = "Snout-Vent Length SVL (mm)"
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

ggsave("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_SVL_male_by_zhixi.pdf",
       plot = p, width = 7, height = 9)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)

scores <- read.csv("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/all_data_imputed.csv")

male_scores <- scores %>% filter(SEX == "M")

groups <- unique(male_scores$zhixi)
if (length(groups) != 2) {
  stop("zhixi groups are not two, cannot perform two-group comparison!")
}

shapiro_results <- male_scores %>%
  group_by(zhixi) %>%
  summarise(
    shapiro_p = shapiro.test(SVL)$p.value,
    .groups = "drop"
  )

levene_result <- leveneTest(SVL ~ zhixi, data = male_scores)

if (all(shapiro_results$shapiro_p > 0.05)) {
  if (levene_result$`Pr(>F)`[1] > 0.05) {
    test_used <- "t-test (equal variance)"
    stat_result <- t.test(SVL ~ zhixi, data = male_scores, var.equal = TRUE)
  } else {
    test_used <- "t-test (unequal variance)"
    stat_result <- t.test(SVL ~ zhixi, data = male_scores, var.equal = FALSE)
  }
} else {
  test_used <- "Mann-Whitney U test"
  stat_result <- wilcox.test(SVL ~ zhixi, data = male_scores)
}

result_table <- data.frame(
  Group = shapiro_results$zhixi,
  Shapiro_p = shapiro_results$shapiro_p
)
result_table <- rbind(
  result_table,
  data.frame(Group = "Levene's Test", Shapiro_p = levene_result$`Pr(>F)`[1])
)
result_table <- rbind(
  result_table,
  data.frame(Group = test_used, Shapiro_p = stat_result$p.value)
)

write.csv(result_table, "G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/SVL_test_results_male.csv", row.names = FALSE)

p <- ggplot(male_scores, aes(x = zhixi, y = SVL, fill = zhixi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  stat_compare_means(method = ifelse(test_used == "Mann-Whitney U test", "wilcox.test", "t.test"),
                     label = "p.format", size = 14,
                     label.y = max(male_scores$SVL, na.rm = TRUE) * 1.05,
                     label.x = 1.5) +
  scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
  theme_minimal() +
  labs(
    title = "SVL Comparison of Male Lizards by Ecotype",
    x = "Ecotype (zhixi)",
    y = "Snout-Vent Length SVL (mm)"
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

ggsave("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAanalysis/PCA_results6/Boxplot_SVL_male_by_zhixi.pdf", 
       plot = p, width = 7, height = 9)