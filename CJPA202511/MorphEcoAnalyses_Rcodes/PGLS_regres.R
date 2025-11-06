library(ape)
library(caper)
library(nlme)
library(ggplot2)
library(geiger)

tree <- read.tree("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/diplo_timetree.tre")
data <- read.csv("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/SVL_PCA_ENV.csv")

name.check(tree, data, data.names = data$species)
data <- data[data$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
rownames(data) <- data$species

comp_data <- comparative.data(tree, data, names.col = "species", vcv = TRUE, warn.dropped = TRUE)

model_lambda <- pgls(Female_PC1 ~ bio_15, data = comp_data, lambda = "ML")
lambda_est <- model_lambda$param["lambda"]
cat("Estimated lambda =", lambda_est, "\n")

cor_pagel <- corPagel(value = lambda_est, phy = tree, fixed = TRUE)
model_gls <- gls(Female_PC1 ~ bio_15, data = data, correlation = cor_pagel, method = "ML")

summary_model <- summary(model_gls)$tTable
estimate <- round(summary_model[2, 1], 4)
p_value <- signif(summary_model[2, 4], 4)

p <- ggplot(data, aes(x = bio_15, y = Female_PC1)) +
  geom_point(size = 4, color = "#000000", alpha = 0.8) +
  stat_smooth(method = "lm", se = TRUE,
              color = "#ED6547",
              fill = "#CDDBEB",
              alpha = 0.3,
              linetype = "solid",
              linewidth = 1.5) +
  labs(
    title = "PGLS Regression: Female_PC1 ~ bio_15",
    subtitle = paste("λ =", round(lambda_est, 3), ", β =", estimate, ", p =", p_value),
    x = "bio_15",
    y = "Female_PC1"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black")
  )

print(p)

library(ape)
library(caper)
library(nlme)
library(ggplot2)
library(geiger)
library(patchwork)

tree <- read.tree("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/diplo_timetree.tre")
data <- read.csv("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/SVL_PCA_ENV.csv")

data <- data[data$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
rownames(data) <- data$species

comp_data <- comparative.data(tree, data, names.col = "species", vcv = TRUE, warn.dropped = TRUE)

pgls_plot <- function(var_name) {
  model_lambda <- pgls(as.formula(paste("Female_PC1 ~", var_name)),
                       data = comp_data, lambda = "ML")
  lambda_est <- model_lambda$param["lambda"]
  
  cor_pagel <- corPagel(value = lambda_est, phy = tree, fixed = TRUE)
  model_gls <- gls(as.formula(paste("Female_PC1 ~", var_name)),
                   data = data, correlation = cor_pagel, method = "ML")
  
  summary_model <- summary(model_gls)$tTable
  estimate <- round(summary_model[2, 1], 4)
  p_value <- signif(summary_model[2, 4], 4)
  
  sig <- ifelse(p_value < 0.001, "***",
                ifelse(p_value < 0.01, "**",
                       ifelse(p_value < 0.05, "*", "")))
  
  ggplot(data, aes(x = .data[[var_name]], y = Female_PC1)) +
    geom_point(size = 4, color = "#000000", alpha = 0.8) +
    stat_smooth(method = "lm", se = TRUE,
                color = "#ED6547",
                fill = "#CDDBEB",
                alpha = 0.3,
                linetype = "solid",
                linewidth = 1.5) +
    labs(
      title = paste("PGLS Regression: Female_PC1 ~", var_name),
      subtitle = paste("λ =", round(lambda_est, 3),
                       ", β =", estimate,
                       ", p =", p_value, sig),
      x = var_name,
      y = "Female_PC1"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(color = "black")
    )
}

vars <- c("Elevation", "bio_1", "bio_4", "bio_12", "bio_15", "NPP")

plots <- lapply(vars, pgls_plot)

final_plot <- wrap_plots(plots, ncol = 2)

print(final_plot)

library(caper)
library(nlme)
library(ggplot2)
library(geiger)
library(patchwork)

tree <- read.tree("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/diplo_timetree.tre")
data <- read.csv("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/SVL_PCA_ENV.csv")

data <- data[data$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
rownames(data) <- data$species

comp_data <- comparative.data(tree, data, names.col = "species", vcv = TRUE, warn.dropped = TRUE)

pgls_plot_and_summary <- function(var_name) {
  model_lambda <- pgls(as.formula(paste("Female_PC1 ~", var_name)),
                       data = comp_data, lambda = "ML")
  lambda_est <- model_lambda$param["lambda"]
  
  cor_pagel <- corPagel(value = lambda_est, phy = tree, fixed = TRUE)
  model_gls <- gls(as.formula(paste("Female_PC1 ~", var_name)),
                   data = data, correlation = cor_pagel, method = "ML")
  
  summary_model <- summary(model_gls)$tTable
  estimate <- round(summary_model[2, 1], 4)
  se <- round(summary_model[2, 2], 4)
  t_value <- round(summary_model[2, 3], 4)
  p_value <- signif(summary_model[2, 4], 4)
  
  rss <- sum(residuals(model_gls)^2)
  tss <- sum((data$Female_PC1 - mean(data$Female_PC1))^2)
  R_squared <- round(1 - rss/tss, 4)
  
  result <- data.frame(
    Variable = var_name,
    Lambda = round(lambda_est, 4),
    Beta = estimate,
    SE = se,
    t_value = t_value,
    p_value = p_value,
    R_squared = R_squared
  )
  
  p <- ggplot(data, aes(x = .data[[var_name]], y = Female_PC1)) +
    geom_point(size = 4, color = "#000000", alpha = 0.8) +
    stat_smooth(method = "lm", se = TRUE,
                color = "#ED6547",
                fill = "#CDDBEB",
                alpha = 0.3,
                linetype = "solid",
                linewidth = 1.5) +
    labs(
      title = paste("PGLS Regression: Female_PC1 ~", var_name),
      subtitle = paste("λ =", round(lambda_est, 3),
                       ", β =", estimate,
                       ", p =", p_value),
      x = var_name,
      y = "Female_PC1"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(color = "black")
    )
  
  return(list(plot = p, result = result))
}

vars <- c("Elevation", "bio_1", "bio_4", "bio_12", "bio_15", "NPP")

all_results <- lapply(vars, pgls_plot_and_summary)

plots <- lapply(all_results, `[[`, "plot")

final_plot <- wrap_plots(plots, ncol = 2)
print(final_plot)

summary_table <- do.call(rbind, lapply(all_results, `[[`, "result"))
print(summary_table)

write.csv(summary_table, "G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/PGLSRegressionGraphs/Female_PGLS_summary.csv", row.names = FALSE)

library(caper)
library(nlme)
library(ggplot2)
library(geiger)
library(patchwork)

tree <- read.tree("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/diplo_timetree.tre")
data <- read.csv("G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/SVL_PCA_ENV.csv")

data <- data[data$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
rownames(data) <- data$species

comp_data <- comparative.data(tree, data, names.col = "species", vcv = TRUE, warn.dropped = TRUE)

pgls_plot_and_summary <- function(var_name) {
  model_lambda <- pgls(as.formula(paste("Male_PC1 ~", var_name)),
                       data = comp_data, lambda = "ML")
  lambda_est <- model_lambda$param["lambda"]
  
  cor_pagel <- corPagel(value = lambda_est, phy = tree, fixed = TRUE)
  model_gls <- gls(as.formula(paste("Male_PC1 ~", var_name)),
                   data = data, correlation = cor_pagel, method = "ML")
  
  summary_model <- summary(model_gls)$tTable
  estimate <- round(summary_model[2, 1], 4)
  se <- round(summary_model[2, 2], 4)
  t_value <- round(summary_model[2, 3], 4)
  p_value <- signif(summary_model[2, 4], 4)
  
  rss <- sum(residuals(model_gls)^2)
  tss <- sum((data$Male_PC1 - mean(data$Male_PC1))^2)
  R_squared <- round(1 - rss/tss, 4)
  
  result <- data.frame(
    Variable = var_name,
    Lambda = round(lambda_est, 4),
    Beta = estimate,
    SE = se,
    t_value = t_value,
    p_value = p_value,
    R_squared = R_squared
  )
  
  p <- ggplot(data, aes(x = .data[[var_name]], y = Male_PC1)) +
    geom_point(size = 4, color = "#000000", alpha = 0.8) +
    stat_smooth(method = "lm", se = TRUE,
                color = "#ED6547",
                fill = "#CDDBEB",
                alpha = 0.3,
                linetype = "solid",
                linewidth = 1.5) +
    labs(
      title = paste("PGLS Regression: Male_PC1 ~", var_name),
      subtitle = paste("λ =", round(lambda_est, 3),
                       ", β =", estimate,
                       ", p =", p_value),
      x = var_name,
      y = "Male_PC1"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(color = "black")
    )
  
  return(list(plot = p, result = result))
}

vars <- c("Elevation", "bio_1", "bio_4", "bio_12", "bio_15", "NPP")

all_results <- lapply(vars, pgls_plot_and_summary)

plots <- lapply(all_results, `[[`, "plot")

final_plot <- wrap_plots(plots, ncol = 2)
print(final_plot)

summary_table <- do.call(rbind, lapply(all_results, `[[`, "result"))
print(summary_table)

write.csv(summary_table, "G:/longxi/longxibioinfo/PhenoEnvCorr/PLGSAnalysis/PGLSRegressionGraphs/Male_PGLS_summary.csv", row.names = FALSE)