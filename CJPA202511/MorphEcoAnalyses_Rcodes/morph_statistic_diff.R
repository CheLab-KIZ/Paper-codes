library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)
library(readr)

setwd("G:/longxi/longxibioinfo/PhenotypeStatistics/PCAnalysis")
data <- read_csv("all_data_ratio_group_imputed.csv") %>%
  mutate(
    Clade = factor(Clade, levels = c("no-Radiation", "Radiation")),
    SEX = factor(SEX, levels = c("F", "M"))
  )

traits <- c("TL", "TRL", "HL", "HW", "HH", "SEL", "FLL", "HLL", "T4L")

create_blackdot_plot <- function(trait, sex) {
  plot_data <- data %>% filter(SEX == sex)
  
  ggplot(plot_data, aes(x = Clade, y = !!sym(trait), fill = Clade)) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.7,
      outlier.shape = NA,
      color = "black",
      size = 1
    ) +
    geom_jitter(
      width = 0.15,
      size = 1.2,
      alpha = 0.6,
      shape = 16,
      color = "black"
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 5,
      label.y = max(plot_data[[trait]], na.rm = TRUE) * 1.1,
      label.x = 1.5
    ) +
    scale_fill_manual(values = c("no-Radiation" = "#24A1C2", "Radiation" = "#BD213B")) +
    labs(
      x = NULL,
      y = trait,
      title = paste(ifelse(sex == "F", "Female", "Male"), trait)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

generate_final_plot <- function(sex) {
  plots <- map(traits, ~ create_blackdot_plot(.x, sex))
  wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = paste(ifelse(sex == "F", "Female", "Male"), "Morphological Traits Comparison (Raw Data)"),
      theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold", margin = margin(b = 20))
      )
    )
}

female_plot <- generate_final_plot("F")
male_plot <- generate_final_plot("M")

print(female_plot)
print(male_plot)

wilcox_results <- map_dfr(traits, function(trait) {
  bind_rows(
    data %>% 
      filter(SEX == "F") %>%
      wilcox_test(as.formula(paste(trait, "~ Clade"))) %>%
      mutate(sex = "F", trait = trait, .before = 1),
    data %>% 
      filter(SEX == "M") %>%
      wilcox_test(as.formula(paste(trait, "~ Clade"))) %>%
      mutate(sex = "M", trait = trait, .before = 1)
  ) %>%
    select(sex, trait, group1, group2, n1, n2, statistic, p, method, alternative) %>%
    rename(Z = statistic)
})

output_dir <- "Wilcoxon_BlackDots_Results"
dir.create(output_dir, showWarnings = FALSE)

write_csv(wilcox_results, file.path(output_dir, "wilcoxon_test_results_withZ.csv"))

ggsave(
  file.path(output_dir, "female_blackdots.pdf"),
  female_plot,
  width = 14,
  height = 11,
  dpi = 300
)

ggsave(
  file.path(output_dir, "male_blackdots.pdf"),
  male_plot,
  width = 14,
  height = 11,
  dpi = 300
)