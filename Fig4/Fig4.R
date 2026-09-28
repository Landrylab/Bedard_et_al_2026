# Loading libraries
rm(list = ls())
library(ggplot2)
library(dplyr)

# Data preparation
data <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig4/Expression_Measurements.csv')
data$GRN.FSC <-data$GRN.B.HLin/data$FSC.HLin
data$log_GRN.FSC <-log10(data$GRN.FSC)

data <- data[data$ORG.G.HLin > 4, ]

df_summary <- data %>%
  group_by(gene, type, rep) %>%
  summarise(
    median_GRN.FSC = median(GRN.FSC, na.rm = TRUE),
    median_log_GRN.FSC = median(log_GRN.FSC, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_filtered <- df_summary %>% dplyr::filter(gene != "Negative")
df_summary_filtered <- mutate_at(df_summary_filtered, vars(gene, type, rep), list(factor))

df_summary_filtered$type <- factor(df_summary_filtered$type, levels = c("single","dup", "99"))
levels(df_summary_filtered$type) <- c("single","duplication", "pSyn99")

df_summary_filtered <- df_summary_filtered %>%
  group_by(gene) %>%
  mutate(
    median_fluo = median(median_GRN.FSC[type == "single"], na.rm = TRUE),
    normalized_fluo = median_GRN.FSC / median_fluo
  ) %>%
  ungroup()

# Fig 4
Fig4 <- ggplot(df_summary_filtered, aes(x = gene, y = median_GRN.FSC, fill = type)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  labs(x = "Gene", y = "Fluorescence (AU)", fill = "Strain") +
  scale_fill_manual(values = c(
    "single" = "#00BA38",
    "duplication" = "#56B1F7",
    "pSyn99" = "#FE7B09"
  )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, color = "black", face="italic"),
        axis.text.y = element_text(size = 20, color = "black"),
        legend.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size =20, margin = margin(r = 113)),
        legend.position = c(0.55, 0.9),
        legend.direction = "horizontal")
Fig4

# Save
path_fig = 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Figures'
#ggsave("Fig4.png", plot = Fig4, path = path_fig, width = 4000, height = 1800, units = "px", dpi = 300)
#ggsave("Fig4.svg", plot = Fig4, path = path_fig, width = 4000, height = 1800, units = "px", dpi = 300)
