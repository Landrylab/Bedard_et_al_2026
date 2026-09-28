# Loading libraries
rm(list = ls())
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(patchwork)

#### Fig 2a ####

# Data preparation
data_Fig2a <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/pSyn_Activity.csv')
data_Fig2a <- data_Fig2a[data_Fig2a$ORG.G.HLin > 4, ]
data_Fig2a$ratio.GRNFSC <- data_Fig2a$GRN.B.HLin/data_Fig2a$FSC.HLin
data_Fig2a$logratio.GRNFSC <- log10(data_Fig2a$ratio.GRNFSC)

data_exp <- data_Fig2a %>% filter(prom_nb != "ctrl")
data_exp$prom_nb <- as.numeric(as.character(data_exp$prom_nb)) 

# Summarize by replicate
df_data_exp <- data_exp %>%
  group_by(prom_nb, rep) %>%
  summarise(
    median_log_ratio_GRNFSC = median(logratio.GRNFSC, na.rm = TRUE),
    .groups = "drop"
  )

#Add promoter estimates to dataframe
prom_estimates <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/DeBoer_estimates.csv', sep = ";")

df_promoters <- merge(df_data_exp, prom_estimates, by = "prom_nb", all = TRUE)
df_promoters$prom_name <- paste0("pSyn", df_promoters$prom_nb)

# Spearman correlation
df_spearman <- df_promoters %>%
  group_by(prom_nb, expr_deboer) %>%
  summarise(
    medianGRNFSC_per_prom = median(median_log_ratio_GRNFSC, na.rm = TRUE),
    .groups = "drop"
  )

spearman_uncorr <- cor.test(df_spearman$medianGRNFSC_per_prom, df_spearman$expr_deboer, method = "spearman")
print(spearman_uncorr)

# Figure preparation

# dataframe preparation for graph
df_summary <- df_promoters %>%
  group_by(expr_deboer, prom_name) %>%
  summarise(
    mean_median_log_ratio_GRNFSC = mean(median_log_ratio_GRNFSC),
    sd_median_log_ratio_GRNFSC = sd(median_log_ratio_GRNFSC),
    n = n(),
    se_median_log_ratio_GRNFSC = sd_median_log_ratio_GRNFSC / sqrt(n))

# Adjust promoter name labels positions for graphic
fit <- lm(median_log_ratio_GRNFSC ~ expr_deboer, data = df_promoters)


y_pred <- predict(fit, newdata = df_summary)

df_labels <- df_summary %>%
  ungroup() %>%
  mutate(
    y_line = y_pred,
    position = ifelse(mean_median_log_ratio_GRNFSC > y_line, "above", "below"),
    nudge_x = ifelse(position == "above", -0.5, 0.15),
    nudge_y = ifelse(position == "above",  0.075, -0.05),
    hjust   = ifelse(position == "above", 0.5, 0),
    vjust   = ifelse(position == "above", 0.5, 0)
  )

#Load plasmid icon
plasmid_png_path <- 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/Syn_prom_icon.png'

img <- readPNG(plasmid_png_path)
grob <- rasterGrob(img, interpolate = TRUE)

# Fig 2a
Fig2a <- ggplot() + 
  geom_point(data = df_summary, aes(x = expr_deboer, y = mean_median_log_ratio_GRNFSC), size = 3) +
  geom_errorbar(data = df_summary, aes(x = expr_deboer, ymin = mean_median_log_ratio_GRNFSC - sd_median_log_ratio_GRNFSC, ymax = mean_median_log_ratio_GRNFSC + sd_median_log_ratio_GRNFSC), width = 0.1) +
  stat_smooth(data = df_promoters, aes(x = expr_deboer, y = median_log_ratio_GRNFSC), formula = y ~ x, method = "lm", col = "darkgreen", se = TRUE) +
  geom_text(
    data = df_labels,
    aes(
      x = expr_deboer + nudge_x,
      y = mean_median_log_ratio_GRNFSC + nudge_y,
      label = prom_name,
      hjust = hjust,
      vjust = vjust
    ),
    size = 4.5
  ) +
  scale_x_continuous(breaks = seq(0, 16, by = 4)) +
  scale_y_continuous(breaks = seq(-3, 0, by = 0.5)) +
  labs(
    x = expression("Expression (DeBoer " * italic("et al.") * " 2020))"),
    y = expression(Log[10] ~ "(Fluorescence (AU))")
  ) +
  annotation_custom(grob, xmin = 1.25, xmax = 9.25, ymin = -2.25, ymax = -1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
  )
Fig2a

# Save figure
path_fig = 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Figures'
#ggsave("Fig2a.png", plot = Fig2a, path = path_fig, width = 2400, height = 1350, units = "px", dpi = 300)

#### Fig 2b ####

# Data preparation
data_Fig2b <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/NativeProm_Activity.csv')
data_Fig2b <- mutate_at(data_Fig2b, vars(Promoter, Type), list(factor))
data_Fig2b <- data_Fig2b[data_Fig2b$ORG.G.HLin > 4, ]
data_Fig2b$ratio_GRN.FSC <-data_Fig2b$GRN.B.HLin/data_Fig2b$FSC.HLin


# Summarize by replicate
data_summary_2 <- data_Fig2b %>%
  group_by(Promoter, Type, rep) %>%
  summarise(
    median_GRN = median(GRN.B.HLin, na.rm = TRUE),
    median_ratio_GRN.FSC = median(ratio_GRN.FSC, na.rm = TRUE),
    log_median_GRNFSC = log10(median_ratio_GRN.FSC),
    .groups = "drop")

#Rename some levels
data_summary_2$Promoter <- factor(
  dplyr::recode(as.character(data_summary_2$Promoter),
                "pKB134" = "Empty",
                "2" = "pSyn2",
                "100" = "pSyn100")
)

# Order the promoter levels
data_summary_2$Promoter <- factor(data_summary_2$Promoter,
                                  levels = c("Empty", "pSyn2", "pSyn100", "HIP1", "GCR1", "CDC25", "WBP1", "CFT1",
                                             "SPC98", "PDC2", "TTI1", "SYF1", "ESP1",
                                             "CCA1"))

#calculate mean for each promoter
means_data <- data_summary_2 %>%
  group_by(Promoter) %>%
  summarise(mean_log_GRNFSC = mean(log_median_GRNFSC, na.rm = TRUE)) %>%
  ungroup()

median_2 <- median(data_summary_2$log_median_GRNFSC[data_summary_2$Promoter == "pSyn2"])
median_100 <- median(data_summary_2$log_median_GRNFSC[data_summary_2$Promoter == "pSyn100"])

#Load plasmid icon
plasmid_png_path_2 <- 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/WT_prom_icon.png'

img_2 <- readPNG(plasmid_png_path_2)
grob_2 <- rasterGrob(img_2, interpolate = TRUE)

# Fig 2b
Fig2b <- ggplot(data_summary_2, aes(x = Promoter, y = log_median_GRNFSC)) +
  geom_boxplot(width = 0.5, fill="gray", outlier.shape = NA, size = 0.7) +
  scale_y_continuous(breaks = seq(-3, -1.5, by = 0.5)) +
  geom_hline(yintercept=median_2, linetype='dashed', col = 'black')+
  geom_hline(yintercept=median_100, linetype='dashed', col = 'black')+
  labs(x = "Promoter", y = expression(Log[10] ~ "(Fluorescence (AU))")) +
  annotation_custom(grob_2, xmin = 1, xmax = 6, ymin = -2.5, ymax = -2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1,size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size =16),
        legend.position = c(0.9, 0.75))+
  geom_jitter()
Fig2b

#ggsave("Fig2b.png", plot = Fig2b, path = path_fig, width = 2400, height = 1350, units = "px", dpi = 300)

#### Final figure ####

final_plot <- (Fig2a / Fig2b)
final_plot

#ggsave("Figure2_V6.png", plot = final_plot, path = path_fig, width = 2400, height = 2700, units = "px", dpi = 300)
#ggsave("Figure2_V6.svg", plot = final_plot, path = path_fig, width = 2400, height = 2700, units = "px", dpi = 300)
