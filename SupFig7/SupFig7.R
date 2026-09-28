#Loading libraries
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(scales)
library(dplyr)
library(gridExtra)

#### SupFig 7a: mCherry Filtering (for protein abundance/cell size measurements) ####

data_mcherry <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig2/pSyn_Activity.csv')
data_mcherry <- mutate_at(data_mcherry, vars(prom_nb, type, rep, pos), list(factor))

q <- ggplot(data_mcherry[data_mcherry$ORG.G.HLin > 0, ], aes(x = ORG.G.HLin)) +
  geom_histogram(bins = 100, fill = "#FE7B09", color = "black") +
  labs(x = "ORG.G", y = "Count", tag = "A") +
  geom_segment(data = data.frame(x = 4), aes(x = x, xend = x, y = 0, yend = Inf),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 1.2, color = "black") +
  scale_x_log10(labels = label_number(big.mark = "")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.1, 2000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        plot.tag = element_text(size = 24, face = "bold"),
        plot.tag.position = c(0.03, 0.97)
  )

q

path_fig = 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Figures'

#ggsave("SupFig7a.png", plot = q, path = path_fig, width = 1650, height = 1350, units = "px", dpi = 300)

#### SupFig 7b: FSC/SSC filtering (for ploidy measurements) ####

data_FSCSSC <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig5/ESP1_DNA_Content.csv')

data_FSCSSC$sample <-paste0(as.character(data_FSCSSC$plasmid),"_rep", as.character(data_FSCSSC$rep))

data_FSCSSC <- mutate_at(data_FSCSSC, vars(strain, plasmid, plate, rep, pos, sample), list(factor))

#Order samples
data_FSCSSC$plasmid <- factor(data_FSCSSC$plasmid, levels = c("stop_alpha10", "stop", "#2", "#30", "#65", "WT", "#100", "Empty", "stop_noRNA", "stop_nocol")) #Reorder levels

vertices <- data.frame(
  X = c(100, 100, 350, 600, 625, 325),
  Y = c(0, 250, 700, 700, 250, 50)
)

vertices <- rbind(vertices, vertices[1, ])

set.seed(1805)
sample_names <- unique(as.character(data_FSCSSC$sample))
sample_random <- sample(sample_names, 1)

data_random <- data_FSCSSC %>%
  dplyr::filter(sample == sample_random)

r <- ggplot(data_random, aes(x = FSC.HLin, y = SSC.HLin)) +
  geom_bin2d(bins = 250) +
  scale_fill_continuous(low="#FFCCA1",high="#FE7B09") +
  labs(x = "FSC", y = "SSC", fill = "Count", tag = "B") +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1000)) +
  geom_polygon(
    data = vertices,
    aes(x = X, y = Y),
    linewidth = 1,
    linetype = "dashed",
    color = "black",
    alpha = 0
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        plot.tag = element_text(size = 24, face = "bold"),
        plot.tag.position = c(0.03, 0.97),
        legend.position = c(0.9, 0.23)
  )

r

#ggsave("SupFig7b.png", plot = r, path = path_fig, width = 1650, height = 1350, units = "px", dpi = 300)

#### Final figure ####
s <- grid.arrange(q, r, ncol = 1)

#ggsave("SupFig7.png", plot = s, path = path_fig, width = 1650, height = 2700, units = "px", dpi = 300)
