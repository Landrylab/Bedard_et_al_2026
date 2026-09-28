####Load libraries####
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(RColorBrewer)

setwd('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/SupFig8')

plate_plan_file<-"./Plan_Growth_Curves.xlsx"
plate_file<-"./Growth_Curves_Conditions.xlsx"
output_file<-"./stats_growth_curve_2025_11_12.csv"

#Format df info
info<-readxl::read_xlsx(path=plate_plan_file)%>%
  mutate(replicate=as.factor(replicate),
         drug_conc=as.factor(drug_conc),
         media=as.factor(media),
         drug=as.factor(drug))%>%
  mutate(well=paste0(row, column))

# Format absorbance measurements
growth_data<-readxl::read_xlsx(plate_file, skip=24)%>%
  .[,-1]%>%
  slice(1:179) %>%
  filter(!(is.na(as.numeric(Time))))%>%
  mutate(Time=as.numeric(Time)*24,
         plate=1)%>%
  dplyr::select(-"T° 600")%>%
  mutate(across(everything(), as.numeric))%>%
  pivot_longer(cols=-c(Time, plate), names_to = "well", values_to = "abs")%>%
  left_join(., info, by=c("plate", "well"))

#Separate data
growth_data_media <- growth_data %>% filter(drug == "None")
growth_data_Benomyl <- growth_data[growth_data$drug %in% c("Benomyl", "DMSO"), ]
growth_data_MMS <- growth_data[growth_data$drug %in% c("MMS", "DMSO"), ]

# SupFig 8a
my_colors <- brewer.pal(9, "YlOrBr")[3:9]

q <- growth_data_Benomyl %>%
  ggplot(aes(x = Time,
             y = abs,
             colour = factor(drug_conc))) +
  
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    aes(fill = factor(drug_conc)),
    alpha = 0.08,
    size = 1.2
  ) +
  labs(x = "Time (h)", y = "Optical density", tag="A") +
  scale_colour_manual(
    values = my_colors,
    name = "Concentration (μg/mL)"
  ) +
  scale_fill_manual(
    values = my_colors,
    name = "Concentration (μg/mL)"
  ) +
  
  ggtitle("Benomyl") +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    legend.position = c(0.05, 0.95),
    legend.title = element_text(size=12),
    legend.text = element_text(size =12),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.03, 0.97)
  )
q

#### SupFig 8b ####
my_colors <- brewer.pal(9, "YlOrBr")[3:7]


r <- growth_data_MMS %>%
  ggplot(aes(x = Time,
             y = abs,
             colour = factor(drug_conc))) +
  
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    aes(fill = factor(drug_conc)),
    alpha = 0.08,
    size = 1.2
  ) +
  labs(x = "Time (h)", y = "Optical density", tag= "B") +
  scale_colour_manual(
    values = my_colors,
    name = "Concentration (μg/mL)"
  ) +
  scale_fill_manual(
    values = my_colors,
    name = "Concentration (μg/mL)"
  ) +
  ggtitle("MMS") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    legend.position = c(0.05, 0.95),
    legend.title = element_text(size=12),
    legend.text = element_text(size =12),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.03, 0.97)
  )
r

# Grid
s <- grid.arrange(q, r, ncol = 1)

# Save
path_fig = 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Figures'
#ggsave("SupFig8.png", plot = s, path = path_fig, width = 2100, height = 2700, units = "px", dpi = 300)
