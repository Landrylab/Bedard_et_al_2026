# Loading libraries
rm(list = ls())
library(ggplot2)
library(dplyr)
library(stringr)
library(car)
library(emmeans)
library(mixtools)
library(tidyverse)

#### Fig 5b ####

data_5b <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig5/ESP1_Expression.csv')
data_5b$GRN.FSC <-data_5b$GRN.B.HLin/data_5b$FSC.HLin
data_5b <- mutate_at(data_5b, vars(plasmid, rep, plate, pos), list(factor))

data_5b <- data_5b[data_5b$ORG.G.HLin > 4, ]

# Removing the empty plasmid (a control) level
data_5b_exp <- data_5b %>% filter(plasmid != "Empty")

df_summary_5b <- data_5b_exp %>%
  group_by(plasmid, rep) %>%
  summarise(
    median_GRN.FSC = median(GRN.FSC, na.rm = TRUE),
    .groups = "drop"
  )

#Remove weird replicate
df_summary_5b <- df_summary_5b %>%
  filter(!(plasmid == "stop_notag" & rep == "7"))


df_summary_5b$plasmid <- factor(df_summary_5b$plasmid, levels = c("stop_notag","stop", "#2", "#30", "#65", "WT", "#100"))
levels(df_summary_5b$plasmid) <- c("neg.", "ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100")

#Add State
df_summary_5b <- df_summary_5b %>%
  mutate(
    State = case_when(
      plasmid == "neg." ~ "Negative",
      plasmid == "ctrl" ~ "Single",
      plasmid == "WT" ~ "Duplicated",
      str_starts(as.character(plasmid), "pSyn") ~ "Synthetic",
      )
  )
`
#t-test`
df_ttest <- df_summary_5b %>%
  filter(plasmid %in% c("pSyn65", "WT"))

t.test(median_GRN.FSC ~ plasmid, data = df_ttest, var.equal = TRUE)

#conditions
by(df_ttest$median_GRN.FSC, df_ttest$plasmid, shapiro.test)
leveneTest(median_GRN.FSC ~ plasmid, data = df_ttest)

# Fig5b
Fig5b <- ggplot(df_summary_5b, aes(x = plasmid, y = median_GRN.FSC, fill = State)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, size = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "Negative" = "grey",
    "Single" = "#00BA38",
    "Duplicated" = "#56B1F7",
    "Synthetic" = "#FE7B09"
  ),
  breaks = c("Negative", "Single", "Duplicated", "Synthetic")) +
  labs(
    x = "Promoter of plasmid gene copy",
    y = "Fluorescence (AU)") +
  geom_segment(aes(x = 5, xend = 6, y = 0.00207, yend = 0.00207), inherit.aes = FALSE, linewidth = 0.7) +
  geom_text(aes(x = 5.5, y = 0.002085, label = "*"), inherit.aes = FALSE, size = 9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        axis.text.x = element_text(size = 16.5, color = "black"),
        axis.text.y = element_text(size = 16.5, color = "black"),
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        legend.position = c(0.15, 0.75))
Fig5b

path_fig = 'C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Figures'
#ggsave("Figure5b.png", plot = Fig5b, path = path_fig, width = 2250, height = 1350, units = "px", dpi = 300)

#### Fig 5d ####

data_5d <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig5/ESP1_FSC.csv')
data_5d <- mutate_at(data_5d, vars(strain, plasmid, plate, rep, pos), list(factor))

# Data preparation
data_5d <- data_5d[data_5d$ORG.G.HLin > 4, ] 

# Removing the empty plasmid (control) level
data_5d_exp <- data_5d %>% filter(plasmid != "Empty")

# Summarize
df_summary_5d <- data_5d_exp %>%
  group_by(plasmid, rep) %>%
  summarise(
    median_FSC = median(FSC.HLin, na.rm = TRUE),
    median_GRN = median(GRN.B.HLin, na.rm = TRUE),
    median_ORG = median(ORG.G.HLin, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_5d <- mutate_at(df_summary_5d, vars(plasmid, rep), list(factor))

df_summary_5d$plasmid <- factor(df_summary_5d$plasmid, levels= c("stop", "#2", "#30", "#65", "WT", "#100"))
levels(df_summary_5d$plasmid) <- c("ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100")

df_summary_5d <- df_summary_5d %>%
  mutate(
    State = case_when(
      plasmid == "ctrl" ~ "Single",
      plasmid == "WT" ~ "Duplicated",
      str_starts(as.character(plasmid), "pSyn") ~ "Synthetic",
    )
  )

# ANOVA
model <- lm(median_FSC ~ plasmid, data = df_summary_5d)
anova(model)

# Conditions
leveneTest(median_FSC ~ plasmid, data = df_summary_5d)

layout(matrix(c(1:5, 5), 2, 3))
plot(model, c(1, 3:5), pch = 20, cex = 1.5)
qqPlot(resid(model), "norm", pch = 20, cex = 1.5, line = "quartiles")
layout(1)

# Tukey test
TCM <- emmeans(model, ~ plasmid)
cld_result <- multcomp::cld(TCM, Letters = LETTERS) 

#Pairwise comparisons
tukey_results <- pairs(TCM, adjust = "none")
summary(tukey_results) |>
  subset(contrast == "ctrl - pSyn100")

letters_df <- data.frame(
  plasmid = c("ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100"),
  label = c("a", "a", "a", "b", "c", "d"),
  y = 845 
)



Fig5d <- ggplot(df_summary_5d, aes(x = plasmid, y = median_FSC, fill = State)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, size = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "Single" = "#00BA38",
    "Duplicated" = "#56B1F7",
    "Synthetic" = "#FE7B09"
  ),
  breaks = c("Single", "Duplicated", "Synthetic")) +
  scale_y_continuous(limits = c(650, 850)) +
  labs(x = "Promoter of plasmid gene copy", y = "FSC (AU)") +
  geom_text(data = letters_df, aes(x = plasmid, y = y, label = label),
            inherit.aes = FALSE, size = 7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        axis.text.x = element_text(size = 16.5, color = "black"),
        axis.text.y = element_text(size = 16.5, color = "black"),
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        legend.position = c(0.85, 0.75))
Fig5d

#ggsave("Figure5d.png", plot = Fig5d, path = path_fig, width = 2250, height = 1350, units = "px", dpi = 300)

#### Fig 5f ####

data_5f <- read.csv('C:/Users/fredb/OneDrive - Université Laval/LandryLab/Projet/Duplication_dosage_paper/Data/Fig5/ESP1_DNA_Content.csv')

# Data preparation
data_5f$sample <-paste0(as.character(data_5f$plasmid),"_rep", as.character(data_5f$rep))
data_5f <- mutate_at(data_5f, vars(strain, plasmid, plate, rep, pos, sample), list(factor))

data_5f$plasmid <- factor(data_5f$plasmid, levels = c("stop_alpha10", "stop", "#2", "#30", "#65", "WT", "#100", "Empty", "stop_noRNA", "stop_nocol")) #Reorder levels

# FSC & SSC filtering
vertices <- data.frame(
  X = c(100, 100, 350, 600, 625, 325),
  Y = c(0, 250, 700, 700, 250, 50)
)

vertices <- rbind(vertices, vertices[1, ])


data_5f$inside_polygon <- sp::point.in.polygon(
  point.x = data_5f$FSC.HLin,
  point.y = data_5f$SSC.HLin,
  pol.x = vertices$X,
  pol.y = vertices$Y
)

data_5f_gated <- data_5f[data_5f$inside_polygon > 0, ]

#Remove the controls for analysis
data_5f_exp <- data_5f_gated %>%
  filter(!plasmid %in% c("stop_alpha10", "Empty", "stop_noRNA", "stop_nocol"))

data_5f_exp <- mutate_at(data_5f_exp, vars(plasmid), list(factor))

data_5f_exp$plasmid <- factor(data_5f_exp$plasmid, levels = c("stop", "#2", "#30", "#65", "WT", "#100"))
levels(data_5f_exp$plasmid) <- c("ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100")

data_5f_exp$sample <-paste0(as.character(data_5f_exp$plasmid),"_rep", as.character(data_5f_exp$rep))
data_5f_exp <- mutate_at(data_5f_exp, vars(sample), list(factor))

# Gaussian mixture model
reorder_by_mu <- function(fit) {
  ord <- order(fit$mu) 
  
  fit$mu     <- fit$mu[ord]
  fit$sigma  <- fit$sigma[ord]
  fit$lambda <- fit$lambda[ord]
  
  fit
}

set.seed(3)


k_val <- 3
mu_init <- c(85, 150, 250)

fits <- data_5f_exp %>%
  group_by(plasmid, sample) %>%
  summarise(
    fit = list(reorder_by_mu(
      normalmixEM(
        GRN.B.HLin,
        k = k_val,
        mu = mu_init
      )
    )),
    .groups = "drop"
  )

#Create the labels for lambda
lambda_labels <- fits %>%
  mutate(
    mu     = map(fit, ~ .x$mu),
    lambda = map(fit, ~ .x$lambda)
  ) %>%
  unnest(c(mu, lambda)) %>%
  group_by(plasmid, sample) %>%
  mutate(
    component = paste0("C", row_number()),
    label = sprintf("%s: μ = %.1f | λ = %.2f", component, mu, lambda),
    x = 250,
    y = seq(0.01, 0.01 - 0.002 * (n() - 1), length.out = n())
  ) %>%
  ungroup()

# Grid
x_grid <- seq(0, 400, length.out = 500)

densities <- fits %>%
  group_by(plasmid, sample) %>%
  group_modify(~ {
    fit <- .x$fit[[1]]
    k <- length(fit$mu)
    
    x_grid <- seq(0, 400, length.out = 500)
    
    comps <- map_dfc(seq_len(k), function(i) {
      fit$lambda[i] * dnorm(x_grid, fit$mu[i], fit$sigma[i])
    })
    
    colnames(comps) <- paste0("C", seq_len(k))
    
    bind_cols(
      tibble(x = x_grid),
      comps,
      total = rowSums(comps)
    ) %>%
      pivot_longer(
        cols = starts_with("C"),
        names_to = "curve",
        values_to = "density"
      )
  }) %>%
  ungroup()

params_df <- fits %>%
  mutate(
    mu     = map(fit, "mu"),
    sigma  = map(fit, "sigma"),
    lambda = map(fit, "lambda")
  ) %>%
  select(sample, mu, sigma, lambda) %>%
  unnest(c(mu, sigma, lambda)) %>%
  group_by(sample) %>%
  mutate(
    component = paste0("C", row_number())
  ) %>%
  ungroup() %>%
  select(sample, component, mu, sigma, lambda)

lambda_df <- params_df %>%
  select(sample, component, lambda) %>%
  distinct() %>%
  pivot_wider(
    names_from  = component,
    values_from = lambda
  )

lambda_df$C1_C2 <- lambda_df$C1/lambda_df$C2
lambda_df$C2_C1 <- lambda_df$C2/lambda_df$C1
lambda_df$plasmid <- sub("_.*", "", lambda_df$sample)

lambda_df <- mutate_at(lambda_df, vars(plasmid), list(factor))
lambda_df$plasmid <- factor(lambda_df$plasmid, levels = c("ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100")) #Reorder levels

lambda_df <- lambda_df %>%
  mutate(
    State = case_when(
      plasmid == "ctrl" ~ "Single",
      plasmid == "WT" ~ "Duplicated",
      str_starts(as.character(plasmid), "pSyn") ~ "Synthetic",
    )
  )

# ANOVA
model <- lm(C2_C1 ~ plasmid, data = lambda_df)
anova(model)

# Conditions
leveneTest(C2_C1 ~ plasmid, data = lambda_df)


layout(matrix(c(1:5, 5), 2, 3))
plot(model, c(1, 3:5), pch = 20, cex = 1.5)
qqPlot(resid(model), "norm", pch = 20, cex = 1.5, line = "quartiles")
layout(1)


# Tukey test
TCM <- emmeans(model, ~ plasmid)
cld_result <- multcomp::cld(TCM, Letters = LETTERS)

# Pairwise comparisons
tukey_results <- pairs(TCM, adjust = "none")
summary(tukey_results) |>
  subset(contrast == "ctrl - pSyn100")

letters_df <- data.frame(
  plasmid = c("ctrl", "pSyn2", "pSyn30", "pSyn65", "WT", "pSyn100"),
  label = c("a", "a", "a", "ab", "ab", "b"),
  y = 3
)

#Fig5f 

Fig5f <- ggplot(lambda_df, aes(x = plasmid, y = C2_C1, fill = State)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.7) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c(
    "Single" = "#00BA38",
    "Duplicated" = "#56B1F7",
    "Synthetic" = "#FE7B09"
  ),
  breaks = c("Single", "Duplicated", "Synthetic")) +
  labs(x = "Promoter of plasmid gene copy", y = "2C/1C") +
  geom_text(data = letters_df, aes(x = plasmid, y = y, label = label),
            inherit.aes = FALSE, size = 7) +
  scale_y_continuous(limits = c(1, 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        axis.text.x = element_text(size = 16.5, color = "black"),
        axis.text.y = element_text(size = 16.5, color = "black"),
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        legend.position = c(0.9, 0.75))
Fig5f

#ggsave("Fig5f.png", plot = Fig5f, path = path_fig, width = 2250, height = 1350, units = "px", dpi = 300)


#### Fig 5e ####
data_5f_subset <- data_5f_exp %>%
  filter(plasmid == "pSyn100" & rep == "3")

# Fit mixture model (k = 3)
fit <- reorder_by_mu(
  normalmixEM(
    data_5f_subset$GRN.B.HLin,
    mu = c(80, 150, 250),  # initial guesses (adjust if needed)
    k = 3
  )
)

# Create lambda labels
lambda_labels <- tibble(
  component = paste0("C", 1:3),
  lambda = fit$lambda
) %>%
  mutate(
    label = sprintf("%s: λ = %.2f", component, lambda),
    x = 300,
    y = seq(0.018, 0.018 - 0.002 * (n() - 1), length.out = n())
  )

# Grid for densities
x_grid <- seq(0, 400, length.out = 500)

# Compute densities
densities <- {
  k <- length(fit$mu)
  
  comps <- map_dfc(seq_len(k), function(i) {
    fit$lambda[i] * dnorm(x_grid, fit$mu[i], fit$sigma[i])
  })
  
  colnames(comps) <- paste0("C", seq_len(k))
  
  bind_cols(
    tibble(x = x_grid),
    comps,
    total = rowSums(comps)
  ) %>%
    pivot_longer(-x, names_to = "curve", values_to = "density")
}

# Fig 5e
Fig5e <- ggplot(data_5f_subset, aes(x = GRN.B.HLin)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 350, length.out = 100),
    fill = "darkgrey",
    color = "white"
  ) +
  geom_line(
    data = densities %>% filter(curve != "total"),
    aes(x = x, y = density, color = curve),
    linewidth = 3
  ) +
  scale_color_manual(values = c(
    "C1" = "#199c75ff",
    "C2" = "#d95e00ff",
    "C3" = "#7570B3"
  )) +
  scale_x_continuous(limits = c(0, 350)) +
  coord_cartesian(ylim = c(0, 0.0125)) +
  labs(x = "Fluorescence (AU)", y = "Density") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 1)
  )
Fig5e

# Save
#ggsave("Fig5e.png", plot = Fig5e, path = path_fig, width = 2010, height = 1050, units = "px", dpi = 300)
