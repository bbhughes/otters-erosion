library(ggplot2)
library(dplyr)
library(glmmTMB)
dir.create("figs", showWarnings = FALSE)

source("analysis/theme_sleek.R")
theme_set(theme_sleek())

sig_size <- 5

#### Experimental Data analyses####

# Load data
Otter_marsh_data <- read.table(file = "data/Otters_Marsh_experimental_data_R.csv", header = TRUE, sep = ",")
# str(Otter_marsh_data)
# View(Otter_marsh_data)

#### Feeding Assay####
Otter_marsh_data$AG_Consum_Crab <- as.numeric(as.character(Otter_marsh_data$AG_Consum_Crab))
Otter_marsh_data$AG_Consum_NoCrab <- as.numeric(as.character(Otter_marsh_data$AG_Consum_NoCrab))
Otter_marsh_data$BG_Consum_Crab <- as.numeric(as.character(Otter_marsh_data$BG_Consum_Crab))
Otter_marsh_data$BG_Consum_NoCrab <- as.numeric(as.character(Otter_marsh_data$BG_Consum_NoCrab))

# Aboveground feeding trial
var.test(Otter_marsh_data$AG_Consum_Crab, Otter_marsh_data$AG_Consum_NoCrab, data = two)
t.test(Otter_marsh_data$AG_Consum_Crab, Otter_marsh_data$AG_Consum_NoCrab, data = two, var.equal = T)
# t = -0.30593, df = 18, p-value = 0.7632, not significant

# Belowground feeding trial
var.test(Otter_marsh_data$BG_Consum_Crab, Otter_marsh_data$BG_Consum_NoCrab, data = two)
t.test(Otter_marsh_data$BG_Consum_Crab, Otter_marsh_data$BG_Consum_NoCrab, data = two, var.equal = T)
# t = 3.3884, df = 18, p-value = 0.003275, significant

# Plot the feeding assay data
# Load Figure data
Otter_marsh_data_assay_fig <- read.table(file = "data/Otters_Marsh_2015_R_assay_fig_data.csv", header = TRUE, sep = ",")
# str(Otter_marsh_data_assay_fig)

# plot the data
assay_plot <- ggplot(Otter_marsh_data_assay_fig, aes(x = Marsh_type, y = biomass_consumed, fill = Treatment, group = paste(Treatment, Marsh_type))) +
  geom_violin(position = position_dodge(1)) +
  geom_dotplot(
    binaxis = "y", stackdir = "center",
    position = position_dodge(1)
  ) +
  # theme_bw() +
  theme(legend.position = c(0.15, 0.85)) +
  annotate("text", x = 0.75, y = 0.35, label = "a") +
  annotate("text", x = 1.25, y = 0.45, label = "a") +
  annotate("text", x = 1.75, y = 0.37, label = "x") +
  annotate("text", x = 2.15, y = 0.5, label = "y") +
  xlab("") +
  ylab(expression("Marsh biomass consumed (g (fw) per trial)"))
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) +
  # theme(axis.text.x = element_text())
assay_plot
# ggsave("figs/extended-fig4.jpg", width = 5, height = 4)
ggsave("figs/extended-fig4.png", width = 5, height = 4)
ggsave("figs/extended-fig4.pdf", width = 5, height = 4)
ggsave("figs/extended-fig4.eps", width = 5, height = 4)
# Fig. S2

#### Field Experimental Data analyses 2014-2016####
# First, run preliminary data analysis comparing burrows (all burrows greater than 1 cm) and marsh percent cover between treatments
# burrows
prelim_exp_data <- read.table(file = "data/Otter_experiment_preliminary_2013.csv", header = TRUE, sep = ",")
# str(prelim_exp_data)

pre_burrow_glmm_nb2 <- glmmTMB(burrows_1_3cm_pre ~ Treatment + (1 | plot_id) + (1 | Site), data = prelim_exp_data, family = nbinom2())
summary(pre_burrow_glmm_nb2)
# P = 0.623, dispersion paramter = 0.48, df residual = 46, no difference in burrows for pre-experimental test

# marsh_percent_cover
# Run LMM
fit_pre_percent_cover <- glmmTMB(
  asin(sqrt(per_cover_marsh)) ~ Treatment + (1 | plot_id) + (1 | Site),
  data = prelim_exp_data,
  family = gaussian()
  # verbose = TRUE
)
summary(fit_pre_percent_cover)
# dispersion parameter for gaussian (sigma^2) = 0.0338, df residuals = 45, p-value =0.488, no preliminary differences.
# r_pre_percent_cover <- DHARMa::simulateResiduals(fit_pre_percent_cover)
# plot(r_pre_percent_cover)

# Run 2014-2016 Experimental Analysis

# subset data to select caged (no otter) and procedural controls (otter)
Otter_marsh_data2 <- subset(Otter_marsh_data, treat == "no_otter" | treat == "otter")
# str(Otter_marsh_data2)

Otter_marsh_data2$AGMass_kg_msq <- as.numeric(as.character(Otter_marsh_data2$AGMass_kg_msq))
Otter_marsh_data2$BGMass_kg_msq <- as.numeric(as.character(Otter_marsh_data2$BGMass_kg_msq))
Otter_marsh_data2$BulkDens_kg_msq <- as.numeric(as.character(Otter_marsh_data2$BulkDens_kg_msq))
Otter_marsh_data2$crab_density_msq <- as.integer(as.character(Otter_marsh_data2$crab_density_msq))
Otter_marsh_data2$Burrows_change <- as.integer(as.character(Otter_marsh_data2$Burrows_change))
Otter_marsh_data2$sed_accrete_g <- as.numeric(as.character(Otter_marsh_data2$sed_accrete_g))
Otter_marsh_data2$treat <- as.factor(as.character(Otter_marsh_data2$treat))
Otter_marsh_data2$Site <- as.factor(as.character(Otter_marsh_data2$Site))

# View(Otter_marsh_data2)

# glmm for change in burrows from 2014-2015 with site and plot ID as a random factor
burrow_lme <- lme4::lmer(Burrows_change ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2)
summary(burrow_lme)

burrow_lme <- glmmTMB::glmmTMB(Burrows_change ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2)
summary(burrow_lme)

# AIC(burrow_glm, burrow_glm_nb1, burrow_glm_p, burrow_glm_qp, burrow_glm_gaus)
# summary(burrow_glm_gaus)
# non Significant Otter effect, P = 0.0164
# decrease in crab burrows with sea otters present

# r_burrow_glm <- DHARMa::simulateResiduals(burrow_glm_gaus)
# plot(r_burrow_glm)

ggplot(Otter_marsh_data2, aes(x = treat, y = Burrows_change)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

tapply(Otter_marsh_data2$Burrows_change, Otter_marsh_data2$treat, FUN = mean)
crab_dec <- ((2.12 - 1.04) / ((2.12 + 1.04) / 2)) * 100
crab_dec # 68.4% increase in crabs without otters

# glmm for crab density (collected in 2016) with site and plot as random factors
crab_glmm <- glmmTMB(crab_density_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = nbinom2(), na.action = na.omit)
summary(crab_glmm) # significant otter effect, P=0.0157

# decrease in crab density with sea otters present
tapply(Otter_marsh_data2$crab_density_msq, Otter_marsh_data2$treat, FUN = mean)
crab_dec <- ((2.12 - 1.04) / ((2.12 + 1.04) / 2)) * 100
crab_dec # 68.4% increase in crabs without otters

# r_crab_glm <- DHARMa::simulateResiduals(crab_glm)
# plot(r_crab_glm)

ggplot(Otter_marsh_data2, aes(x = treat, y = crab_per_trap)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

# glmm for bulk density collected in 2015
bulkdens_glmm_gamma <- glmmTMB(BulkDens_kg_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = Gamma(link = "log"), na.action = na.omit)
summary(bulkdens_glmm_gamma)

bulkdens_glmm_tweedie <- glmmTMB(BulkDens_kg_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = tweedie(link = "log"), na.action = na.omit)
summary(bulkdens_glmm_tweedie)

bulkdens_lmm <- glmmTMB(BulkDens_kg_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = gaussian(link = "identity"), na.action = na.omit)
summary(bulkdens_lmm)

AIC(bulkdens_glmm_gamma, bulkdens_glmm_tweedie, bulkdens_lmm)
# bulk_dens_lmm is the top selected model
# P = 0.0082

# r_bulkdens_lmm <- DHARMa::simulateResiduals(bulkdens_lmm)
# plot(r_bulkdens_lmm)
# Not normally distributed (KS Test, p = 0.0366), however still the best fitting model that meets most of the assumptions
# transformations do not improve the model

tapply(Otter_marsh_data2$BulkDens_kg_msq, Otter_marsh_data2$treat, FUN = mean, na.rm = TRUE)
bulk_inc <- ((126.0 - 136.5) / ((126.0 + 136.5) / 2)) * 100
bulk_inc
# 8% difference in bulk density

ggplot(Otter_marsh_data2, aes(x = treat, y = BulkDens_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

# glmm for belowground biomass

ggplot(Otter_marsh_data2, aes(x = treat, y = BGMass_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

BG_mass_glmm <- glmmTMB(BGMass_kg_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = Gamma(link = "log"), na.action = na.omit)
summary(BG_mass_glmm)

# r_BG_mass_glmm <- DHARMa::simulateResiduals(BG_mass_glmm)
# plot(r_BG_mass_glmm)

# Significant otter effect, P = 0.035
# Otters had 0.46 kg more BG biomass per msq
# Difference in BG biomass with sea otters present
tapply(Otter_marsh_data2$BGMass_kg_msq, Otter_marsh_data2$treat, FUN = mean, na.rm = TRUE)
BG_inc <- ((2.98 - 3.46) / ((2.98 + 3.46) / 2)) * 100
BG_inc # 14.9% Difference in BG biomass

# glmm for aboveground biomass

ggplot(Otter_marsh_data2, aes(x = treat, y = AGMass_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

AG_mass_glmm <- glmmTMB(AGMass_kg_msq ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = Gamma(link = "log"), na.action = na.omit)
summary(AG_mass_glmm)

# r_AG_mass_glm <- DHARMa::simulateResiduals(AG_mass_glm)
# plot(r_AG_mass_glm)

# Significant otter effect, P = 0.000905
# Otters had 0.37 kg more AG biomass per msq
# difference in bulk density with sea otters present
tapply(Otter_marsh_data2$AGMass_kg_msq, Otter_marsh_data2$treat, FUN = mean, na.action = na.omit)
AG_inc <- ((0.65 - 1.06) / ((0.65 + 1.06) / 2)) * 100
AG_inc # Difference 48.0%

# Plot the data
# Burrows bargraph

segment_dat <- Otter_marsh_data2 |>
  select(Burrows_change, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = Burrows_change)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
Otter_marsh_data2$x <- c(segment_dat$x, segment_dat$xend)
plot(Otter_marsh_data2$x, Otter_marsh_data2$treat)

# generate standard errors: -------------------------------
burrow_change_summary <- Rmisc::summarySE(Otter_marsh_data2, measurevar = "Burrows_change", groupvars = c("treat"))
burrow_change_summary

# make the plot: ------------------------------------------
p1 <- ggplot() +
  geom_bar(
    data = burrow_change_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = Burrows_change, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data2, mapping = aes(x = x, y = Burrows_change),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  ylab(expression(atop("2014-2015 change", "in burrows" %.% m^-2))) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.04, .04))) +
  ggtitle("C") +
  annotate("text", x = 0.85, y = 40, label = "*", col = "grey30", size = sig_size) +
  theme(legend.position = "none")  + theme(axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank()) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50")
p1

# crab density bargraph
segment_dat <- Otter_marsh_data2 |>
  select(crab_density_msq, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = crab_density_msq)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.04)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.04)
Otter_marsh_data2$x <- c(segment_dat$x, segment_dat$xend)
Otter_marsh_data2$crab_density_msq <- jitter(Otter_marsh_data2$crab_density_msq, amount = 0.1)

# generate standard errors: -------------------------------
crab_density_msq_summary <- Rmisc::summarySE(Otter_marsh_data2, measurevar = "crab_density_msq", groupvars = c("treat"))
crab_density_msq_summary

p2 <- ggplot() +
  geom_bar(
    data = crab_density_msq_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = crab_density_msq, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data2, mapping = aes(x = x, y = crab_density_msq),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("No Otter", "Otter")
  ) +
  xlab("") +
  ylab(expression(Crabs %.% "2 m"^-2)) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("A") +
  annotate("text", x = 2.0, y = 8, label = "*", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter", "Otter"))
p2

# Aboveground biomass, barplot
# remove one row with NA and the other plot in the block
Otter_marsh_data3 <- Otter_marsh_data2[-c(25, 50), ]
# str(Otter_marsh_data3)
Otter_marsh_data3$AGMass_kg_msq <- as.numeric(as.character(Otter_marsh_data3$AGMass_kg_msq))

Otter_marsh_data3$x <- NULL # in case you've run the code below already
Otter_marsh_data3$xend <- NULL
segment_dat <- Otter_marsh_data3 |>
  select(AGMass_kg_msq, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = AGMass_kg_msq)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
Otter_marsh_data3$x <- c(segment_dat$x, segment_dat$xend)
plot(Otter_marsh_data3$x, Otter_marsh_data3$treat)

# generate standard errors: -------------------------------
AGMass_kg_msq_summary <- Rmisc::summarySE(Otter_marsh_data3, measurevar = "AGMass_kg_msq", groupvars = c("treat"))
AGMass_kg_msq_summary

p3 <- ggplot() +
  geom_bar(
    data = AGMass_kg_msq_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = AGMass_kg_msq, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data3, mapping = aes(x = x, y = AGMass_kg_msq),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  xlab("") +
  ylab(expression(atop("Aboveground mass", (kg(dw) %.% m^-2)))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("D") +
  annotate("text", x = 0.85, y = 2.2, label = "**", col = "grey30", size = sig_size) +
  theme(legend.position = "none")  + theme(axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank())
p3


# belowground biomass, bargraph
Otter_marsh_data2$BGMass_kg_msq <- as.numeric(as.character(Otter_marsh_data2$BGMass_kg_msq))

segment_dat <- Otter_marsh_data2 |>
  select(BGMass_kg_msq, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = BGMass_kg_msq)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
Otter_marsh_data2$x <- c(segment_dat$x, segment_dat$xend)
plot(Otter_marsh_data2$x, Otter_marsh_data2$treat)

# generate standard errors: -------------------------------
BGMass_kg_msq_summary <- Rmisc::summarySE(Otter_marsh_data2, measurevar = "BGMass_kg_msq", groupvars = c("treat"))
BGMass_kg_msq_summary

p4 <- ggplot() +
  geom_bar(
    data = BGMass_kg_msq_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = BGMass_kg_msq, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data2, mapping = aes(x = x, y = BGMass_kg_msq),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  xlab("") +
  ylab(expression(atop("Belowground mass", (kg(dw) %.% m^-2)))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("E") +
  annotate("text", x = 0.85, y = 6, label = "*", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("No Otter", "Otter"))
p4

# Bulk density, bargraph
Otter_marsh_data2$BulkDens_kg_msq <- as.numeric(as.character(Otter_marsh_data2$BulkDens_kg_msq))

segment_dat <- Otter_marsh_data2 |>
  select(BulkDens_kg_msq, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = BulkDens_kg_msq)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
Otter_marsh_data2$x <- c(segment_dat$x, segment_dat$xend)
plot(Otter_marsh_data2$x, Otter_marsh_data2$treat)

# generate standard errors: -------------------------------
BulkDens_kg_msq_summary <- Rmisc::summarySE(Otter_marsh_data2, measurevar = "BulkDens_kg_msq", groupvars = c("treat"))
BulkDens_kg_msq_summary

p5 <- ggplot() +
  geom_bar(
    data = BulkDens_kg_msq_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = BulkDens_kg_msq, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data2, mapping = aes(x = x, y = BulkDens_kg_msq),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  xlab("") +
  ylab(expression(atop("Bulk density", (kg(dw) %.% m^-2)))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("B") +
  annotate("text", x = 0.75, y = 150, label = "*", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter", "Otter"))
p5

# Make Experimental figure
cowplot::plot_grid(p1, p3, p4, ncol = 1, align = "v", rel_heights = c(1, 1, 1.20))
ggsave("figs/experiment_figure_aligned.pdf", width = 3.2, height = 6)

# extended fig ----------------------------------------------------------------

# Comparison of burrows in procedural controls (Otter) with undisturbed transect plots at each experimental creek
# subset data to select procedural controls (otter) and transect plots
Otter_marsh_data4 <- subset(Otter_marsh_data, treat == "otter_transect" | treat == "otter")
# str(Otter_marsh_data4)

Otter_marsh_data4$Burrows_msq_2015 <- as.integer(as.character(Otter_marsh_data4$Burrows_msq_2015))
transect_burrow_glmm_nb2 <- glmmTMB(Burrows_msq_2015 ~ treat + (1 | Site), data = Otter_marsh_data4, family = nbinom2())
summary(transect_burrow_glmm_nb2)

# variance components collapsed to zero
# diagnose(burrow_glmm)

transect_burrow_glm <- glmmTMB(Burrows_msq_2015 ~ treat, data = Otter_marsh_data4, family = nbinom2())
summary(transect_burrow_glm)

transect_burrow_glm_nb1 <- glmmTMB(Burrows_msq_2015 ~ treat, data = Otter_marsh_data4, family = nbinom1())
summary(transect_burrow_glm_nb1)

transect_burrow_glm_p <- glmmTMB(Burrows_msq_2015 ~ treat, data = Otter_marsh_data4, family = poisson())
summary(transect_burrow_glm_p)

transect_burrow_glm_qp <- glm(Burrows_msq_2015 ~ treat, data = Otter_marsh_data4, family = quasipoisson())
summary(transect_burrow_glm_qp)

AIC(transect_burrow_glm, transect_burrow_glm_nb1, transect_burrow_glm_p, transect_burrow_glm_qp)
summary(transect_burrow_glm)
# no significant differences, P = 0.337

# Bargraph, Burrow comparison between procedural controls and transects at the 5 experimental sites

# Otter_marsh_data4 <- na.omit(Otter_marsh_data4)
# Otter_marsh_data4 <- dplyr::filter(Otter_marsh_data4, !is.na(Burrows_msq_2015), !is.na(treat))
# str(Otter_marsh_data4)

# Otter_marsh_data4 %>%
#   dplyr::group_by(treat) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L)

# segment_dat <- Otter_marsh_data4 |>
#   select(Burrows_msq_2015, treat, Site, Block) |>
#   tidyr::pivot_wider(names_from = treat, values_from = Burrows_msq_2015)
#
# segment_dat <- na.omit(segment_dat)

# jitter the x values a bit: ------------------------------
# set.seed(1)
# # control jitter amount here:
# segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
# segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
# Otter_marsh_data4$x <- c(segment_dat$x, segment_dat$xend)

# d <- tidyr::pivot_longer(segment_dat, cols = c(otter, otter_transect), values_to = "Burrows_msq_2015", names_to = "treat")

# generate standard errors: -------------------------------

Otter_marsh_data4 <- filter(Otter_marsh_data4, !is.na(Burrows_msq_2015))
burrow_compar_summary <- Rmisc::summarySE(Otter_marsh_data4, measurevar = "Burrows_msq_2015", groupvars = c("treat"))
burrow_compar_summary

# make the plot: ------------------------------------------
p6 <- ggplot() +
  geom_bar(
    data = burrow_compar_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("otter_transect", "otter"))),
      y = Burrows_msq_2015, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#3182bd", "#3182bd")) +
  geom_point(
    data = Otter_marsh_data4, mapping = aes(x = treat, y = Burrows_msq_2015),
    pch = 21, fill = "grey85", colour = "grey40", size = 2,
    position = position_jitter(width = 0.02, height = 0)
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  ylab(expression("Crab burrows" %.% m^-2)) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("D") +
  annotate("text", x = 0.75, y = 65, label = "NS", col = "grey30") +
  theme(legend.position = "none") +
  scale_x_discrete(limit = c("otter", "otter_transect"), labels = c("Otter\n(procedural control)", "Otter\n(no manipulation)"))
p6

# Comparison of sediment accretion in procedural controls (Otter) with undisturbed transect plots at each experimental creek
# First run a linear model
accret_lm <- lm(sed_accrete_g ~ treat, data = Otter_marsh_data2, na.action = na.omit)
# plot(accret_lm)
# data appear normally distributed and equal variances

# Run full linear mixed model
accret_lmm <- glmmTMB(sed_accrete_g ~ treat + (1 | plot_id) + (1 | Site), data = Otter_marsh_data2, family = gaussian(), na.action = na.omit)
summary(accret_lmm)
# P = 0.877, no significant differences in sediment accretion rates

# sedimentation, bargraph

Otter_marsh_data2a <- filter(Otter_marsh_data2, !is.na(sed_accrete_g))

segment_dat <- Otter_marsh_data2a |>
  select(sed_accrete_g, treat, plot_id) |>
  tidyr::pivot_wider(names_from = treat, values_from = sed_accrete_g)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
Otter_marsh_data2a$x <- c(segment_dat$x, segment_dat$xend)
plot(Otter_marsh_data2a$x, Otter_marsh_data2a$treat)

# generate standard errors: -------------------------------
# Otter_marsh_data2 <- na.omit(Otter_marsh_data2)
sedimentation_summary <- Rmisc::summarySE(Otter_marsh_data2a, measurevar = "sed_accrete_g", groupvars = c("treat"))
sedimentation_summary

p7 <- ggplot() +
  geom_bar(
    data = sedimentation_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(treat, level = c("no_otter", "otter"))),
      y = sed_accrete_g, fill = treat
    )
  ) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = no_otter, yend = otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = Otter_marsh_data2a, mapping = aes(x = x, y = sed_accrete_g),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  xlab("") +
  ylab(expression(atop("Sediment accretion", "(g DW per year)"))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("C") +
  annotate("text", x = 0.75, y = max(segment_dat$no_otter) * 0.9, label = "NS", col = "grey30") +
  theme(legend.position = "none") +
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter", "Otter"))
p7

# Make Extended Experimental figure

half_line <- 11 / 2
m <- theme(plot.margin = margin(t = half_line - 3, r = half_line, b = half_line - 15, l = half_line - 1))

cowplot::plot_grid(
  p2 + m,
  p7 + m,
  p5 + m,
  p6 + m,
  ncol = 2, nrow = 2, align = "hv") + theme(plot.margin = margin(0, 0, -5, 0))

# ggsave("figs/extended-fig1.eps", width = 6.5, height = 5)
ggsave("figs/extended-fig1.jpg", width = 6.5, height = 5, quality = 85)
ggsave("figs/extended-fig1.pdf", width = 6.5, height = 5)
