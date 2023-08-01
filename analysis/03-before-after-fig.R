### Data analysis for testing sea otter impacts on crabs, salt marsh and tidal creek bank erosion####
# By Brent B. Hughes et al. 2023#

# Load libraries
library(Rmisc)
library(ggplot2)
library(dplyr)
library(glmmTMB)

dir.create("figs", showWarnings = FALSE)

sig_size <- 5

#### Before/After Comparison during the most recent otter expansion in Elkhorn Slough####
#### GLMM tweedie analysis of sea otters in creeks, crab consumption, and erosion across 13 creeks####
before_after_sea_otter <- read.csv("data/before_after_sea_otter.csv")
str(before_after_sea_otter)

# convert year to factor
before_after_sea_otter$Period <- as.factor(as.character(before_after_sea_otter$Period))

# Sea Otters

# Run tweedie analysis
fit_otter <- glmmTMB(
  Otter_per_ha ~ Period + (1 | Creek),
  data = before_after_sea_otter,
  # dispformula = ~ Period,
  family = tweedie(),
  # verbose = TRUE
)
summary(fit_otter)
# dispersion parameter = 0.078, df residuals = 21, p-value < 0.0005, significant difference between before/after sea otter expansion.
r_otter <- DHARMa::simulateResiduals(fit_otter)
plot(r_otter)

# figure of otters in creeks before and after sea otter colonization of marshes

# bar graph of Before/After Otters
# run this once at the top to set the theme:
# from remotes::install_github('seananderson/ggsidekick')

source("analysis/theme_sleek.R")
theme_set(theme_sleek())

before_after_sea_otter$x <- NULL # in case you've run the code below already
before_after_sea_otter$xend <- NULL
segment_dat <- before_after_sea_otter |>
  # keep your variable + period + creek:
  select(Otter_per_ha, Period, Creek) |>
  # change values_from to your variable:
  tidyr::pivot_wider(names_from = Period, values_from = Otter_per_ha)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
before_after_sea_otter$x <- c(segment_dat$x, segment_dat$xend)

# generate standard errors: -------------------------------
baso_summary <- summarySE(before_after_sea_otter, measurevar = "Otter_per_ha", groupvars = c("Period"))
baso_summary

# make the plot: ------------------------------------------
baso_plot <- ggplot() +
  geom_bar(
    data = baso_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(Period, level = c("Before_otter", "After_otter"))),
      y = Otter_per_ha, fill = Period
    )
  ) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = Before_otter, yend = After_otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = before_after_sea_otter, mapping = aes(x = x, y = Otter_per_ha),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  labs(y = "Sea otters\nper ha of tidal creek", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("A") +
  annotate("text", x = 0.84, y = 0.55, label = "***", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  theme(axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank())
baso_plot

# ggsave("figs/saso-bar.png", width = 3.2, height = 3)

# Pachygrapsus consumed per day
fit_pachy <- glmmTMB(
  (Pachy_consumed_per_day) ~ Period + (1 | Creek),
  data = before_after_sea_otter,
  family = tweedie()
  # verbose = TRUE
)
summary(fit_pachy)
# dispersion parameter = 0.678, df residuals = 21, p-value < 0.0005, significant difference between before and after sea otter expansion.
r_pachy <- DHARMa::simulateResiduals(fit_pachy)
plot(r_pachy)

# bar graph of Before/After Pachygrapsus consumed
before_after_sea_otter$x <- NULL # in case you've run the code below already
before_after_sea_otter$xend <- NULL
segment_dat <- before_after_sea_otter |>
  # keep your variable + period + creek:
  select(Pachy_consumed_per_day, Period, Creek) |>
  # change values_from to your variable:
  tidyr::pivot_wider(names_from = Period, values_from = Pachy_consumed_per_day)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
before_after_sea_otter$x <- c(segment_dat$x, segment_dat$xend)

# generate standard errors: -------------------------------
baso_pachy_summary <- summarySE(before_after_sea_otter, measurevar = "Pachy_consumed_per_day", groupvars = c("Period"))
baso_pachy_summary

# make the plot: ------------------------------------------
baso_pachy_plot <- ggplot() +
  geom_bar(
    data = baso_pachy_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(Period, level = c("Before_otter", "After_otter"))),
      y = Pachy_consumed_per_day, fill = Period
    )
  ) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = Before_otter, yend = After_otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = before_after_sea_otter, mapping = aes(x = x, y = Pachy_consumed_per_day),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("", "")
  ) +
  labs(y = "Crab consumed\nper ha per day", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  ggtitle("B") +
  annotate("text", x = 0.84, y = 40, label = "***", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  theme(axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank())
baso_pachy_plot

# ggsave("figs/pachy_eaten-bar.png", width = 3.2, height = 3)

# Before After comparison of annual erosion in 13 creek banks

# Run LMM
fit_erosion <- glmmTMB(
  erosion_m_yr ~ Period + (1 | Creek),
  data = before_after_sea_otter,
  family = gaussian()
  # verbose = TRUE
)
summary(fit_erosion)
# dispersion parameter for gaussian (sigma^2) = 0.0044, df residuals = 22, p-value < 0.0005, significant difference between 2011 and 2017.
r_erosion <- DHARMa::simulateResiduals(fit_erosion)
plot(r_erosion)

# bar graph of Before/After creek bank erosion
before_after_sea_otter$x <- NULL # in case you've run the code below already
before_after_sea_otter$xend <- NULL
segment_dat <- before_after_sea_otter |>
  # keep your variable + period + creek:
  select(erosion_m_yr, Period, Creek) |>
  # change values_from to your variable:
  tidyr::pivot_wider(names_from = Period, values_from = erosion_m_yr)

# jitter the x values a bit: ------------------------------
set.seed(1)
# control jitter amount here:
segment_dat$x <- jitter(rep(1, nrow(segment_dat)), amount = 0.01)
segment_dat$xend <- jitter(rep(2, nrow(segment_dat)), amount = 0.01)
before_after_sea_otter$x <- c(segment_dat$x, segment_dat$xend)

# generate standard errors: -------------------------------
baso_erosion_summary <- summarySE(before_after_sea_otter, measurevar = "erosion_m_yr", groupvars = c("Period"))
baso_erosion_summary

# make the plot: ------------------------------------------
ba_erosion_plot <- ggplot() +
  geom_bar(
    data = baso_erosion_summary, stat = "identity", width = .5,
    mapping = aes(
      x = as.numeric(factor(Period, level = c("Before_otter", "After_otter"))),
      y = erosion_m_yr, fill = Period
    )
  ) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  geom_segment(data = segment_dat, aes(
    x = x, xend = xend,
    y = Before_otter, yend = After_otter
  ), inherit.aes = FALSE, alpha = 0.15) +
  geom_point(
    data = before_after_sea_otter, mapping = aes(x = x, y = erosion_m_yr),
    pch = 21, fill = "grey85", colour = "grey40", size = 2
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Otter\nExpansion", "Post-Otter\nExpansion")
  ) +
  labs(y = "Bank erosion\n(m per year)", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.04, .04))) +
  ggtitle("C") +
  annotate("text", x = 0.84, y = 0.4, label = "***", col = "grey30", size = sig_size) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, lty = 2, col = "grey40") # SA added
ba_erosion_plot

# ggsave("figs/ba_erosion-bar.png", width = 3.2, height = 3)

# Save Fig. Before After
# pdf("ba_plot.pdf", width=6, height=13)
# gridExtra::grid.arrange(baso_plot, baso_pachy_plot, ba_erosion_plot, ncol = 1, nrow = 3)
# dev.off()
cowplot::plot_grid(baso_plot, baso_pachy_plot, ba_erosion_plot, ncol = 1, align = "v", rel_heights = c(1, 1, 1.25))
ggsave("figs/ba_plot.pdf", width = 3.2, height = 6)
ggsave("figs/ba_plot.png", width = 3.2, height = 6, dpi = 300)
