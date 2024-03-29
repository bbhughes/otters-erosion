library(ggplot2)
library(dplyr)
library(mgcv)
library(here)
theme_set(ggsidekick::theme_sleek()) # https://github.com/seananderson/ggsidekick

ribbon_col <- RColorBrewer::brewer.pal(5, "Blues")[5]

erosion <- readr::read_csv(here("data", "erosion_evd.csv"))
otter <- readr::read_csv(here("data", "sea_otter_population_data.csv"))
names(erosion) <- tolower(names(erosion))
names(otter) <- tolower(names(otter))
erosion$id <- as.factor(as.character(erosion$id))
erosion$section <- as.factor(as.character(erosion$section))
RE <- TRUE

if (RE) {
  fit_erosion <- gamm(erosion_m_y ~ s(year),
    random = list(section = ~1),
    family = gaussian(link = "identity"),
    data = erosion
  )
  fit_erosion
  plot(fit_erosion$gam)
}

if (!RE) {
  fit_erosion <- gam(erosion_m_y ~ s(year),
    family = gaussian(link = "identity"),
    data = erosion
  )
  plot(fit_erosion)
}

fit_otter <- gam(sea_otter_no ~ s(year),
  family = nb(link = "log"),
  data = otter
)
plot(fit_otter)

nd <- data.frame(year = seq(min(otter$year), max(otter$year), length.out = 300))

if (!RE) {
  pred_erosion <- predict(fit_erosion, newdata = nd)
} else {
  pred_erosion <- predict(fit_erosion$gam, newdata = nd)
}
pred_otter <- predict(fit_otter, newdata = nd)
pred <- data.frame(year = nd$year, erosion = pred_erosion, otter = pred_otter)
nd_int <- data.frame(year = otter$year)

if (!RE) {
  pred_erosion <- predict(fit_erosion, newdata = nd_int)
} else {
  pred_erosion <- predict(fit_erosion$gam, newdata = nd_int)
}
pred_otter <- predict(fit_otter, newdata = nd_int)
pred_obs <- data.frame(
  year = nd_int$year,
  erosion = pred_erosion, otter = pred_otter
)

if (!RE) {
  se_erosion <- predict(fit_erosion, newdata = nd_int, se.fit = TRUE)$se.fit
} else {
  se_erosion <- predict(fit_erosion$gam, newdata = nd_int, se.fit = TRUE)$se.fit
}
lwr_erosion <- pred_erosion - 2 * se_erosion
upr_erosion <- pred_erosion + 2 * se_erosion

se_otter <- predict(fit_otter, newdata = nd_int, se.fit = TRUE)$se.fit
lwr_otter <- exp(pred_otter - 2 * se_otter)
upr_otter <- exp(pred_otter + 2 * se_otter)

dir.create("figs", showWarnings = FALSE)

nd2 <- data.frame(year = seq(min(erosion$year), max(erosion$year), length.out = 300))
if (!RE) {
  pe <- predict(fit_erosion, newdata = nd2, se.fit = TRUE)
} else {
  pe <- predict(fit_erosion$gam, newdata = nd2, se.fit = TRUE)
}

ee <- pe$fit
elwr <- pe$fit - qnorm(0.975) * pe$se.fit
eupr <- pe$fit + qnorm(0.975) * pe$se.fit
elwr2 <- pe$fit - qnorm(0.75) * pe$se.fit
eupr2 <- pe$fit + qnorm(0.75) * pe$se.fit
de <- data.frame(year = nd2$year, est = ee, lwr = elwr, upr = eupr, lwr2 = elwr2, upr2 = eupr2)

set.seed(1)

make_ts_plot <- function(dat, point_dat, xlim, ylim, xlab = "Year", ylab = "", point_y = "erosion_m_y") {
  half_line <- 11 / 2
  ggplot(dat, aes(year, y = est, ymin = lwr, ymax = upr)) +
    geom_point(aes_string(x = "year", y = point_y), data = point_dat, inherit.aes = FALSE, alpha = 1, position = position_jitter(width = 0.80, height = 0), pch = 21, fill = "white", colour = "grey50", size = 1) +
    geom_ribbon(alpha = 0.3, fill = ribbon_col) +
    geom_ribbon(aes(ymin = lwr2, ymax = upr2), alpha = 0.4, fill = ribbon_col) +
    geom_line(lwd = 0.6, colour = ribbon_col) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ylab(ylab) +
    xlab(xlab) +
    theme(plot.margin = margin(
      t = half_line,
      r = half_line + 5, b = half_line, l = half_line
    ))
}

nd3 <- data.frame(year = seq(min(otter$year), max(otter$year), length.out = 300))
po <- predict(fit_otter, newdata = nd3, se.fit = TRUE)
oe <- exp(po$fit)
olwr <- exp(po$fit - qnorm(0.975) * po$se.fit)
oupr <- exp(po$fit + qnorm(0.975) * po$se.fit)
olwr2 <- exp(po$fit - qnorm(0.75) * po$se.fit)
oupr2 <- exp(po$fit + qnorm(0.75) * po$se.fit)
do <- data.frame(year = nd3$year, est = oe, lwr = olwr, upr = oupr, lwr2 = olwr2, upr2 = oupr2)

g1 <- make_ts_plot(de, erosion,
  xlim = c(min(erosion$year) - 1, 2020),
  ylim = quantile(erosion$erosion_m_y, probs = c(0.025, 0.975)),
  point_y = "erosion_m_y", ylab = "Erosion (m/year)"
)

mx <- pred_obs[pred_obs$erosion == max(pred_obs$erosion), ]

g3 <- ggplot() +
  geom_linerange(aes(x = exp(otter), ymin = lwr_erosion, ymax = upr_erosion, colour = year), data = pred_obs, alpha = 0.4) +
  geom_linerange(aes(xmin = lwr_otter, xmax = upr_otter, y = erosion, colour = year), data = pred_obs, alpha = 0.4) +
  scale_colour_viridis_c() +
  xlab("Sea otters") +
  ylab("Erosion (m/year)") +
  labs(colour = "Year") +
  geom_path(aes(exp(otter), erosion, colour = year), data = pred, lwd = 1.5) +
  geom_point(aes(exp(otter), erosion), data = pred_obs, inherit.aes = FALSE, size = 2, colour = "grey30") +
  geom_point(aes(exp(otter), erosion), data = pred_obs, inherit.aes = FALSE, size = 1, colour = "white") +
  # theme(legend.position = c(0.8, 0.8), legend.background = element_blank()) +
  guides(colour = "none") +
  annotate("text", x = 128, y = 0.15, label = "2018") +
  annotate("text", x = 17, y = 0.032, label = "1980s") +
  annotate("text", y = mx$erosion + 0.04, x = 31, label = mx$year) +
  coord_cartesian(expand = FALSE, xlim = c(-0.5, 138))

# mcmc panel --------------------------------------------------------------

mcmc <- readRDS("analysis/gompertz-mcmc.rds")
YEAR_START <- 8
o <- otter[seq((YEAR_START - 1), nrow(otter) - 1), "sea_otter_no", drop = TRUE]
o[is.na(o)] <- mean(c(51, 67)) # FIXME integrate over in Stan? TODO NOTE
ots <- seq(min(o), max(o), length.out = 300)
b <- lapply(mcmc, `[[`, "b") %>% unlist()
r <- lapply(mcmc, `[[`, "r") %>% unlist()
mm <- purrr::map_dfr(ots, function(.o) {
  # data.frame(otters = .o, effect = r - r * (1 - exp(-b * .o)))
  data.frame(otters = .o, effect = exp(r - r * (1 - exp(-b * .o))) - 1)
})

g_mcmc <- mm %>%
  group_by(otters) %>%
  summarise(
    lwr = quantile(effect, probs = 0.025),
    upr = quantile(effect, probs = 0.975),
    lwr2 = quantile(effect, probs = 0.25),
    upr2 = quantile(effect, probs = 0.75),
    med = quantile(effect, probs = 0.5)
  ) %>%
  ggplot(aes(otters, med)) +
  geom_line(colour = ribbon_col) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, fill = ribbon_col) +
  geom_ribbon(aes(ymin = lwr2, ymax = upr2), alpha = 0.4, fill = ribbon_col) +
  ylab("Expected erosion rate") +
  xlab("Sea otters") +
  coord_cartesian(expand = FALSE, ylim = c(0, 0.42))

# 1985-2008: Recolonization
# 2009-2012: Pre-otter expansion of tidal creeks
# 2013-2015: Otter expansion of tidal creeks
# 2016-2018: Post-otter expansion of tidal creeks

half_line <- 11 / 2
vline_col <- "#00000060"
lab_y <- 5
hjust <- 0
text_col <- "grey30"
lab_size <- 2.5
nudge <- 0.5
g2 <- ggplot(otter, aes(year, y = sea_otter_no)) +
  # geom_vline(xintercept = 1985, col = vline_col) +
  geom_vline(xintercept = 2008.5, col = vline_col) +
  geom_vline(xintercept = 2012.5, col = vline_col) +
  geom_vline(xintercept = 2015.5, col = vline_col) +
  annotate("text", x = 1998, y = lab_y + 2, hjust = hjust, angle = 0, col = text_col, label = "Recolonization", size = lab_size) +
  annotate("text", x = 2010 - nudge, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "Pre-otter", size = lab_size) +
  annotate("text", x = 2011.2 - nudge, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "expansion", size = lab_size) +
  annotate("text", x = 2014 - nudge - 0.2, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "Otter", size = lab_size) +
  annotate("text", x = 2015.2 - nudge - 0.2, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "expansion", size = lab_size) +
  annotate("text", x = 2017.5 - nudge, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "Post-otter", size = lab_size) +
  annotate("text", x = 2018.7 - nudge, y = lab_y, hjust = hjust, angle = 90, col = text_col, label = "expansion", size = lab_size) +
  geom_ribbon(aes(year, ymin = lwr, ymax = upr), fill = ribbon_col, data = do, inherit.aes = FALSE, alpha = 0.3) +
  geom_ribbon(aes(year, ymin = lwr2, ymax = upr2), fill = ribbon_col, data = do, inherit.aes = FALSE, alpha = 0.4) +
  geom_line(aes(year, y = est), colour = ribbon_col, lwd = 0.8, data = do, inherit.aes = FALSE) +
  geom_point(aes_string(x = "year", y = "sea_otter_no"), data = otter, inherit.aes = FALSE, alpha = 1, pch = 21, fill = "white", colour = "grey30", size = 1.5) +
  coord_cartesian(xlim = c(min(otter$year) - 1, 2020), ylim = c(0, max(otter$sea_otter_no)), expand = FALSE) +
  ylab("Sea otters") +
  xlab("Year") +
  theme(plot.margin = margin(t = half_line, r = half_line + 5, b = half_line, l = half_line))
g2

g <- cowplot::plot_grid(g2, g1, g3, g_mcmc, ncol = 2)
ggsave("figs/fig1-lower.png", width = 6.6, height = 4.4, dpi = 300)
ggsave("figs/fig1-lower.pdf", width = 6.6, height = 4.4)

# Values for paper --------------------------------------------------------

# expected rate of widening at 0 and 100:
quantile(r - r * (1 - exp(-b * 0)), probs = c(0.05, 0.5, 0.95))
quantile(r - r * (1 - exp(-b * 100)), probs = c(0.05, 0.5, 0.95))

# true percentage: (exp(x) - 1)
round(quantile(exp(r - r * (1 - exp(-b * 0))) - 1, probs = c(0.05, 0.5, 0.95)), 2)
round(quantile(exp(r - r * (1 - exp(-b * 100))) - 1, probs = c(0.05, 0.5, 0.95)), 2)

round(quantile(exp(r - r * (1 - exp(-b * 0))) - 1, probs = c(0.025, 0.5, 0.975)), 2)
round(quantile(exp(r - r * (1 - exp(-b * 100))) - 1, probs = c(0.025, 0.5, 0.975)), 2)
