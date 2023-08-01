library(ggplot2)
library(dplyr)
library(mgcv)
library(readr)
library(lubridate)
dir.create("figs", showWarnings = FALSE)

source("analysis/theme_sleek.R")
theme_set(theme_sleek())

d <- readr::read_csv("data/crab_change_2014_to_2016.csv", na = "na")

names(d) <- tolower(names(d))
d$sample_date <- lubridate::mdy(d$date)
d$sample_date_group <- paste(lubridate::year(d$sample_date),
  lubridate::month(d$sample_date),
  sep = "-"
)
d <- dplyr::rename(d, crabs = total_number_of_crabs_trapped)

ggplot(d, aes(sample_date_group, crabs, colour = treatment)) +
  geom_boxplot()
ggplot(d, aes(sample_date_group, change_in_crab, colour = treatment)) +
  geom_boxplot()

d$block <- as.factor(d$block)
d$site <- as.factor(d$site)
d$plot_id <- factor(paste(d$site, d$block, d$treatment, sep = "-"))
d$block_id <- factor(paste(d$site, d$block, sep = "-"))

# sanity check:
n_rows <- nrow(d)
n_plots <- length(unique(d$plot_id))
n_dates <- length(unique(d$sample_date_group))
stopifnot(n_rows == n_plots * n_dates)

d$treatment <- factor(d$treatment, levels = c("Otter", "No Otter"))
d$julian_day <- d$julian_day - min(d$julian_day) + 0
d$month <- lubridate::month(d$sample_date)
d$dum <- 1

d <- ungroup(d) |>
  group_by(sample_date_group) |>
  mutate(average_julian_date = mean(julian_day)) |>
  ungroup()

init <- dplyr::filter(d, sample_date_group == "2014-5") |>
  dplyr::rename(crabs_init = crabs) |>
  select(treatment, block_id, crabs_init) |>
  distinct()

d <- left_join(d, init)
d <- mutate(d, change_in_crab2 = crabs - crabs_init)

fit <- gam(
  change_in_crab2 ~
    s(julian_day, k = 3) +
    treatment +
    s(month, bs = "cc", k = 4) +
    s(block_id, bs = "re", by = dum),
  knots = list(month = seq(0.5, 12.5, length.out = 4)),
  data = d,
  family = gaussian(), method = "REML"
)

summary(fit)
par(mfrow = c(2, 2))
plot(fit)

par(mfrow = c(2, 2))
gam.check(fit)

par(mfrow = c(1, 1))
x <- plot(fit, select = 1)
dat <- data.frame(est = x[[1]]$fit, julian_day = x[[1]]$x, se = x[[1]]$se)
trans <- as.numeric
dat$lwr <- trans(dat$est - dat$se * 1.96)
dat$upr <- trans(dat$est + dat$se * 1.96)
dat$est <- trans(dat$est)
dplot <- filter(d, sample_date_group != "2014-5")

gday <- ggplot(dat, aes(julian_day, est, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) + #+ geom_hline(yintercept = 0, lty = 1)
  ylab("s(Day)") +
  xlab("Sample date") +
  theme(axis.text.x.bottom = element_text(angle = 90)) +
  scale_x_continuous(
    breaks = unique(dplot$average_julian_date),
    labels = unique(dplot$sample_date_group),
    expand = expansion(mult = c(0, 0))
  ) +
  theme(axis.title.y = ggtext::element_markdown(), axis.text.x.bottom = element_text(size = 8))
gday

x <- plot(fit, select = 2, xlim = c(1, 12))
datm <- data.frame(est = x[[2]]$fit, month = x[[2]]$x, se = x[[2]]$se)
trans <- as.numeric
datm$lwr <- trans(datm$est - datm$se * 2)
datm$upr <- trans(datm$est + datm$se * 2)
datm$est <- trans(datm$est)

gmonth <- ggplot(datm, aes(month, est, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  ylab("s(month)") +
  xlab("Sample date (arranged by month)") +
  scale_x_continuous(
    breaks = unique(dplot$month),
    labels = unique(dplot$sample_date_group),
    expand = expansion(mult = c(0, 0))
  ) +
  geom_vline(xintercept = 1:12, col = "grey50", alpha = 0.2) +
  theme(axis.text.x.bottom = element_text(angle = 90, size = 8))
gmonth
nd <- filter(d, sample_date_group != "2014-5") |>
  select(month, julian_day, treatment) |>
  distinct()

nd$block_id <- unique(d$block_id)[1]
nd$dum <- 0
p <- predict(fit, se.fit = TRUE, newdata = nd)
nd$est <- as.numeric(p$fit)
nd$lwr <- as.numeric(p$fit - 2 * p$se.fit)
nd$upr <- as.numeric(p$fit + 2 * p$se.fit)

gdat_day <- ggplot(dplot, aes(julian_day, change_in_crab, colour = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 50, dodge.width = 100), alpha = 0.6, pch = 21) +
  geom_line(aes(group = paste(block_id, treatment)), alpha = 0.2, lwd = 0.3) +
  scale_x_continuous(breaks = unique(dplot$julian_day), labels = unique(dplot$sample_date_group)) +
  theme(axis.text.x.bottom = element_text(angle = 90)) +
  scale_colour_brewer(palette = "Set2") +
  geom_pointrange(
    data = nd, mapping = aes(x = julian_day, y = est, ymin = lwr, ymax = upr, colour = treatment),
    position = position_dodge(width = 100)
  ) +
  labs(colour = "Treatment") +
  ylab("Change in crabs per 2m<sup>2</sup>") +
  xlab("Sample date") +
  theme(legend.position = c(0.85, 0.8)) +
  theme(axis.title.y = ggtext::element_markdown(), axis.text.x.bottom = element_text(size = 8))
gdat_day

g1 <- cowplot::plot_grid(gday, gmonth,
  ncol = 1,
  align = "v", labels = c("B", "C"), label_fontface = "plain"
)
g1
g <- cowplot::plot_grid(gdat_day, g1,
  ncol = 2, labels = "A",
  label_fontface = "plain"
) +
  theme(plot.margin = margin(2, 6, 1, 1))
g

ggsave("figs/crabs-change.pdf", width = 8.75, height = 4.25)
ggsave("figs/crabs-change.png", width = 8.75, height = 4.25)
