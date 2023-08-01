library(ggplot2)
library(dplyr)
library(mgcv)
library(here)
theme_set(theme_light())
library(rstan)
rstan_options(auto_write = TRUE)
library(future)
# future::plan(multisession)

# Year to start erosion modelling (with respect to otter data; first year 1985)
# so, year start 8 means starting in (1985 + 8 - 1) = 1992 (erosion) and a 1 year lag (so 1991) for otters
YEAR_START <- 8

# Data --------------------------------------------------------------------

erosion <- readr::read_csv(here("data", "erosion_evd.csv"))
otter <- readr::read_csv(here("data", "sea_otter_population_data.csv"))
names(erosion) <- tolower(names(erosion))
names(otter) <- tolower(names(otter))
erosion$X <- NULL
erosion$Y <- NULL
erosion$id <- as.factor(as.character(erosion$id))

ggplot(otter, aes(year, sea_otter_no)) +
  geom_point() +
  geom_smooth(
    method = "gam", formula = y ~ s(x),
    method.args = list(family = nb(link = "log")), colour = "red"
  ) +
  ylab("Sea otter abundance") +
  xlab("Year")

# Explore -----------------------------------------------------------------

# - turn erosion back into width
# - model rate of growth of width from when it starts
# - model reduction in rate of growth per otter

fit_erosion <- gamm(erosion_m_y ~ s(year),
  random = list(id = ~1),
  family = gaussian(link = "identity"),
  data = erosion
)
nd_int <- data.frame(year = otter$year)
pred_erosion <- predict(fit_erosion$gam, newdata = nd_int)
se <- predict(fit_erosion$gam, newdata = nd_int, se.fit = TRUE)
png("figs/gomp-gam-erosion-year-choice.png", width = 4, height = 3.5, units = "in", res = 300)
par(mfrow = c(2, 1), cex = 0.7, mar = c(2, 4, 1, 1), oma = c(1, 0, 1, 1))
plot(nd_int$year, pred_erosion, type = "l", ylab = "Erosion (m/year)", xlab = "Year")
points(nd_int$year, pred_erosion)
abline(h = 0, lty = 2)
abline(v = YEAR_START + min(nd_int$year) - 1, lty = 2)

plot(otter$year, otter$sea_otter_no, type = "o", xlim = range(nd_int$year), ylab = "Sea otters", xlab = "")
abline(v = YEAR_START + min(nd_int$year) - 1, lty = 2)

dev.off()

n_t <- length(pred_erosion)
width_obs <- numeric(length = length(pred_erosion))
width_obs[1] <- 1
for (i in 2:length(width_obs)) {
  width_obs[i] <- width_obs[i - 1] + pred_erosion[i] * width_obs[i - 1]
}
plot(width_obs)
abline(v = YEAR_START)

# Lagged by a year:

# no lag:
# otters  erosion
# 1991    1991
# 1992    1992
# 1993    1993

# one year lag:
# otters  erosion
# 1991    1992
# 1992    1993
# 1993    1994

w <- width_obs[YEAR_START:length(width_obs)]
o <- otter[seq((YEAR_START - 1), nrow(otter) - 1), "sea_otter_no", drop = TRUE]
plot(o, w)
o[is.na(o)] <- mean(c(51, 67))
plot(w)

n_t <- length(o)
w_hat <- numeric(length = n_t)
w_hat[1] <- w[1]
r <- 0.2
for (i in seq(2, n_t)) {
  w_hat[i] <- w_hat[i - 1] * exp(r)
}

yrs <- nd_int$year[YEAR_START:length(nd_int$year)]
plot(yrs, w)
lines(yrs, w_hat)

plot(yrs, w)
lines(yrs, w_hat)
yrs_data <- intersect(nd_int$year, unique(erosion$year))
abline(v = yrs_data)
yrs_i <- which(yrs_data %in% yrs)

# Gompertz ----------------------------------------------------------------

# <https://en.wikipedia.org/wiki/Gompertz_function#Gompertz_curve>

gompertz <- function(x, limit, b) {
  limit * (1 - exp(-b * x))
}
plot(o, gompertz(o, 0.2, 0.05))
gompertz(0, 0.2, 0.05)
x <- seq(0, 130, length.out = 300)
plot(x, gompertz(x, 0.2, 0.05))

calc_width <- function(x) {
  width_obs <- numeric(length = length(x))
  width_obs[1] <- 1
  for (i in 2:length(width_obs)) {
    width_obs[i] <- width_obs[i - 1] + x[i] * width_obs[i - 1]
  }
  width_obs
}

width_df <- group_by(erosion, id) %>%
  arrange(id, year) %>%
  group_split() %>%
  purrr::map_dfr(function(.x) {
    data.frame(
      id = .x$id[1],
      year = .x$year, width = calc_width(.x$erosion_m_y)
    )
  })

width_df %>%
  filter(year >= min(nd_int$year)) %>%
  ggplot(aes(year, width, group = id)) +
  geom_line(alpha = 0.1) +
  coord_cartesian(xlim = c(1994, 2018), ylim = c(-10, 30))

# Do it in Stan once ------------------------------------------------------

stan_dat <- list(
  N = length(w),
  width = w,
  otters = o
)

fit <- stan(
  file = here("analysis/erosion-gompertz.stan"),
  data = stan_dat,
  chains = 4L,
  iter = 2000L,
  cores = 1L,
  control = list(adapt_delta = 0.95)
)
pars <- c("r", "a", "b", "sigma")
print(fit, pars = pars)

png("figs/gomp-mcmc-hex.png", width = 8, height = 8, units = "in", res = 200)
bayesplot::mcmc_pairs(fit, pars = pars, off_diag_fun = "hex")
dev.off()

e <- extract(fit)
yr_lu <- data.frame(year_i = seq_along(yrs), year = yrs)
obs <- data.frame(year = yrs, otters = o, width = w)
post <- extract(fit)
hist(post$r)
set.seed(1028)
post_mid <- e$mid
draws <- sample(seq_len(nrow(post$eta)), 100L)
eta <- post$eta[draws, ]
g <- reshape2::melt(eta) %>%
  rename(year_i = Var2, w_hat = value) %>%
  left_join(yr_lu) %>%
  ggplot(aes(year, w_hat, group = iterations)) +
  geom_line(alpha = 0.05) +
  geom_point(aes(year, width), data = obs, inherit.aes = FALSE) +
  geom_vline(xintercept = yrs_data) +
  coord_cartesian(xlim = range(yrs)) +
  geom_vline(xintercept = post_mid + min(yrs) - 1)
ggsave("figs/gomp-fitted-once.png", width = 6, height = 4)

g <- g + scale_y_log10()
ggsave("figs/gomp-fitted-once-log10.png", width = 6, height = 4)

hist(e$r)
hist(e$b)

plot(1, 1, type = "n", xlim = range(o), ylim = c(0, 0.5))
plot(o, gompertz(o, limit = e$r[1], b = e$b[1]))

fit_erosion <- gam(erosion_m_y ~ s(year),
  family = gaussian(link = "identity"),
  data = erosion
)
plot(fit_erosion)

# Simulate from the GAM posterior -----------------------------------------

# simulate directly from Gaussian approximate posterior...
br <- gam.mh(fit_erosion, thin = 2, ns = 2000, rw.scale = .15)
X <- predict(fit_erosion, newdata = nd_int, type = "lpmatrix")
p <- X %*% t(br$bs)

png("figs/gomp-gam-mh.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1, 2), cex = 0.7)
matplot(p, lty = 1, type = "l", col = "#00000010", xlab = "Year step", ylab = "Erosion (m/year)")
abline(v = YEAR_START)

n_t <- nrow(nd_int)
n_samps <- 200L
width_obs <- matrix(nrow = n_samps, ncol = n_t)
width_obs[, 1] <- 1 # FIXME?
set.seed(1028234)
draws <- sample(seq_len(ncol(p)), n_samps)
for (j in seq_along(draws)) {
  for (i in 2:ncol(width_obs)) {
    width_obs[j, i] <- width_obs[j, i - 1] + p[i, draws[j]] * width_obs[j, i - 1]
  }
}

matplot(t(width_obs), type = "l", log = "y", xlab = "Year step", ylab = "Width assuming starting value of 1", lty = 1, col = "#00000050")
abline(v = YEAR_START)
dev.off()

matplot(t(width_obs), type = "l", col = "black")
abline(v = YEAR_START)
w <- width_obs[, YEAR_START:ncol(width_obs)]
plot(apply(w, 2, mean))

# Fit to each sample ----------------------------------------------------

o <- otter[seq((YEAR_START - 1), nrow(otter) - 1), "sea_otter_no", drop = TRUE]
o[is.na(o)] <- mean(c(51, 67)) # FIXME integrate over in Stan? TODO NOTE
out <- purrr::map(1:nrow(w), function(i) {
# out <- furrr::future_map(1:nrow(w), function(i) {
  stan_dat <- list(
    N = ncol(w),
    width = w[i, ],
    otters = o
  )
  fit <- stan(
    file = here("analysis/erosion-gompertz.stan"),
    data = stan_dat,
    chains = 1L,
    iter = 500L,
    cores = 1L,
    # pars = c("f", "r", "sigma", "a"),
    control = list(adapt_delta = 0.95)
  )
  e <- extract(fit)
  e
})

png("figs/gomp-b-r-hist.png", width = 7, height = 4, units = "in", res = 200)
par(mfrow = c(1, 2), cex = 0.85)
b <- lapply(out, `[[`, "b") %>% unlist()
hist(b, main = "b posterior")

r <- lapply(out, `[[`, "r") %>% unlist()
hist(r, main = "r posterior")
dev.off()

zz <- purrr::map_dfr(out, function(.x) {
  reshape2::melt(.x$eta) %>%
    rename(year_i = Var2, w_hat = value)
})

med <- apply(w, 2, median)
lwr <- apply(w, 2, quantile, probs = 0.1)
upr <- apply(w, 2, quantile, probs = 0.9)
wq <- data.frame(med = med, lwr = lwr, upr = upr, year_i = seq_along(med))

zz %>%
  group_by(year_i) %>%
  summarise(
    lwr = quantile(w_hat, probs = 0.1),
    upr = quantile(w_hat, probs = 0.9),
    med = quantile(w_hat, probs = 0.5)
  ) %>%
  ggplot(aes(year_i, med)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  geom_linerange(aes(ymin = lwr, ymax = upr, x = year_i),
    colour = "blue", data = wq
  ) +
  geom_point(aes(y = med, x = year_i),
    colour = "blue", alpha = 0.9, data = wq
  ) +
  scale_y_log10() +
  ylab("Width (relative to starting value of 1)") +
  xlab("Year increment")
ggsave("figs/gomp-fitted-predicted-width.png", width = 5, height = 4)

aa <- purrr::map_dfr(out, function(.x) {
  reshape2::melt(.x$reff) %>%
    rename(year_i = Var2, reff = value)
})

aa %>%
  group_by(year_i) %>%
  summarise(
    lwr = quantile(reff, probs = 0.025),
    upr = quantile(reff, probs = 0.975),
    med = quantile(reff, probs = 0.5)
  ) %>%
  # ggplot(aes(year_i, exp(med) - 1)) +
  ggplot(aes(year_i, med)) +
  # geom_ribbon(aes(ymin = exp(lwr) - 1, ymax = exp(upr) - 1), alpha = 0.4) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  xlab("Year increment") +
  ylab(expression(Annual ~ rate ~ of ~ erosion ~ accounting ~ "for" ~ otters ~ (r[eff])))
ggsave("figs/gomp-reff-ts.png", width = 5, height = 4)

data.frame(b = e$b) %>%
  ggplot(aes(b)) +
  geom_density()

# ~2.5 % per otter at max initial point

ots <- seq(min(o), max(o), length.out = 200)

mm <- purrr::map_dfr(ots, function(.o) {
  data.frame(otters = .o, effect = (1 - exp(-b * .o)))
})

mm %>%
  group_by(otters) %>%
  summarise(
    lwr = quantile(effect, probs = 0.025),
    upr = quantile(effect, probs = 0.975),
    med = quantile(effect, probs = 0.5)
  ) %>%
  ggplot(aes(otters, med)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  ylab("Proportion reduction in base rate of erosion") +
  xlab("Otters")
ggsave("figs/gomp-gompertz-curve.png", width = 5, height = 4)

mm <- purrr::map_dfr(ots, function(.o) {
  data.frame(otters = .o, effect = r - r * (1 - exp(-b * .o)))
})

mm %>%
  group_by(otters) %>%
  summarise(
    lwr = quantile(effect, probs = 0.025),
    upr = quantile(effect, probs = 0.975),
    lwr2 = quantile(effect, probs = 0.25),
    upr2 = quantile(effect, probs = 0.75),
    med = quantile(effect, probs = 0.5)
  ) %>%
  ggplot(aes(otters, med)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_ribbon(aes(ymin = lwr2, ymax = upr2), alpha = 0.4) +
  ylab("Expected erosion rate") +
  xlab("Otters") +
  coord_cartesian(expand = FALSE) +
  ggsidekick::theme_sleek()

ggsave("figs/gomp-expected-erosion-otters.png", width = 5, height = 4)

future::plan(sequential)
saveRDS(out, file = "analysis/gompertz-mcmc.rds")
