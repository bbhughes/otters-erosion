library(ggplot2)
library(dplyr)
dir.create("figs", showWarnings = FALSE)

source("theme_sleek.R")
theme_set(theme_sleek())
sig_size <- 5

#### TRANSECTS, No assumptions of when sea otters occupied tidal creeks prior to 2013.####
crab_marsh15 <- read.csv("otter_marsh_transects.csv")
str(crab_marsh15)
crab_marsh15$Year <- as.factor(as.integer(crab_marsh15$Year))
crab_marsh15$Otters_ha <- as.numeric(as.character(crab_marsh15$usgs_oph_2015))
crab_marsh15$bg_biomass_2015_kg_msq <- as.numeric(as.character(crab_marsh15$bg_biomass_2015_kg_msq))
crab_marsh15$bank_retreat_2013_2015_m_yr <- as.numeric(as.character(crab_marsh15$bank_retreat_2013_2015_m_yr))
str(crab_marsh15)

# GLM model using log link
crab_marsh15_glm1 <- glm(Density ~ usgs_oph_2015, data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))
crab_marsh15_glm2 <- glm(Density ~ log(usgs_oph_2015), data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))
AIC(crab_marsh15_glm1, crab_marsh15_glm2)
summary(crab_marsh15_glm1)
# P = 0.0027

min_x <- min(crab_marsh15$usgs_oph_2015)
max_x <- max(crab_marsh15$usgs_oph_2015)
min_y <- 0
max_y <- 5
crab_otter15_SI_plot <- ggplot(data = crab_marsh15, aes(x = usgs_oph_2015, y = Density)) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5, position = position_jitter(w = 0, h = 0.02)) +
  geom_smooth(
    method = "glm", formula = y ~ x, se = TRUE, colour = "black",
    alpha = 0.3, method.args = list(family = Gamma(link = "log"))
  ) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
  scale_y_continuous(limits = c(min_y, max_y), expand = expansion(mult = c(0, .04))) +
  annotate("text", x = 0.05, y = 4.5, label = "*", col = "grey30", size = sig_size) +
  labs(y = "Crab density per trap", x = "Otters per ha")
crab_otter15_SI_plot

# Comparison of otter abundance and Pachygrapsus consumed per day
# GLM model using log tranformation
pachy_otter_glm <- glm(usgs_pachy_eaten_ha_day ~ usgs_oph_2015, data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))
summary(pachy_otter_glm)

pachy_otter_glm2 <- glm(usgs_pachy_eaten_ha_day ~ log(usgs_oph_2015), data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))

AIC(pachy_otter_glm, pachy_otter_glm2)

summary(pachy_otter_glm2)
# P=0.0043

# Plot the data
pachy_otter_plot <- ggplot(data = crab_marsh15, aes(x = usgs_oph_2015, y = (usgs_pachy_eaten_ha_day))) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5, position = position_jitter(w = 0, h = 0)) +
  geom_smooth(method = "glm", formula = y ~ log(x), se = TRUE, colour = "black", alpha = 0.3, method.args = list(family = Gamma(link = "log"), na.action = na.omit)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04)), limits = c(0, NA)) +
  annotate("text", x = 0.005, y = 17, label = "**", col = "grey30", size = sig_size) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) +
  labs(y = expression(atop("Crabs consumed", "per day per ha")), x = "Otters per ha")
pachy_otter_plot

# Comparison of bank retreat and Pachygrapsus consumed per day,
# LM model

pachy_retreat_lm <- lm(bank_retreat_2013_2015_m_yr ~ usgs_pachy_eaten_ha_day, data = crab_marsh15, na.action = na.omit, family = gaussian(link = "identity"))
summary(pachy_retreat_lm)

# P= 0.089

r_pachy_retreat_glm <- DHARMa::simulateResiduals(pachy_retreat_glm)
plot(r_pachy_retreat_glm)

# P= 0.0091

# Plot the data
min_y <- min(crab_marsh15$bank_retreat_2013_2015_m_yr)
max_y <- max(crab_marsh15$bank_retreat_2013_2015_m_yr)
min_x <- min(crab_marsh15$usgs_pachy_eaten_ha_day)
max_x <- max(crab_marsh15$usgs_pachy_eaten_ha_day)
pachy_retreat_plot <- ggplot(data = crab_marsh15, aes(x = (usgs_pachy_eaten_ha_day), y = bank_retreat_2013_2015_m_yr)) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5, position = position_jitter(w = 0.02, h = 0)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", alpha = 0.3, method.args = list(family = gaussian(link = "identity"))) +
  scale_x_continuous(limits = c(NA, NA), expand = expansion(mult = c(0, .03))) +
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0, .04))) +
  annotate("text", x = 8, y = 0.5, label = "P = 0.089", col = "grey30", size = 3) +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0)) +
  labs(y = expression(atop("Bank erosion", "(m per year)")), x = "Crabs consumed per ha per day")
pachy_retreat_plot

# Comparison of marsh biomass and Pachygrapsus consumed per day,
# GLM model using log tranformation
pachy_bg_glm <- glm(bg_biomass_2015_kg_msq ~ usgs_pachy_eaten_ha_day, data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))
summary(pachy_bg_glm)
# P = 0.174

pachy_bg_glm2 <- glm(bg_biomass_2015_kg_msq ~ log(usgs_pachy_eaten_ha_day), data = crab_marsh15, na.action = na.omit, family = Gamma(link = "log"))
AIC(pachy_bg_glm, pachy_bg_glm2)
summary(pachy_bg_glm2)
# P = 0.0329

# Plot the data
min_y <- min(crab_marsh15$bg_biomass_2015_kg_msq)
max_y <- max(crab_marsh15$bg_biomass_2015_kg_msq)
min_x <- min(crab_marsh15$usgs_pachy_eaten_ha_day)
max_x <- max(crab_marsh15$usgs_pachy_eaten_ha_day)
pachy_bg_plot <- ggplot(data = crab_marsh15, aes(x = (usgs_pachy_eaten_ha_day), y = bg_biomass_2015_kg_msq)) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5, position = position_jitter(w = 0, h = 0)) +
  geom_smooth(method = "glm", formula = y ~ log(x), se = TRUE, colour = "black", alpha = 0.3, method.args = list(family = Gamma(link = "log"))) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  annotate("text", x = 0.5, y = 3, label = "*", col = "grey30", size = sig_size) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0)) +
  xlab("Crabs consumed per ha per day") +
  ylab(expression(atop("Belowground mass ", (kg(dw) %.% m^-2))))
pachy_bg_plot

# Save Spatial Analysis Figure
# pdf("Spatial_Figure.pdf", width=7, height=15)
# gridExtra::grid.arrange(pachy_otter_plot,pachy_bg_plot,pachy_retreat_plot, ncol = 1, nrow = 3)
# dev.off()

half_line <- 11 / 2
m <- theme(plot.margin = margin(t = half_line - 0.5, r = half_line, b = half_line - 4, l = half_line - 1))
cowplot::plot_grid(
  pachy_otter_plot + m,
  pachy_bg_plot + m,
  pachy_retreat_plot + m,
  ncol = 1, nrow = 3, align = "v"
) + theme(plot.margin = margin(0, 3, 0, 0))
ggsave("figs/Spatial_Figure.pdf", width = 3.25, height = 6.5)

