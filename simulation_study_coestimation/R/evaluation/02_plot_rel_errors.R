#!/usr/bin/Rscript --no-site-file

library(ggplot2)
library(ggthemes)
library(dplyr)
library(scales)
library(cowplot)

load("cache/sim-data-df.Rda")
options(stringsAsFactors = FALSE)

my_theme <- theme_cowplot() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.75, linetype='dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey90', size = 2, linetype='solid'),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90", size = 2, linetype='solid'),
        axis.title = element_text(face = "plain"),
        plot.title = element_text(vjust = 2))

facet_labels <- list(
  model = list(
    IWM = "IWM",
    Bottleneck = "Bottleneck",
    Growth = "Growth"
  ),
  par = list(
    m = "m",
    tau = expression(tau),
    theta = expression(theta),
    rho = expression(rho),
    q = "q",
    s = "s"
  )
)

facet_labeller <- function(variable, value) {
  if (any(value == 1 || value == "1")) stop("Invalid values")
  facet_labels[[variable]][value]
}


# ----------------------------------------------------------------
# Overview plots for Appendix
# ----------------------------------------------------------------

plt_dm_neutr <- ggplot(sim_data %>% filter(par != "s" & selection == "Neutral"),
              aes(s, estimate / true_value, fill = s)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  facet_grid("model ~ par", labeller = facet_labeller) +
  xlab("True Fraction of Selected Loci") +
  scale_y_continuous("Estimate : True Value",
                     breaks=c(1/4, 0.5, 1, 2, 4),
                     label=c("1:4", "1:2", "1:1", "2:1", "4:1"),
                     limits = c(.3, 2.2),
                     trans=log2_trans()) +
  guides(fill=FALSE) +
  my_theme + scale_fill_tableau() +
  ggtitle("Neutral Analysis")

ggsave(plt_dm_neutr, file = "plots/demography_neutral_estimates.pdf", width = 10, height = 12)


plt_dm_sel <- ggplot(sim_data %>% filter(par != "s" & selection == "Selection"),
                 aes(s, estimate / true_value, fill = s)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  facet_grid("model ~ par", labeller = facet_labeller) +
  xlab("True Fraction of Selected Loci") +
  scale_y_continuous("Estimate : True Value",
                     breaks=c(1/4, 0.5, 1, 2, 4),
                     limits = c(.3, 2.2),
                     label=c("1:4", "1:2", "1:1", "2:1", "4:1"),
                     trans=log2_trans()) +
  guides(fill=FALSE) +
  my_theme + scale_fill_tableau() +
  ggtitle("Co-estimation of Demography and Selection")

ggsave(plt_dm_sel, file = "plots/demography_coestimation_estimates.pdf", width = 10, height = 12)


plt_s <- ggplot(subset(sim_data, par == "s"),
                aes(as.factor(true_value), estimate,
                    fill = as.factor(true_value))) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  facet_grid(". ~ model", labeller = facet_labeller) +
  xlab("True Fraction of Selected Loci s") +
  scale_y_continuous("Estimated Fraction s",
                     breaks = c(0.05, 0.1, 0.25, 0.5, 0.8)) +
  guides(fill=FALSE) +
  my_theme + scale_fill_tableau()

ggsave(plt_s, file = "plots/s_estimates.pdf", width = 10, height = 6)



# ----------------------------------------------------------------
# Bottleneck Model for manuscript
# ----------------------------------------------------------------
plt_dm_neutr <- ggplot(sim_data %>% filter(par != "s" & model == "Bottleneck" & selection == "Neutral"),
                 aes(s, estimate / true_value, fill = s)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  #geom_point(size = 3) +
  facet_grid(". ~ par", labeller = facet_labeller) +
  xlab("True Fraction of Selected Loci s") +
  scale_y_continuous("Estimate : True Value",
                     breaks = c(1/4, 0.5, 1, 2, 4),
                     label = c("1:4", "1:2", "1:1", "2:1", "4:1"),
                     limits = c(.3, 2.2),
                     trans = log2_trans()) +
  guides(fill = FALSE, col = FALSE) +
  my_theme + scale_fill_tableau() +
  ggtitle("Neutral Analysis")


plt_dm_sel <- ggplot(sim_data %>% filter(par != "s" & model == "Bottleneck" & selection == "Selection"),
                       aes(s, estimate / true_value, fill = s)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  facet_grid(". ~ par", labeller = facet_labeller) +
  xlab("True Fraction of Selected Loci s") +
  scale_y_continuous("Estimate : True Value",
                     breaks=c(1/4, 0.5, 1, 2, 4),
                     label=c("1:4", "1:2", "1:1", "2:1", "4:1"),
                     limits = c(.3, 2.2),
                     trans=log2_trans()) +
  guides(fill = FALSE, col = FALSE) +
  my_theme + scale_fill_tableau() +
  ggtitle("Co-estimation of Demography and Selection")


plt_s <- ggplot(subset(sim_data, par == "s" & model == "Bottleneck"),
                aes(as.factor(true_value), estimate,
                    fill = as.factor(true_value))) +
  geom_boxplot(position = "dodge", outlier.size = 0) +
  xlab("True Fraction of Selected Loci") +
  scale_y_continuous("Estimated Fraction s",
                     breaks = c(0.05, 0.1, 0.25, 0.5, 0.8)) +
  guides(fill = FALSE, col = FALSE) +
  my_theme + scale_fill_tableau()
ggsave(plt_s, file = "plots/bn_s_estimates.pdf", width = 10, height = 12)


plt_both <- plot_grid(plt_dm_neutr, plt_dm_sel, labels = c("A", "B"),
                      rel_heights = c(1, 1), nrow = 2)
ggsave(plt_both, file = "plots/bn_estimates.pdf", width = 10, height = 12)

