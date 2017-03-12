#!/usr/bin/Rscript --no-save --no-site-file

library(ggplot2)
library(dplyr)
library(magrittr)
library(ggthemes)
library(scales)
library(cowplot)


my_theme <- theme_cowplot() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.75),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90", size = 2),
        axis.title = element_text(face = "plain"),
        plot.title = element_text(vjust = 2))


roc_values <- NULL
load("cache/roc_values.Rda")
stopifnot(!is.null(roc_values))

roc_values %<>% filter(study != "unphased")

plt <- ggplot(roc_values, aes(false_pos, true_pos, lty = stat, color = stat)) + 
  geom_line() +
  scale_y_continuous(limit = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = c(0, .05, .25, .5, .75, 1),
                     labels = c("0", "0.05", "0.25", "0.5", "0.75", "1")) +
  scale_linetype(breaks = c("mcmf", "omega", "ihs"), 
                 labels = c("MCMF", "Omega", "iHS")) +
  scale_color_tableau(breaks = c("mcmf", "omega", "ihs"), 
                       labels = c("MCMF", "Omega", "iHS")) +
  labs(y = "Sensitivity", x = "False positive rate", 
       color = "Summary Statistic", lty = "Summary Statistic") + 
  my_theme
print(plt)

ggsave("figures/roc_curve.pdf", plt)

x11()
plot(roc_values$false_pos[roc_values$stat=="mcmf"],
     roc_values$true_pos[roc_values$stat=="mcmf"],t="l",
     xlab="False positive rate",ylab="Sensitivity",
     xlim=c(0,1),ylim=c(0,1),lwd=2)
lines(roc_values$false_pos[roc_values$stat=="omega"],
     roc_values$true_pos[roc_values$stat=="omega"],lty=2,lwd=2)
lines(roc_values$false_pos[roc_values$stat=="ihs"],
     roc_values$true_pos[roc_values$stat=="ihs"],lty=3,lwd=2)
abline(h=c(0,1))
abline(v=c(0,1))
legend("bottomright",bg="white",legend=c("MCSF","Omega","iHS"),lwd=2,lty=1:3)
## dev.copy2pdf(file="figures/roc_curve_bw.pdf")
