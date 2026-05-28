#bsp-L3
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("/Users/zhuchendi/Desktop/TB234/BEAST/")
df <- read.csv("L3iall.csv", header = TRUE, fill = TRUE)

clade_colors <- c(
  "L3.1.1.1" = "#F5CE5A",
  "L3.1.1.2" = "#8B78BC",
  "L3.1.1.3" = "#AE5D85",
  "L3.1.1.4" = "#74AE5A",
  "L3.1.1.5" = "#E99C4C",
  "L3.1.1.6" = "#CF0909",
  "other" = "#A9A9A9"
)

p <- ggplot(df, aes(x = year, y = log, color = sublineage)) +
  geom_line(size = 1.5) +
  scale_x_continuous(
    name = "Year",
    limits = c(min(df$year), 2010),
    breaks = seq(floor(min(df$year)), ceiling(max(df$year)), by = 100)
  ) +
  scale_color_manual(values = clade_colors) +
  theme_bw() +
  ylab(expression(log(Ne))) +
  theme(
    axis.text = element_text(size = 25, color="black"),
    axis.title = element_text(size = 25, color="black"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

p
ggsave("BSP_all_sublineages.pdf", p, width = 12, height = 8)
