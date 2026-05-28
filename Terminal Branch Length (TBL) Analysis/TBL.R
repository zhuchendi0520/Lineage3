setwd("/Users/zhuchendi/Desktop/TB234/tree/")
# TBL
# raincloud
library(ggplot2)
library(ggdist)
library(dplyr)
library(gghalves)
data <- read.csv("all.csv")
data$Tip.Name <- factor(data$Tip.Name, levels = c(
  "L2.1","L2.2","L2.3",
  "L3.1",
  "L4.1","L4.2","L4.3","L4.4","L4.5","L4.6","L4.7","L4.8","L4.9","L4.11"
))

plot <- ggplot(data, aes(x = Tip.Name, y = TBL, colour = Lineage, fill = Lineage)) +
  stat_halfeye(
    adjust = 25, 
    width = 0.65, 
    justification = -0.12, 
    .width = 0, 
    point_colour = NA, 
    alpha = 0.6
  ) +
  
  geom_boxplot(
    width = 0.15,
    fill ="white",
    outlier.shape = NA, 
    alpha = 0.6
  ) +
  geom_point(
    aes(x = as.numeric(factor(Tip.Name)) - 0.2, y = TBL, color = Lineage),
    position = position_jitter(width = 0.05),
    size = 0.2, shape = 20,alpha = 0.3
  ) +
  
  scale_color_manual(values = c("#3d85c6","#674ea7","#cc0000")) +
  scale_fill_manual(values = c("#3d85c6","#674ea7","#cc0000")) + 
  
  labs(x = NULL, y = "TBL") +
  theme_classic() +
  theme(
    axis.text.x = element_text( size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks = element_line(size = 0.5, colour = "black"),
    axis.ticks.length = unit(0.3, "cm")
  ) +coord_flip() +
  ylim(0,200)

plot


png(filename = "raincloud_subLineage.png",width = 8,height = 4,units = "in",bg = "white",res = 300)
plot
dev.off()



# violin-plot
setwd("/Users/zhuchendi/Desktop/TB234/tree/")
library(ggplot2)
library(ggdist)
library(dplyr)
library(gghalves)
data <- read.csv("Lineage4.9_TBL.csv")

plot <- ggplot(data, aes(x = Type, y = TBL,fill = Type)) +
  geom_violin(adjust =10)+
  geom_boxplot(
    width = 0.15,   
    fill ="white",
    outlier.shape = NA, 
    alpha = 0.9
  ) +
  
  scale_color_manual(values = c("#3d85c6","#674ea7","#cc0000")) +
  scale_fill_manual(values = c("#3d85c6","#674ea7","#cc0000")) + 
  
  labs(x = NULL, y = "TBL") +
  theme_classic() +
  theme(
    axis.text.x = element_text( size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks = element_line(size = 0.5, colour = "black"),
    axis.ticks.length = unit(0.3, "cm")
  ) 
plot
ttest_result <- t.test(TBL ~ Type, data = data)
ttest_result
png(filename = "raincloud.png",width = 8,height = 4,units = "in",bg = "white",res = 300)
plot
dev.off()



#heatmap
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

anova_result <- aov(TBL ~ Tip.Name, data = data)
summary(anova_result)
res <- as.data.frame(TukeyHSD(anova_result)$Tip.Name)
res$Comparison <- rownames(res)

order_lineage <- c("L2.1","L2.2","L2.3",
                   "L4.1","L4.2","L4.3","L4.4",
                   "L4.5","L4.6","L4.7","L4.8",
                   "L4.9","L4.11","L3.1")
df <- res %>%
  mutate(
    g1_raw = word(Comparison, 1, sep = "-"),
    g2_raw = word(Comparison, 2, sep = "-"),
    idx1 = match(g1_raw, order_lineage),
    idx2 = match(g2_raw, order_lineage),
    
    g1 = ifelse(idx1 < idx2, g1_raw, g2_raw),
    g2 = ifelse(idx1 < idx2, g2_raw, g1_raw),
    diff = ifelse(idx1 < idx2, diff, -diff),
    `p adj` = `p adj`
  ) %>%
  filter(!is.na(g1) & !is.na(g2)) %>%
  distinct(g1, g2, .keep_all = TRUE) %>%
  mutate(
    g1 = factor(g1, levels = order_lineage),
    g2 = factor(g2, levels = order_lineage),
    significance = case_when(
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01  ~ "**",
      `p adj` < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )


diag_df <- data.frame(
  g1 = order_lineage,
  g2 = order_lineage,
  diff = 0,
  significance = ""
)

plot <-ggplot(df, aes(x = g1, y = g2, fill = diff)) +
  geom_tile(color = "black", size = 0.5) +  
  geom_text(aes(label = significance), size = 4, color = "black") +
  scale_fill_gradient2(low = "#4575b4",
                       mid = "white",
                       high = "#d73027",
                       midpoint = 0,
                       name = "") +
  scale_y_discrete(limits = rev(order_lineage)) +
  geom_tile(data = diag_df, aes(x = g1, y = g2),
            fill = "white", color = "black", size = 0.5) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  scale_x_discrete(limits = order_lineage) +    
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
plot
png(filename = "heatmap_subLineage.png",width = 6,height =6,units = "in",bg = "white",res = 300)
plot
dev.off()