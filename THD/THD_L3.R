library(thd)
library(ape)
library(readxl)
library(data.table)
library(stringr)
library(lmerTest)
library(beeswarm)

# Load FASTA then build SNP distance matrix
setwd("/Users/zhuchendi/Desktop/TB234/dr/THD/L3i6")
dna <- read.dna("L3i6.fa", "fasta")
dna_len <- 1978
dna_names <- sort(labels(dna))

# Hamming distance matrix
H <- dist.dna(dna, "raw", pairwise.deletion = T, as.matrix = T) * dna_len

# Load metadata
read.csv("L3i6.csv") ->d
d <- as.data.table(d)
H <- H[ d$Names, d$Names ]

# truncation limit p
dtgeom <- function(x, alpha, p) {
  k <- (1 - alpha) / (1 - alpha^(p+1)) * alpha^x
  k[!is.finite(k)] <- 1/p
  return(k)
}
# See https://www.nature.com/articles/srep45326 for mathematical details
gkdew <- function(h, alpha = 0.5, p = max(h), weights, from, to, skipself = T) {
  if(missing(from) & missing(to)) {
    # Pairwise densities
    n <- nrow(h)
    k <- NULL 
    if(skipself) 
      k <- sapply(1:n, function(i) {
        weighted.mean(dtgeom(h[,i][-i], alpha, p), weights[-i])
      })
    else
      k <- sapply(1:n, function(i) {
        weighted.mean(dtgeom(h[,i], alpha, p), weights[-i])
      })
    return(k)		
  } else {
    # Densities from one group to another
    # Chekc whether from = to and skip diagonal
    if(missing(from)) stop("Argument 'to' cannot be provided alone")
    if(missing(to)) to <- from
    from[is.na(from)] <- FALSE
    to[is.na(to)]     <- FALSE
    skipdiag <- all(from == to)
    hh <- h[to, from] # Scan columns rather than rows for efficiency
    if(sum(to) == 1)   hh <- matrix(hh, 1, sum(from))
    if(sum(from) == 1) hh <- matrix(hh, 1, sum(to))
    n <- ncol(hh)
    kdefunc <- if(skipdiag) {
      function(i, hh, alpha, p) {weighted.mean(dtgeom(hh[,i][-i], alpha, p), weights[-i])}
    } else {
      function(i, hh, alpha, p) {weighted.mean(dtgeom(hh[,i], alpha, p), weights[-i])}
    }
    k <- sapply(1:n, kdefunc, hh = hh, alpha = alpha, p = p)
    return(k)
  }
}

# THD parameters
# Effective genome size (H37RV size with 10% masked)
effsize <- 3969000
# Mutation rate per site per year
mu <- 4.6e-8
# 10y timescale, recent epidemic
timescale <- 10
bandwidth <- thd.bandwidth(timescale, effsize, mu, 1/2)

# Unweighted THD
thd <- thd::thd(H, timescale, effsize, mu)
d[, thd := thd]

# Use display form of THD, x100
d[, thd100 := thd * 300]

# Figure 5 panel a

# Helper function to control color transparency
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

d <- d[!is.na(Names)]
write.csv(d,"L3i6.csv")

table(d$group)

read.csv("all_THD.csv")->d
d<-read.csv("L3i6.csv")
library(data.table)
d <- as.data.table(d)

d$group <- interaction(d$type, d$possible_compensation, drop = TRUE)
order_groups <- c(
  "Sensitive.0",
  "other_DR.0","MDR.0","pre-XDR.0",
  "other_DR.>=1","MDR.>=1","pre-XDR.>=1"
)

d$group <- factor(d$group, levels = order_groups)
d <- d[!is.na(d$group), ]


set.seed(123)

d0 <- d[possible_compensation == 0]
d1 <- d[!(possible_compensation == 0)]
d0 <- d0[thd100 >= 0.00001]
d0_sub <- d0[sample(.N, size = floor(.N * 0.4))]


d_plot <- rbind(d0_sub, d1)


at_pos <- c(1,2,3,4, 6,7,8)

comp_colors <- c("darkturquoise", "lightcoral")
col_vector <- c(rep(add.alpha(comp_colors[1], .25), 4),
                rep(add.alpha(comp_colors[2], .25), 3))

png(
  filename = "thd_boxplot.png",
  width = 7.5, height = 6,
  units = "in", res = 300, bg = "white"
)

par(mar=c(6,6,2,2))
boxplot(thd100 ~ group, data=d,
        at = at_pos,
        col = col_vector,
        xaxt = "n", outline = FALSE,
        ylab = "THD success index x300",
        xlab = "", las = 2)

beeswarm(thd100 ~ group, data=d_plot,
         at = at_pos,
         spacing=0.8, 
         pch=19, cex=0.6,
         add=TRUE, method="center",
         col=c(rep(add.alpha(comp_colors[1], .75), 4),
               rep(add.alpha(comp_colors[2], .75), 3)))

axis(1, at=at_pos,
     labels=c("S","other_DR","MDR","Pre-XDR",
              "other_DR","MDR","Pre-XDR"),
     las=1, cex.axis=0.9)


legend("topleft", bty="n",
       title="Possible compensatory mutations",
       legend=c("0",">=1"), fill=comp_colors)

dev.off()

#DR vs S
comp_colors <- c("darkturquoise", "lightcoral")
at_pos <- c(1,2)
png(
  filename = "thd_dr_boxplot.png",
  width = 4, height = 6,
  units = "in", res = 300, bg = "white"
)

par(mar=c(6,6,2,2))
boxplot(thd100 ~ type1, data=d,
        at = at_pos,
        col = add.alpha(comp_colors, alpha=0.4),
        xaxt = "n", outline = FALSE,
        ylab = "THD success index x300",
        xlab = "", las = 2)

axis(1, at=at_pos,
     labels=c("S","DR"),
     las=1)

dev.off()
t.test(thd100 ~ type1, data = d)


library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
df <- fread("
Sublineage	other_DR	MDR	pre_XDR	XDR	S
L3.1.1.1	94	114	16	0	790
L3.1.1.2	152	126	71	3	618
L3.1.1.3	140	86	84	1	618
L3.1.1.4	178	95	93	2	894
L3.1.1.5	32	45	9	0	326
L3.1.1.6	98	90	78	0	320
")


# long format
df_long <- melt(df,
                id.vars = "Sublineage",
                variable.name = "Type",
                value.name = "Count")

setDT(df_long) 
df_long[, Percent := Count / sum(Count), by = Sublineage]
my_colors <- c(
  "S"        = "#C1A872",
  "other_DR" = "#558295",
  "MDR"      = "#75638D",
  "pre_XDR"  = "#CD6B6C",
  "XDR"      = "#6DAA65"
)
df_long$Type <- factor(df_long$Type, levels = c("XDR","pre_XDR","MDR", "other_DR", "S" ))

p <- ggplot(df_long, aes(x = Sublineage, y = Percent, fill = Type)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = my_colors,
                    breaks = names(my_colors)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    panel.border = element_rect(color = "black", size = 1, fill = NA)
  ) +
  
  labs(
    x = "",
    y = "Percentage of strains",
    fill = "Resistance type"
  )

png(
  filename = "perentage_dr.png",
  width = 8, height = 6,
  units = "in", res = 300, bg = "white"
)
p
dev.off()
