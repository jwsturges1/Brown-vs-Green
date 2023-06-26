getwd()
setwd("/Users/mack/Desktop/RESEARCH/FCE FOOD WEBS/FOOD WEBS R")
dat = read.csv("fce_tn and tp_all.csv")
datMEAN = read.csv("np_ind por_fceasm.csv")
datLONG = read.csv("np_ind por_fceasm_long.csv")


plot <- ggplot(datLONG, aes(SITE, NP, group = POR)) +
  geom_point() +
  geom_line( = ) +
  labs(x = "Sampling Site", 
       y = "N:P Mass Ratio (ug)", 
       title = "Spatiotemporal Variability in N:P Ratios: Shark River Slough")

