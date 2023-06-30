library(MixSIAR)
library(forcats)
library(car)
library(scales)
library(rgl)
library(magick)
library(R2jags)
library(rjags)
library(tidyverse)
library(dplyr)

###Boxplots
SRS_sumstats_gb<-read.csv('data/SRSMixout_gb.csv')
SRS_sumstats_gb$site<-fct_relevel(SRS_sumstats_gb$site, "SRS3","RB10", "SRS4","SRS6")
SRS_sumstats_gb$season <- gsub("wet19", "Wet 2019", SRS_sumstats_gb$season, ignore.case = TRUE)
SRS_sumstats_gb$season <- gsub("dry19", "Dry 2019", SRS_sumstats_gb$season, ignore.case = TRUE)
y_label_formatter <- function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

mixoutput_bxplt_gb_SRS<-ggplot(SRS_sumstats_gb,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+
  scale_fill_manual(values=c("saddlebrown",'limegreen'))+ 
  coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 12))+
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()

mixoutput_bxplt_gb_SRS

ggsave("figures/mixoutput_bxplt_gb_SRS.png", width = 10, height = 8, dpi = 300)




TS_sumstats_gb<-read.csv('data/TS_sumstats_gb.csv')
TS_sumstats_gb$site<-fct_relevel(TS_sumstats_gb$site, "TS3","TS7","TS9","TS10","TS11")

mixoutput_bxplt_gb_TS<-ggplot(TS_sumstats_gb,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+
  scale_fill_manual(values=c("saddlebrown",'limegreen'))+ 
  coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 12))+
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()


mixoutput_bxplt_gb_TS
ggsave("figures/mixoutput_bxplt_gb_TS.png", width = 10, height = 8, dpi = 300)


