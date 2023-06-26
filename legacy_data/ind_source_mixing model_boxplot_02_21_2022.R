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

getwd()

setwd("~/Desktop/RESEARCH/FCE FOOD WEBS/FOOD WEBS R")

# ONE LARGE BOXPLOT BELOW

TS_sumstats_gb<-read.csv('TS_sumstats_all.csv')
TS_sumstats_gb$site<-fct_relevel(TS_sumstats_gb$site, "TS3","TS7","TS9","TS10","TS11")

mixoutput_bxplt_gb_TS<-ggplot(TS_sumstats_gb,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw() + coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c("springgreen2","olivedrab", "tan4", "tan3"))+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")

mixoutput_bxplt_gb_TS

# NEED TO GO TRHOUGH AND SUBSET FOR EACH INDIVIDUAL SITE AND MAKE BOXPLOTS

TS3 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS3")

# nothing from TS3 so we dont need to worry about it!

TS7 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS7")

TS9 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS9")

TS10 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS10")

TS11 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS11")

# alright, now the data is all there - we needed to get through and subset each site so we don't have
# a giant, messy boxplot -> need to go through code as above and do it for each site!

# lemme know if you need any help!

# TS7

mixoutput_bxplt_gb_T7<-ggplot(TS7,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw() + coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c("springgreen2","olivedrab", "tan4", "tan3"))
scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")

mixoutput_bxplt_gb_T7


# TS9

mixoutput_bxplt_gb_T9<-ggplot(TS9,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw() + coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c("springgreen2","olivedrab", "tan4", "tan3"))
scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")

mixoutput_bxplt_gb_T9

# TS10

mixoutput_bxplt_gb_T10<-ggplot(TS10,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw() + coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c("springgreen2","olivedrab", "tan4", "tan3"))
scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")

mixoutput_bxplt_gb_T10

#TS11

mixoutput_bxplt_gb_T11<-ggplot(TS11,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw() + coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c("springgreen2","olivedrab", "tan4", "tan3"))
scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")

mixoutput_bxplt_gb_T11

################## UPDATED PLOTS WITH CODE FROM RYAN

TS_sumstats_gb<-read.csv('TS_sumstats_all.csv')
TS_sumstats_gb$site<-fct_relevel(TS_sumstats_gb$site, "TS3","TS7","TS9","TS10","TS11")

# NEED TO GO TRHOUGH AND SUBSET FOR EACH INDIVIDUAL SITE AND MAKE BOXPLOTS

TS3 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS3")

# nothing from TS3 so we dont need to worry about it!

TS7 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS7")

TS9 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS9")

TS10 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS10")

TS11 = TS_sumstats_gb %>% filter(TS_sumstats_gb$site == "TS11")

# START WITH PLOTTING EVERYTHING

bxplt_TSALL<-ggplot(TS_sumstats_gb,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+scale_fill_manual(values=c("limegreen","seagreen2","saddlebrown","wheat3"))+ coord_cartesian(ylim = c( 0,1))+facet_grid(site~.)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")+coord_flip()

bxplt_TSALL

# PLOT SITE TS7

bxplt_TS7<-ggplot(TS7,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+scale_fill_manual(values=c("limegreen","seagreen2","saddlebrown","wheat3"))+ coord_cartesian(ylim = c( 0,1))+facet_grid(site~.)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")+coord_flip()

bxplt_TS7

# PLOT SITE TS9

bxplt_TS9<-ggplot(TS9,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+scale_fill_manual(values=c("limegreen","seagreen2","saddlebrown","wheat3"))+ coord_cartesian(ylim = c( 0,1))+facet_grid(site~.)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")+coord_flip()

bxplt_TS9

# PLOT SITE TS10

bxplt_TS10<-ggplot(TS10,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+scale_fill_manual(values=c("limegreen","seagreen2","saddlebrown","wheat3"))+ coord_cartesian(ylim = c( 0,1))+facet_grid(site~.)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")+coord_flip()

bxplt_TS10

# PLOT SITE TS11

bxplt_TS11<-ggplot(TS11,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
  theme_bw()+scale_fill_manual(values=c("limegreen","seagreen2","saddlebrown","wheat3"))+ coord_cartesian(ylim = c( 0,1))+facet_grid(site~.)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0.01, 0))+labs(y="Dietary contribution")+coord_flip()

bxplt_TS11
