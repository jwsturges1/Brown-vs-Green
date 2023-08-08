# Plot that generates CN and CS biplots for all 9 sites monitored



# read in libraries
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

# create list of sites as i vector
sites = c("SRS3", "SRS4","SRS6", "RB10", "TS3", "TS7", "TS9", "TS10", "TS11")

# produces two biplots per site (one CN and one CS)
for(i in 1:length(sites)) {
  
  fish = read_csv(paste0('data/', sites[i], 'mix.csv'))
  
  sources = read.csv(paste0('data/sources', sites[i], '.csv')) %>% 
    mutate(Meand13C = Meand13C + 1.3,
           Meand15N = Meand15N + 3.3,
           Meand34S = Meand34S + .5)
  
  
  wcn = ggplot(data = sources, aes(Meand13C, Meand15N))+
    geom_point(data = sources, size = 3, pch=c(20))+ 
    geom_errorbar(data = sources, aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0) + 
    geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
    ylab(expression(paste(delta^{15}, "N (\u2030)")))+
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    labs(title = paste("CN ", sites[i])) +  # Add title
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    geom_point(data = fish, aes(x = d13C, y = d15N,color = common_name), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(0,14))+
    theme( legend.title = element_blank(),
           legend.text=element_text(size=12))#,legend.position=c(.85,.15))
  
  ggsave(paste0('figures/biplots/',sites[i], 'CN.pdf'), units="in", width=10, height=6)
  
  # C and S
  wcs = ggplot(data = sources, aes(Meand13C, Meand34S))+
    geom_point(data = fish, aes(x = d13C, y = d34S,color = common_name), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    geom_point(data = sources, size = 3, pch=c(20))+ 
    geom_errorbar(data = sources, aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
    geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
    ylab(expression(paste(delta^{34}, "S (\u2030)")))+
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    labs(title = paste("CS ", sites[i])) +  # Add title
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(-16.5, 24))+
    theme(legend.title = element_blank())#, legend.position=c(.85,.85))
  ggsave(paste0('figures/biplots/',sites[i], 'CS.pdf'), units="in", width=10, height=6)
  
}



# creates the same isotope biplots but colors points by species to look at wet dry trends
for(i in 1:length(sites)) {
  
  fish = read_csv(paste0('data/', sites[i], 'mix.csv'))
  
  sources = read.csv(paste0('data/sources', sites[i], '.csv')) %>% 
    mutate(Meand13C = Meand13C + 1.3,
           Meand15N = Meand15N + 3.3,
           Meand34S = Meand34S + .5)
  
  
  wcn_season = ggplot(data = sources, aes(Meand13C, Meand15N))+
    geom_point(data = sources, size = 3, pch=c(20))+ 
    geom_errorbar(data = sources, aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0) + 
    geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
    ylab(expression(paste(delta^{15}, "N (\u2030)")))+
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    labs(title = paste("CN ", sites[i])) +  # Add title
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    geom_point(data = fish, aes(x = d13C, y = d15N,color = hydroseason), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(0,14))+
    theme( legend.title = element_blank(),
           legend.text=element_text(size=12))#,legend.position=c(.85,.15))
  
  ggsave(paste0('figures/biplots/',sites[i], 'CN_WD.pdf'), units="in", width=10, height=6)
  
  # C and S
  wcs_season = ggplot(data = sources, aes(Meand13C, Meand34S))+
    geom_point(data = fish, aes(x = d13C, y = d34S,color = hydroseason), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    geom_point(data = sources, size = 3, pch=c(20))+ 
    geom_errorbar(data = sources, aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
    geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
    ylab(expression(paste(delta^{34}, "S (\u2030)")))+
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    labs(title = paste("CS ", sites[i])) +  # Add title
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(-16.5, 24))+
    theme(legend.title = element_blank())#, legend.position=c(.85,.85))
  ggsave(paste0('figures/biplots/',sites[i], 'CS_WD.pdf'), units="in", width=10, height=6)
  
}
