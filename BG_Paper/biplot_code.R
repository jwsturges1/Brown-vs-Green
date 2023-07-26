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
fish = read_csv('data/TS11mix.csv')

sources = read.csv('data/sourcesTS11.csv')%>%
  mutate(Meand13C = Meand13C + 1.3,
         Meand15N = Meand15N + 3.3,
         Meand34S = Meand34S + .5)



#Isotope biplots
# C and N
wcn = ggplot(data = sources, aes(Meand13C, Meand15N))+
  geom_point(data = sources, size = 3, pch=c(20))+ 
  geom_errorbar(data = sources, aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0) + 
  geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
  geom_point(data = fish, aes(x = d13C, y = d15N,color = common_name), size=3, pch=c(20))+
  # scale_color_manual(values = cols, drop = F)+
  # scale_x_continuous(limits = c(-27, -6))+
  # scale_y_continuous(limits = c(0,14))+
  theme( legend.title = element_blank(),
         legend.text=element_text(size=12))#,legend.position=c(.85,.15))

#ggsave('flbayCNwet.pdf', units="in", width=10, height=6)

# C and S
wcs = ggplot(data = sources, aes(Meand13C, Meand34S))+
  geom_point(data = fish, aes(x = d13C, y = d34S,color = common_name), size=3, pch=c(20))+
  # scale_color_manual(values = cols, drop = F)+
  geom_point(data = sources, size = 3, pch=c(20))+ 
  geom_errorbar(data = sources, aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
  geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
  # scale_x_continuous(limits = c(-27, -6))+
  # scale_y_continuous(limits = c(-16.5, 24))+
  theme(legend.title = element_blank())#, legend.position=c(.85,.85))
#ggsave('flbayCSwet.pdf', units="in", width=10, height=6)

wcn
wcs


sites = c("SRS3", "SRS4","SRS6", "RB10", "TS3", "TS7", "TS9", "TS10", "TS11")

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
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    geom_point(data = fish, aes(x = d13C, y = d15N,color = common_name), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(0,14))+
    theme( legend.title = element_blank(),
           legend.text=element_text(size=12))#,legend.position=c(.85,.15))
  
  ggsave(paste0(sites[i], 'CN.pdf'), units="in", width=10, height=6)
  
  # C and S
  wcs = ggplot(data = sources, aes(Meand13C, Meand34S))+
    geom_point(data = fish, aes(x = d13C, y = d34S,color = common_name), size=3, pch=c(20))+
    # scale_color_manual(values = cols, drop = F)+
    geom_point(data = sources, size = 3, pch=c(20))+ 
    geom_errorbar(data = sources, aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
    geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
    ylab(expression(paste(delta^{34}, "S (\u2030)")))+
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    theme_classic() + geom_text(data = sources, aes(label = source),hjust=-.1, vjust=-1) +
    # scale_x_continuous(limits = c(-27, -6))+
    # scale_y_continuous(limits = c(-16.5, 24))+
    theme(legend.title = element_blank())#, legend.position=c(.85,.85))
  ggsave(paste0(sites[i], 'CS.pdf'), units="in", width=10, height=6)
  
}

# dry 
not = c('Snook', 'Jack crevale', 'Hard head catfish',
        'Bonefish','Toad fish', 'Blue crab',
        'Gray snapper', 'Spotted seatrout')
fish = read_csv('FLBayMM.csv')%>% 
  filter(Season == 'dry', !(Species %in% not)) %>% 
  mutate(Species = str_replace(Species, "Mojarra", "Silver Jenny mojarra"))

#Isotope biplots
# C and N
dcn = ggplot(data = sources, aes(Meand13C, Meand15N))+
  geom_point(data = sources, size = 3, pch=c(20))+ 
  geom_errorbar(data = sources, aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0) + 
  geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme_classic() + geom_text(data = sources, aes(label = Source),hjust=-.1, vjust=-1) +
  geom_point(data = fish, aes(x = d13C, y = d15N,color = Species), size=3, pch=c(20))+
  scale_color_manual(values = cols, drop = F)+
  scale_x_continuous(limits = c(-27, -6))+
  scale_y_continuous(limits = c(0,14))+
  theme( legend.title = element_blank())#,legend.position=c(.85,.15))

#ggsave('flbayCNdry.pdf', units="in", width=10, height=6)

# C and S
dcs = ggplot(data = sources, aes(Meand13C, Meand34S))+
  geom_point(data = fish, aes(x = d13C, y = d34S,color = Species), size=3, pch=c(20))+
  scale_color_manual(values = cols, drop = F)+
  geom_point(data = sources, size = 3, pch=c(20))+ 
  geom_errorbar(data = sources, aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
  geom_errorbarh(data = sources, aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme_classic() + geom_text(data = sources, aes(label = Source),hjust=-.1, vjust=-1) +
  scale_x_continuous(limits = c(-27, -6))+
  scale_y_continuous(limits = c(-16.5, 24))+
  theme(legend.title = element_blank())#, legend.position=c(.85,.85))

ggarrange(wcn, dcn, wcs, dcs, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          legend = 'bottom', common.legend = T)
ggsave("figs/FLbayBiplotsDO_pom.tiff", units="in", width=10, height=8, dpi=600,compression = 'lzw')
#ggsave('flbayCSdry.pdf', units="in", width=10, height=6)
