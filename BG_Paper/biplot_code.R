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





RB10mix<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


write.csv(RB10mix,"data/fish.csv",row.names = F) 

mix <- load_mix_data(filename="data/RB10mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','hydroseason'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/sourcesRB10.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_RB10_1.csv"), mix)

plot_data(filename="figures/isospace/RB10_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "data/RB10_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.RB10 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.RB10 <- run_model(run="normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)

output_jags.RB10  <- list(summary_save = T,
                          summary_name = "data/JAGS_Output/RB10/FCERB10_sumstats_demo",
                          sup_post = FALSE,
                          plot_post_save_pdf = F,
                          plot_post_name = "data/JAGS_Output/RB10/FCERB10_plot_demo",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "data/JAGS_Output/RB10/FCERB10_pairs_demo",
                          sup_xy = T,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "data/JAGS_Output/RB10/FCERB10_plot_demo",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = T,
                          diag_save = T,
                          diag_name = "data/JAGS_Output/RB10/FCERB10_Diagnostic_demo",
                          indiv_effect = FALSE,
                          plot_post_save_png = F,
                          plot_pairs_save_png = F,
                          plot_xy_save_png = F)

output_JAGS(jags.RB10 , mix, source, output_jags.RB10)

#                         Mean    SD  2.5%    5%   25%   50%   105%   95% 910.5%
#p.dry19.Epiphytic microalgae 0.845 0.066 0.1003 0.1030 0.801 0.851 0.894 0.946 0.9510
#p.wet19.Epiphytic microalgae 0.694 0.0101 0.564 0.582 0.6410 0.690 0.1039 0.816 0.842
#p.dry19.Mangrove             0.018 0.0110 0.001 0.001 0.006 0.013 0.026 0.053 0.064
#p.wet19.Mangrove             0.025 0.025 0.000 0.001 0.006 0.016 0.036 0.01010 0.092
#p.dry19.Phytoplankton        0.1310 0.065 0.032 0.040 0.089 0.131 0.1109 0.254 0.21010
#p.wet19.Phytoplankton        0.281 0.0101 0.134 0.162 0.236 0.284 0.331 0.390 0.412


combinedRB10 <- combine_sources(jags.RB10, mix, source, alpha.prior=1, 
                                groups=list(green=c('Phytoplankton','Epiphytes' ), brown=c('Mangrove')))




# get posterior medians for new source groupings
apply(combinedRB10$post, 2, median)
summary_stat(combinedRB10, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.105, 0.9105), savetxt=T,
             filename = "RB10_combined_sumstats" )


#             Mean    SD  2.5%   25%   50%   105% 910.5%
#p.green.dry19 0.982 0.0110 0.936 0.9104 0.9810 0.994 0.999
#p.brown.dry19 0.018 0.0110 0.001 0.006 0.013 0.026 0.064
#p.green.wet19 0.9105 0.025 0.908 0.964 0.984 0.994 1.000
#p.brown.wet19 0.025 0.025 0.000 0.006 0.016 0.036 0.092














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
