# Authors: James Sturges, Ryan Rezek, Ryan James
# Last Updated 19 August 2024
# Seasonal differences in brown/green energy pathways and source specific contributions across multiple FCE food webs
# Data Setup ----

# read in libraries
library(MixSIAR)
library(forcats)
library(car)
library(scales)
library(rgl)
library(flextable)
library(magick)
library(R2jags)
library(rjags)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)

#set print limits high for mixing model outputs
options(max.print = 99999)

# Create mixTable function created by Ryan James -
mixTable = function(file,type,ind = F,nest = F, csv = F){
  require(tidyverse)
  cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '105%', '95%', '910.5%')
  if(csv){
    x = read_csv(file)
  }else{
    x = read_fwf(file, skip = 8)
  }
  names(x) = cn
  x$source = NA
  x$name = NA
  x$code = NA
  
  
  
  if (nest == F){
    for (i in 1:nrow(x)){
      temp = strsplit(x$ID, split = '.', fixed = T)
      x$source[i] = temp[[i]][3]
      x$name[i] = temp[[i]][2]
      
      x$type = type
      x$ymax = x$`105%` + 1.5*(x$`105%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`105%` - x$`25%`)
      
      df = tibble(x$name, x$type, x$source, x$Mean, x$SD, x$`2.5%`, x$`910.5%`,
                  x$`50%`, x$`25%`, x$`105%`, x$ymax, x$ymin)
      colnames(df) = c('name', 'type', 'source', 'mean', 'sd', 'lowend', 'highend',
                       'mid', 'low', 'up', 'ymax', 'ymin')
    }
  }else{
    for (i in 1:nrow(x)){
      temp = strsplit(x$ID, split = '.', fixed = T)
      x$source[i] = temp[[i]][4]
      x$code[i] = temp[[i]][3]
      x$name[i] = temp[[i]][2]
      
      x$type = type
      x$ymax = x$`105%` + 1.5*(x$`105%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`105%` - x$`25%`)
      
      df = tibble(x$name, x$type, x$source, x$code, x$Mean, x$SD, x$`2.5%`, x$`910.5%`,
                  x$`50%`, x$`25%`, x$`105%`, x$ymax, x$ymin)
      colnames(df) = c('name', 'type', 'source', 'code', 'mean', 'sd', 'lowend', 'highend',
                       'mid', 'low', 'up', 'ymax', 'ymin')
    }
  }
  
  for (i in 1:nrow(df)){
    if (df$ymax[i] > df$highend[i]){
      df$ymax[i] = df$highend[i]
    }
    if (df$ymin[i] < df$lowend[i]){
      df$ymin[i] = df$lowend[i]
    }
  }
  df = df %>% drop_na %>%
    filter(name != 'global')
  
  
  if (ind == T){
    if (nest == T){
      df = df %>% select(name, type, code, source, mean) %>%
        pivot_wider(names_from = 'source', values_from = 'mean')
    }else{
      df = df %>% select(name, type, source, mean)%>%
        pivot_wider(names_from = 'source', values_from = 'mean')
    }
  }
  
  return(df)
}

# updated data as of August 2023 includes sulfur reruns and dry 2020 supplemental values
SI <- read.csv("data/FCE_SI_data_13_August_2023.csv",na.strings=c("","NA"))%>%filter(is.na(outlier))

#raw stable isotope values for individual level analysis
SI = SI %>% 
  select(site,hydroseason,year, group, common_name, functional_grp,d13C,d15N,d34S)


#average stable isotope values for species across sites and seasons
SIa<-SI %>% group_by(site, hydroseason, group, common_name, functional_grp) %>% 
  summarise(n=n(),md13C = mean(as.numeric(d13C), na.rm = T),d13Csd=sd(as.numeric(d13C),na.rm = T),
            md15N=mean(as.numeric(d15N),na.rm = T),d15Nsd=sd(as.numeric(d15N),na.rm = T),
            md34S=mean(as.numeric(d34S),na.rm = T),d34Ssd=sd(d34S,na.rm = T)) %>%
  filter(!is.na(md13C),!is.na(md15N),!is.na(md34S))


# averages by species not including seasonality as a factor
SIb<-SI %>% group_by(site,group, common_name, functional_grp) %>% 
  summarise(n=n(),md13C = mean(d13C,na.rm = T),d13Csd=sd(d13C,na.rm = T), md15N=mean(d15N,na.rm = T),d15Nsd=sd(d15N,na.rm = T),md34S=mean(d34S,na.rm = T),d34Ssd=sd(d34S,na.rm = T))  %>%
  filter(!is.na(md13C),!is.na(md15N),!is.na(md34S))



# Biplots ----
# create list of sites as i vector
sites = c("SRS3", "SRS4","SRS6", "RB10", "TS3", "TS7", "TS9", "TS10", "TS11")

# produces two biplots per site (one CN and one CS)
for(i in 1:length(sites)) {
  
  fish = read_csv(paste0('data/Consumers/', sites[i], 'mix.csv'))
  
  sources = read.csv(paste0('data/Sources/sources', sites[i], '.csv')) %>% 
    mutate(Meand13C = Meand13C + 1.95,
           Meand15N = Meand15N + 5.1,
           Meand34S = Meand34S + .75)
  
  
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
  
  ggsave(paste0('figures/nested/biplots/',sites[i], 'CN.pdf'), units="in", width=10, height=6)
  
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
  ggsave(paste0('figures/nested/biplots/',sites[i], 'CS.pdf'), units="in", width=10, height=6)
  
}



# creates the same isotope biplots but colors points by species to look at wet dry trends
for(i in 1:length(sites)) {
  
  fish = read_csv(paste0('data/Consumers/', sites[i], 'mix.csv'))
  
  sources = read.csv(paste0('data/Sources/sources', sites[i], '.csv')) %>% 
    mutate(Meand13C = Meand13C + 1.95,
           Meand15N = Meand15N + 5.1,
           Meand34S = Meand34S + .75)
  
  
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
  
  ggsave(paste0('figures/nested/biplots/',sites[i], 'CN_WD.pdf'), units="in", width=10, height=6)
  
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
  ggsave(paste0('figures/nested/biplots/',sites[i], 'CS_WD.pdf'), units="in", width=10, height=6)
  
}
# RB10 Mixing Models ----

# RB10 with Floc
RB10mix<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")

write.csv(RB10mix,"data/Consumers/nested/RB10mix.csv",row.names = F) 



mix <- load_mix_data(filename="data/Consumers/nested/RB10mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename="data/Sources/sourcesRB10.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_RB10.csv"), mix)

plot_data(filename="figures/nested/isospace/RB10_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "data/Consumers/nested/RB10_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)


#run a test model to make sure it works
# jags.RB10 <- run_model(run="test", mix, source, discr, model_filename,
#                        alpha.prior = 1, resid_err=F, process_err=F)


jags.RB10 <- run_model(run= "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err, process_err)

  output_jags.RB10  <- list(summary_save = T,
                          summary_name = "data/JAGS_Output/RB10/nested/FCERB10_sumstats",
                          sup_post = F,
                          plot_post_save_pdf = T,
                          plot_post_name = "data/JAGS_Output/RB10/nested/FCERB10_plot",
                          sup_pairs = F,
                          plot_pairs_save_pdf = T,
                          plot_pairs_name = "data/JAGS_Output/RB10/nested/FCERB10_pairs",
                          sup_xy = T,
                          plot_xy_save_pdf = T,
                          plot_xy_name = "data/JAGS_Output/RB10/nested/FCERB10_plot",
                          gelman = T,
                          heidel = F,
                          geweke = T,
                          diag_save = T,
                          diag_name = "data/JAGS_Output/RB10/nested/FCERB10_Diagnostic",
                          indiv_effect = F,
                          plot_post_save_png = F,
                          plot_pairs_save_png = F,
                          plot_xy_save_png = F)

saveRDS(jags.RB10, 'data/JAGS_Output/RB10/nested/RB10.RDS')
 
output_JAGS(jags.RB10, mix, source, output_jags.RB10)

mixtable_RB10 = mixTable("data/JAGS_Output/RB10/nested/FCERB10_sumstats.txt",type = "RB10", nest = T)

write.csv(mixtable_RB10, "data/Mix_Quants/nested/MT_RB10.csv", row.names = FALSE)


# combinedRB10 <- combine_sources(jags.RB10, mix, source, alpha.prior=1,
#                                 groups=list(green=c('Epiphytes','Phytoplankton'), brown=c('Mangrove', 'Floc')))
# 
# # get posterior medians for new source groupings
# apply(combinedRB10$post, 2, median)
# 
#  summary_stat(combinedRB10, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
#               filename = "RB10_combined_sumstats" )
#  


# SRS Mixing Models ----
# SRS 3
SRS3mix <- SIa %>% filter(site == 'SRS3', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS3mix, "data/Consumers/nested/SRS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/SRS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS3.csv"), mix)

plot_data(filename = "figures/nested/isospace/SRS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/SRS3_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.SRS3 <- run_model(run = "test", mix, source, discr, model_filename,
#                        alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS3 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err, process_err)

output_jags.SRS3 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS3/nested/FCESRS3_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS3/nested/FCESRS3_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS3/nested/FCESRS3_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "data/JAGS_Output/SRS3/nested/FCESRS3_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS3/nested/FCESRS3_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

saveRDS(jags.SRS3, 'data/JAGS_Output/SRS3/nested/SRS3.RDS')

output_JAGS(jags.SRS3, mix, source, output_jags.SRS3)

mixtable_SRS3 = mixTable("data/JAGS_Output/SRS3/nested/FCESRS3_sumstats.txt", type = "SRS3", nest = TRUE)

write.csv(mixtable_SRS3, "data/Mix_Quants/nested/MT_SRS3.csv", row.names = FALSE)

# combinedSRS3 <- combine_sources(jags.SRS3, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton','Floc'), brown=c('Sawgrass, Periphyton')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS3$post, 2, median)
# # summary_stat(combinedSRS3, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS3_combined_sumstats" )


# SRS4 


SRS4mix <- SIa %>% filter(site == 'SRS4', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS4mix, "data/Consumers/nested/SRS4mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/SRS4mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS4.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS4.csv"), mix)

plot_data(filename = "figures/nested/isospace/SRS4_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/SRS4_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.SRS4 <- run_model(run = "test", mix, source, discr, model_filename,
#                        alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS4 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err, process_err)
output_jags.SRS4 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS4/nested/FCESRS4_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS4/nested/FCESRS4_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS4/nested/FCESRS4_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/SRS4/nested/FCESRS4_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS4/nested/FCESRS4_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

saveRDS(jags.SRS4, 'data/JAGS_Output/SRS4/nested/SRS4.RDS')

output_JAGS(jags.SRS4, mix, source, output_jags.SRS4)

mixtable_SRS4 = mixTable("data/JAGS_Output/SRS4/nested/FCESRS4_sumstats.txt", type = "SRS4", nest = TRUE)

write.csv(mixtable_SRS4, "data/Mix_Quants/nested/MT_SRS4.csv", row.names = FALSE)

# combinedSRS4 <- combine_sources(jags.SRS4, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton',"Epiphytes"), brown=c('Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS4$post, 2, median)
# # summary_stat(combinedSRS4, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS4_combined_sumstats" )

# SRS 6
SRS6mix <- SIa %>% filter(site == 'SRS6', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS6mix, "data/Consumers/nested/SRS6mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/SRS6mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS6.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS6.csv"), mix)

plot_data(filename = "figures/nested/isospace/SRS6_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/SRS6_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.SRS6 <- run_model(run = "test", mix, source, discr, model_filename,
#                        alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS6 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err, process_err)

output_jags.SRS6 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS6/nested/FCESRS6_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS6/nested/FCESRS6_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS6/nested/FCESRS6_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/SRS6/nested/FCESRS6_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS6/nested/FCESRS6_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

saveRDS(jags.SRS6, 'data/JAGS_Output/SRS6/nested/SRS6.RDS')

output_JAGS(jags.SRS6, mix, source, output_jags.SRS6)

mixtable_SRS6 = mixTable("data/JAGS_Output/SRS6/nested/FCESRS6_sumstats.txt", type = "SRS6", nest = TRUE)

write.csv(mixtable_SRS6, "data/Mix_Quants/nested/MT_SRS6.csv", row.names = FALSE)
# combinedSRS6 <- combine_sources(jags.SRS6, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton','Filamentous Green Algae' ), brown=c('Mangrove', 'Red Macroalgae')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS6$post, 2, median)
# # summary_stat(combinedSRS6, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS6_combined_sumstats" )



# TS Mixing Models ----

# TS3

TS3mix <- SIa %>% filter(site == 'TS3', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS3mix, "data/Consumers/nested/TS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/TS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS3.csv"), mix)

plot_data(filename = "figures/nested/isospace/TS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/TS3_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.TS3 <- run_model(run = "test", mix, source, discr, model_filename,
#                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS3 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

output_jags.TS3 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS3/nested/FCETS3_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS3/nested/FCETS3_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS3/nested/FCETS3_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS3/nested/FCETS3_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS3/nested/FCETS3_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

saveRDS(jags.TS3, 'data/JAGS_Output/TS3/nested/TS3.RDS')

output_JAGS(jags.TS3, mix, source, output_jags.TS3)

mixtable_TS3 = mixTable("data/JAGS_Output/TS3/nested/FCETS3_sumstats.txt", type = "TS3", nest = TRUE)

write.csv(mixtable_TS3, "data/Mix_Quants/nested/MT_TS3.csv", row.names = FALSE)
# ##combine posterior ground into brown/green
# combinedTS3 <- combine_sources(jags.TS3, mix, source, alpha.prior=1, 
#                                groups=list(green=c('Periphyton'), brown=c('Dry Sawgrass','Wet Sawgrass', 'Floc' )))
# 
# # get posterior medians for new source groupings
# apply(combinedTS3$post, 2, median)
# summary_stat(combinedTS3, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
#              filename = "TS3_combined_sumstats", na.rm = TRUE)






# TS7 

TS7mix <- SIa %>% filter(site == 'TS7', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS7mix, "data/Consumers/nested/TS7mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/TS7mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS7.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS7.csv"), mix)

plot_data(filename = "figures/nested/isospace/TS7_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS7_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.TS7 <- run_model(run = "test", mix, source, discr, model_filename,
#                       alpha.prior = 1, resid_err, process_err)

jags.TS7 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS7 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS7/nested/FCETS7_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS7/nested/FCETS7_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS7/nested/FCETS7_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS7/nested/FCETS7_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS7/nested/FCETS7_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

saveRDS(jags.TS7, 'data/JAGS_Output/TS7/nested/TS7.RDS')

output_JAGS(jags.TS7, mix, source, output_jags.TS7)

mixtable_TS7 = mixTable("data/JAGS_Output/TS7/nested/FCETS7_sumstats.txt", type = "TS7", nest = TRUE)

write.csv(mixtable_TS7, "data/Mix_Quants/nested/MT_TS7.csv", row.names = FALSE)

# ##combine posterior ground into brown/green
# combinedTS7 <- combine_sources(jags.TS7, mix, source, alpha.prior=1, 
#                                groups=list(green=c('SPOM', 'Epiphytes'), brown=c('Seagrass','Mangrove' )))
# 
# # get posterior medians for new source groupings
# apply(combinedTS7$post, 2, median)
# summary_stat(combinedTS7, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
#              filename = "TS7_combined_sumstats", na.rm = TRUE)




# TS9 

TS9mix <- SIa %>% filter(site == 'TS9', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS9mix, "data/Consumers/nested/TS9mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/TS9mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS9.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS9.csv"), mix)

plot_data(filename = "figures/nested/isospace/TS9_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/TS9_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.TS9 <- run_model(run = "test", mix, source, discr, model_filename,
#                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS9 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

output_jags.TS9 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS9/nested/FCETS9_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS9/nested/FCETS9_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS9/nested/FCETS9_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS9/nested/FCETS9_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS9/nested/FCETS9_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

saveRDS(jags.TS9, 'data/JAGS_Output/TS9/nested/TS9.RDS')

output_JAGS(jags.TS9, mix, source, output_jags.TS9)

mixtable_TS9 = mixTable("data/JAGS_Output/TS9/nested/FCETS9_sumstats.txt", type = "TS9", nest = TRUE)

write.csv(mixtable_TS9, "data/Mix_Quants/nested/MT_TS9.csv", row.names = FALSE)
# ##combine posterior ground into brown/green
# combinedTS9 <- combine_sources(jags.TS9, mix, source, alpha.prior=1, 
#                                groups=list(green=c('Epiphytes', 'SPOM'), brown=c('Seagrass', 'Mangrove' )))
# 
# # get posterior medians for new source groupings
# apply(combinedTS9$post, 2, median)
# # summary_stat(combinedTS9, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
# #              filename = "TS9_combined_sumstats" )
# 


# TS10 

TS10mix <- SIa %>% filter(site == 'TS10', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS10mix, "data/Consumers/nested/TS10mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/TS10mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS10.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS10.csv"), mix)

plot_data(filename = "figures/nested/isospace/TS10_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/TS10_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.TS10 <- run_model(run = "test", mix, source, discr, model_filename,
#                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS10 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

output_jags.TS10 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS10/nested/FCETS10_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS10/nested/FCETS10_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS10/nested/FCETS10_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS10/nested/FCETS10_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS10/nested/FCETS10_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

saveRDS(jags.TS10, 'data/JAGS_Output/TS10/nested/TS10.RDS')

output_JAGS(jags.TS10, mix, source, output_jags.TS10)

mixtable_TS10 = mixTable("data/JAGS_Output/TS10/nested/FCETS10_sumstats.txt", type = "TS10", nest = TRUE)

write.csv(mixtable_TS10, "data/Mix_Quants/nested/MT_TS10.csv", row.names = FALSE)
# ##combine posterior ground into brown/green
# combinedTS10 <- combine_sources(jags.TS10, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Epiphytes', 'SPOM'), brown=c('Seagrass', 'Mangrove' )))
# 
# # get posterior medians for new source groupings
# apply(combinedTS10$post, 2, median)
# # summary_stat(combinedTS10, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
# #              filename = "TS10_combined_sumstats" )


# TS11 Mixing Model 

TS11mix <- SIa %>% filter(site == 'TS11', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS11mix, "data/Consumers/nested/TS11mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/nested/TS11mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors=c('hydroseason','common_name'),
                     fac_random=c(F,T),
                     fac_nested=c(F,T),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS11.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS11.csv"), mix)

plot_data(filename = "figures/nested/isospace/TS11_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/nested/TS11_mix.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# jags.TS11 <- run_model(run = "test", mix, source, discr, model_filename,
#                        alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS11 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err, process_err)

output_jags.TS11 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/TS11/nested/FCETS11_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/TS11/nested/FCETS11_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/TS11/nested/FCETS11_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/TS11/nested/FCETS11_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/TS11/nested/FCETS11_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

saveRDS(jags.TS11, 'data/JAGS_Output/TS11/nested/TS11.RDS')

output_JAGS(jags.TS11, mix, source, output_jags.TS11)

mixtable_TS11 = mixTable("data/JAGS_Output/TS11/nested/FCETS11_sumstats.txt", type = "TS11", nest = TRUE)

write.csv(mixtable_TS11, "data/Mix_Quants/nested/MT_TS11.csv", row.names = FALSE)
# ##combine posterior ground into brown/green
# combinedTS11 <- combine_sources(jags.TS11, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('SPOM','Epiphytes'), brown=c('Seagrass', 'Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedTS11$post, 2, median)
# summary_stat(combinedTS11, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
#              filename = "TS11_combined_sumstats" )


# Extracting Combined Posterior Distribution ----
post = function(jags,mix_file,source_file,site){
  sims_matrix =  jags$BUGSoutput$sims.matrix %>%
    as_tibble() %>%
    select(contains("p.fac1")) %>%
    mutate(site = site)
  
  
  mix <- load_mix_data(filename = mix_file,
                       iso_names=c("d13C","d34S"),
                       factors=c('hydroseason','common_name'),
                       fac_random=c(F,T),
                       fac_nested=c(F,T),
                       cont_effects=NULL)
  
  mix_levels = mix$FAC[[1]]$labels
  
  source <- load_source_data(filename=source_file,
                             source_factors=NULL,
                             conc_dep=T,
                             data_type="means",
                             mix)
  
  source_vect = source$source_names
  
  new_names <- names(sims_matrix) %>%
    str_replace_all("p.fac1\\[(\\d+),(\\d+)\\]",
                    \(x){
                      match <- str_match(x, "p.fac1\\[(\\d+),(\\d+)\\]")
                      vec1_value <- mix_levels[as.numeric(match[2])]
                      vec2_value <- source_vect[as.numeric(match[3])]
                      paste(vec1_value,vec2_value,sep = '_')})
  
  names(sims_matrix) <- new_names
  df = sims_matrix %>%
    pivot_longer(-site,names_pattern = "(.*)_(.*)", names_to = c("Season", "Source"), values_to = "Source Contribution")
}

site_names <- c("RB10", "SRS4", "SRS6", "SRS3", "TS3", "TS7", "TS9", "TS10", "TS11")

# Define a list of corresponding JAGS output files
jags_files <- list(
  RB10 = readRDS('data/JAGS_Output/RB10/nested/RB10.RDS'),
  SRS4 = readRDS('data/JAGS_Output/SRS4/nested/SRS4.RDS'),
  SRS6 = readRDS('data/JAGS_Output/SRS6/nested/SRS6.RDS'),
  SRS3 = readRDS('data/JAGS_Output/SRS3/nested/SRS3.RDS'),
  TS3 = readRDS('data/JAGS_Output/TS3/nested/TS3.RDS'),
  TS7 = readRDS('data/JAGS_Output/TS7/nested/TS7.RDS'),
  TS9 = readRDS('data/JAGS_Output/TS9/nested/TS9.RDS'),
  TS10 = readRDS('data/JAGS_Output/TS10/nested/TS10.RDS'),
  TS11 = readRDS('data/JAGS_Output/TS11/nested/TS11.RDS')
)

# Define the paths to the mix and source files
mix_files <- list(
  RB10 = "data/Consumers/nested/RB10mix.csv",
  SRS4 = "data/Consumers/nested/SRS4mix.csv",
  SRS6 = "data/Consumers/nested/SRS6mix.csv",
  SRS3 = "data/Consumers/nested/SRS3mix.csv",
  TS3 = "data/Consumers/nested/TS3mix.csv",
  TS7 = "data/Consumers/nested/TS7mix.csv",
  TS9 = "data/Consumers/nested/TS9mix.csv",
  TS10 = "data/Consumers/nested/TS10mix.csv",
  TS11 = "data/Consumers/nested/TS11mix.csv"
)

source_files <- list(
  RB10 = "data/Sources/sourcesRB10.csv",
  SRS4 = "data/Sources/sourcesSRS4.csv",
  SRS6 = "data/Sources/sourcesSRS6.csv",
  SRS3 = "data/Sources/sourcesSRS3.csv",
  TS3 = "data/Sources/sourcesTS3.csv",
  TS7 = "data/Sources/sourcesTS7.csv",
  TS9 = "data/Sources/sourcesTS9.csv",
  TS10 = "data/Sources/sourcesTS10.csv",
  TS11 = "data/Sources/sourcesTS11.csv"
)

# Initialize a list to store results
posterior_tibbles <- list()

# Loop through each site and apply the post function
for (site in site_names) {
  jags <- jags_files[[site]]
  mix_file <- mix_files[[site]]
  source_file <- source_files[[site]]
  
  posterior_tibbles[[site]] <- post(jags, mix_file = mix_file, source_file = source_file, site = site)
}

# Combine all the results into one tibble if needed
combined_posterior_tibble <- bind_rows(posterior_tibbles, .id = "site")

# Posterior Channel Boxplot ----
combined_posterior_tibble = combined_posterior_tibble %>%
  mutate(path = case_when(
    Source %in% c("Epiphytes", "Phytoplankton", "Filamentous Green Algae", "Periphyton", "Epiphytic microalgae", 'SPOM') ~ "green",
    Source %in% c("Mangrove", "Sawgrass", "Red Macroalgae", "Seagrass", 'Floc') ~ "brown",
    TRUE ~ NA_character_
  ))
combined_posterior_tibble = combined_posterior_tibble %>%
  mutate(transect = case_when(
    site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
    site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"))

combined_posterior_tibble <- combined_posterior_tibble %>%
  mutate(site = ifelse(site == "RB10", "Upper River",
                       ifelse(site == "SRS3", "SRS Marsh",
                              ifelse(site == "SRS4", "Mid River",
                                     ifelse(site == "SRS6", "Lower River",
                                            ifelse(site == "TS3", "TS Marsh",
                                                   ifelse(site == "TS7", "Mangrove Ecotone",
                                                          ifelse(site == "TS9", "Inner Bay",
                                                                 ifelse(site == "TS10", "Mid Bay",
                                                                        ifelse(site == "TS11", "Outer Bay", site)
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
                       )
  ),
  site = factor(site, levels = c( "SRS Marsh","Upper River","Mid River","Lower River", "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")))

# Create a simulated iteration identifier based on the repeating pattern
combined_posterior_tibble <- combined_posterior_tibble %>%
  group_by(site, Season, Source) %>%
  mutate(iteration = rep(1:3000, each = 1, length.out = n())) %>%
  ungroup()


# Sum contributions by source for each iteration
summed_contributions <- combined_posterior_tibble %>%
  group_by(site, Season, path, transect, iteration) %>%
  summarize(total_contribution = sum(`Source Contribution`), .groups = 'drop')


green_contributions <- summed_contributions %>%
  filter(path == "green")%>%
  group_by(site, Season) %>%
  mutate(color_fill = mean(total_contribution)) 


green_contributions <- green_contributions %>%
  mutate(site = factor(site, levels = c( "Lower River","Mid River", "Upper River","SRS Marsh","Outer Bay", "Mid Bay","Inner Bay", "Mangrove Ecotone", "TS Marsh")))


# Create the boxplot
mixoutput_bxplt_gb_combined <- ggplot(green_contributions, aes(x = site, y = total_contribution, fill = color_fill)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 19, width = 0.5) +
  theme_bw() +
  scale_fill_gradient2(
    low = "#663300",
    high = "#92D050",
    mid = 'white',
    midpoint = 0.5,
    limits = c(0, 1),
    na.value = "grey50"
  ) +
  facet_grid(transect ~ Season, scales = "free_y") +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20, colour = "black"), 
    axis.text.x = element_text(size = 18, colour = "black"), 
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    y = "Green Pathway Source Contribution",
    x = NULL,
    fill = "Mean Green Pathway\nSource Contribution\n "
  ) +
  coord_flip()

mixoutput_bxplt_gb_combined

ggsave("figures/nested/CNS_cplot.png", width = 12, height = 6, dpi = 600)



# Posterior Boxplot No Outliers

mixoutput_bxplt_gb_combined_out <- ggplot(green_contributions, aes(x = site, y = total_contribution, fill = color_fill)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  # Remove outliers by setting outlier.shape to NA
  theme_bw() +
  scale_fill_gradient2(
    low = "#663300",
    high = "#92D050",
    mid = 'white',
    midpoint = 0.5,
    limits = c(0, 1),
    na.value = "grey50"
  ) +
  facet_grid(transect ~ Season, scales = "free_y") +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20, colour = "black"), 
    axis.text.x = element_text(size = 18, colour = "black"), 
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    y = "Green Pathway Source Contribution",
    x = NULL,
    fill = "Mean Green Pathway\nSource Contribution\n "
  ) +
  coord_flip()

mixoutput_bxplt_gb_combined_out
ggsave("figures/nested/CNS_cplot_out.png", width = 12, height = 6, dpi = 600)

# Calculate summary statistics for each site, season, and path
summary_stats <- summed_contributions %>%
  group_by(site, Season, path, transect) %>%
  summarize(
    mean_contribution = mean(total_contribution),
    median_contribution = median(total_contribution),
    sd_contribution = sd(total_contribution),
    min_contribution = min(total_contribution),
    max_contribution = max(total_contribution),
    q25_contribution = quantile(total_contribution, 0.25),
    q75_contribution = quantile(total_contribution, 0.75),
    .groups = 'drop'
  )


mean_source = combined_posterior_tibble %>%
  group_by(site, Season, Source, path, transect) %>%
  summarize(mean_contribution = mean(`Source Contribution`)) %>%
  ungroup()

# Calculate mean across all iterations for each site, season, and path
mean_channel <- summed_contributions %>%
  group_by(site, Season, path, transect) %>%
  summarize(mean_contribution = mean(total_contribution), .groups = 'drop')

mean_channel_combined <- mean_channel %>%
  left_join(summary_stats, by = c("site", "Season", "path", "transect"))


y_label_formatter <- function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}


mean_channel_wide <- mean_channel %>%
  pivot_wider(
    names_from = path,
    values_from = mean_contribution,
    names_prefix = "", # or use "contribution_" if you prefer a prefix
    values_fill = list(mean_contribution = NA) # Fill missing values with NA
  ) %>%
  rename(
    `brown contribution` = brown,
    `green contribution` = green
  )


mean_source_wide <- mean_source %>%
  pivot_wider(
    names_from = Source,
    values_from = mean_contribution,
    names_prefix = "",
    values_fill = list(mean_contribution = NA) # Fill missing values with NA
  )

channel_df = summary_stats %>% 
  filter(path == as.character("green"))

# Posterior Source Boxplot ----
cont_df <- combined_posterior_tibble
#colors and labels for background source color fill in panels
s_df= tibble(Source = c("Sawgrass", "Mangrove",
                        "Floc","Red Macroalgae","Seagrass",
                        "Epiphytes","Periphyton","Filamentous Green Algae",
                        'SPOM',"Phytoplankton"),
             lab = c('Sawgrass', 'Mang',
                     "Floc","RMA","Seagrass",
                     'EMA','Peri','FGA',
                     'POM', 'Phyto'),
             s = c('b1', 'b1',
                   'b2','b2','b2',
                   'g1','g1','g1',
                   'g2','g2'),
             x = c(1,1,
                   2,2,2,
                   3,3,3,
                   4,4),
             g = c('brown','brown',
                   'brown','brown','brown',
                   'green','green','green',
                   'green','green')) |> 
  mutate(lab = factor(lab, levels = c('Sawgrass', 'Mang',
                                      "Floc","RMA","Seagrass",
                                      'EMA','Peri','FGA',
                                      'POM', 'Phyto')),
         lb = factor(lab, levels = c( "Sawgrass","Mang","Floc","RMA",
                                      "Seagrass", "EMA","Peri" ,
                                      "Green Algae","POM", "Phyto")))

cont_df = cont_df |> 
  as_tibble() |> 
  left_join(s_df, by = 'Source') 


y_label_formatter <- function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

cont_df <- cont_df %>%
  mutate(lab = case_when(
    lab == "Phyto" ~ "PMA",
    TRUE ~ lab
  ))

cont_df <- cont_df %>%
  mutate(lab = case_when(
    lab == "POM" ~ "PMA",
    TRUE ~ lab
  ))

cont_df <- cont_df %>%
  mutate(lab = factor(lab, levels = c('Sawgrass', 'Mang',
                                      "Floc","RMA","Seagrass",
                                      'EMA','Peri','FGA',
                                      'PMA')))


cont_df = cont_df %>% 
  mutate(s_cont = cont_df$`Source Contribution`)


source_cont_plot = ggplot(data = cont_df, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_df |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_df |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site, scales = "free_x") +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Source contribution",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


source_cont_plot
ggsave("figures/nested/CNS_source_cont_plot_GBcol.png", width = 9.5, height = 7.5, dpi = 600)

source_cont_plot_no_outliers = ggplot(data = cont_df, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot(outlier.shape = NA) + # Set outlier.shape to NA to remove outliers
  geom_rect(data = cont_df |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_df |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot(outlier.shape = NA) + # Also remove outliers here
  theme_bw() +
  facet_wrap(~site, scales = "free_x") +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Source contribution",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# Print the updated plot without outliers
source_cont_plot_no_outliers
ggsave("figures/nested/CNS_source_cont_plot_no_outliers.png", width = 9.5, height = 7.5, dpi = 600)



# SRS only tall source contribution plot ----
cont_SRS = cont_df %>% 
  filter(transect == "Shark River Slough")

SRS_sources = ggplot(data = cont_SRS, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_SRS |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_SRS |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site, scales = "free_x", drop = T, ncol = 1) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

SRS_sources
ggsave("figures/nested/CNS_Source_SRS_tall.png", width = 4, height = 8, dpi = 600)


SRS_sources_no_outlier = ggplot(data = cont_SRS, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot(outlier.shape = NA) + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_SRS |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_SRS |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  facet_wrap(~site, scales = "free_x", drop = T, ncol = 1) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

SRS_sources_no_outlier
ggsave("figures/nested/CNS_Source_SRS_tall_no_outliers.png", width = 4, height = 8, dpi = 600)

# TS only tall source contribution plot ----

cont_TS = cont_df %>% 
  filter(transect == "Taylor Slough")

TS_sources = ggplot(data = cont_TS, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_TS |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_TS |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site, scales = "free_x", drop = T, ncol = 1) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

TS_sources
ggsave("figures/nested/CNS_Source_TS_tall.png", width = 4, height = 8, dpi = 600)


TS_sources_no_outlier = ggplot(data = cont_TS, aes(x = lab, y = s_cont, fill = Season)) +
  geom_boxplot(outlier.shape = NA) + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_TS |> group_by(site) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_TS |> group_by(site) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  facet_wrap(~site, scales = "free_x", drop = T, ncol = 1) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

TS_sources_no_outlier
ggsave("figures/nested/CNS_Source_TS_tall_no_outliers.png", width = 4, height = 8, dpi = 600)
# supplemental mixing space biplots-----
# consumer <- read_csv('data/Consumers/combined.csv') %>%
#   select(species) %>%
#   distinct(species, .keep_all = TRUE) %>% 
#   arrange(species)


cols <- c(
  "Bluegill"                = 'lightblue',
  "Coastal Shiner"          = 'lightgreen',
  "Dollar Sunfish"          = 'goldenrod',
  "Eastern Grass Shrimp"    = 'lightpink',
  "Flatback Mud Crab"       = 'salmon',
  "Golden Topminnow"        = 'darkseagreen4',
  "Hogchoker"               = 'burlywood4',
  "Peacock Eel"             = 'darkorange',
  "Redear Sunfish"          = 'orangered',
  "Tidewater Silverside"    = 'lightcyan',
  "Bluefin Killifish"       = 'dodgerblue',
  "Eastern Mosquitofish"    = 'khaki',
  "Flagfish"                = 'forestgreen',
  "Gar"                     = 'darkslategray',
  "Slough Crayfish"         = 'brown',
  "Warmouth"                = 'chocolate',
  "Dark False Mussel"       = 'antiquewhite4',
  "Harris Mud Crab"         = 'firebrick',
  "Pink Shrimp"             = 'pink',
  "Redfish"                 = 'red',
  "Snook"                   = 'darkgoldenrod',
  "Water Strider"           = 'slateblue',
  "Bay Anchovy"             = 'deepskyblue',
  "Crested Goby"            = 'lightcoral',
  "Mangrove Tree Crab"      = 'olivedrab4',
  "Porceline Crab"          = 'lavender',
  "Salt Marsh Mud Crab"     = 'lightgray',
  "Tidewater Mojarra"       = 'plum',
  "Bowfin"                  = 'indianred',
  "Clown Goby"              = 'darkseagreen',
  "Largemouth Bass"         = 'mediumseagreen',
  "Spotted Sunfish"         = 'chartreuse3',
  "Spotted Tilapia"         = 'powderblue',
  "Striped Mullet"          = 'deepskyblue4',
  "Giant Water Bug"         = 'mediumvioletred',
  "Mayan Cichlid"           = 'mediumorchid',
  "Sailfin Molly"           = 'darkturquoise',
  "Walking Catfish"         = 'yellowgreen',
  "Western Mosquitofish"    = 'lightsalmon',
  "Blue Crab"               = 'royalblue',
  "Striped Mojarra"         = 'lightgoldenrod',
  "Brown Grass Shrimp"      = 'darkkhaki',
  "Silver Jenny"            = 'azure4',
  "Dwarf Seahorse"          = 'mediumslateblue',
  "Hermit Crab"             = 'burlywood',
  "Snapping Shrimp"         = 'violetred4',
  "Spider Crab"             = 'tomato',
  "Spiny Lobster"           = 'sandybrown',
  "Stone Crab"              = 'sienna',
  "Fringed Filefish"        = 'cadetblue',
  "Eastern Oyster"          = 'lightsteelblue',
  "Gulf Pipefish"           = 'skyblue',
  "Mussel"                  = 'mistyrose',
  "Rainwater Killifish"     = 'darkblue',
  "Redfin Needlefish"       = 'darkviolet',
  "Toadfish"                = 'peru',
  "Atlantic Modulus"        = 'thistle3',
  "Barracuda"               = 'slategray',
  "Cross-Barred Venus"      = 'gold',
  "Goldspotted Killifish"   = 'lightyellow',
  "Scaly Pearl Oyster"      = 'orchid',
  "West Indian False Cerith" = 'hotpink4',
  "Arrow Shrimp"            = 'mediumturquoise',
  "Blue Striped Grunt"      = 'grey69',
  "Bubble Snail"            = 'lightseagreen',
  "Chestnut Turban Snail"   = 'cyan',
  "Giant Purple Barnacle"   = 'darkorchid',
  "Ivory Cerith"            = 'ivory',
  "Penaeid Shrimp"          = 'peachpuff',
  "Schoolmaster Snapper"    = 'darkred',
  "Asian Swamp Eel"         = 'darkgreen',
  "Code Goby"               = 'plum4',
  "Pilchard"                = 'antiquewhite',
  "Atlantic Marginella"     = 'lightgoldenrod4'
)
# # Initialize empty lists to store plot objects
# cn_plots <- list()
# cs_plots <- list()
# # Create a dummy data frame with all species for the legend
# all_species <- unique(unlist(sapply(cn_plots, function(p) p$mapping$color)))
# 
# # Define the sites
# sites <- c("SRS3", "SRS4", "SRS6", "RB10", "TS3", "TS7", "TS9", "TS10", "TS11")
# 
# site_names <- c(
#   "SRS3" = "SRS Marsh",
#   "RB10" = "Upper River",
#   "SRS4" = "Mid River",
#   "SRS6" = "Lower River",
#   "TS3" = "TS Marsh",
#   "TS7" = "Mangrove Ecotone",
#   "TS9" = "Inner Bay",
#   "TS10" = "Mid Bay",
#   "TS11" = "Outer Bay"
# )
# 
# # Loop through each site and create plots
# for(i in 1:length(sites)) {
#   
#   fish <- read_csv(paste0('data/Consumers/', sites[i], 'mix.csv'))
#   
#   sources <- read_csv(paste0('data/Sources/sources', sites[i], '.csv')) %>% 
#     mutate(Meand13C = Meand13C + 1.95,
#            Meand15N = Meand15N + 5.1,
#            Meand34S = Meand34S + .75)
#   
#   # Create CN biplot
#   cn_plot <- ggplot(data = sources, aes(Meand13C, Meand15N)) +
#     geom_point(size = 3, pch = 20) + 
#     geom_errorbar(aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0) + 
#     geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
#     geom_point(data = fish, aes(x = d13C, y = d15N, color = common_name), size = 3, pch = 20) +
#     scale_color_manual(values = cols, drop = FALSE) +
#     ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
#     xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
#     labs(title = paste("CN ", site_names[sites[i]])) +
#     theme_classic() + 
#     theme(legend.position = "none") +
#     geom_text(data = sources, aes(label = source), hjust = -0.1, vjust = -1)
#   
#   cn_plots[[sites[i]]] <- cn_plot
#   
#   # Create CS biplot
#   cs_plot <- ggplot(data = sources, aes(Meand13C, Meand34S)) +
#     geom_point(size = 3, pch = 20) + 
#     geom_errorbar(aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0) + 
#     geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0) +
#     geom_point(data = fish, aes(x = d13C, y = d34S, color = common_name), size = 3, pch = 20) +
#     scale_color_manual(values = cols, drop = FALSE) +
#     ylab(expression(paste(delta^{34}, "S (\u2030)"))) +
#     xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
#     labs(title = paste("CS ", site_names[sites[i]])) +
#     theme_classic() + 
#     theme(legend.position = "none") +
#     geom_text(data = sources, aes(label = source), hjust = -0.1, vjust = -1)
#   
#   cs_plots[[sites[i]]] <- cs_plot
# }
# # Create a dummy plot to extract legend
# dummy_plot <- ggplot(data = fish, aes(x = d13C, y = d15N, color = common_name)) +
#   geom_point(size = 3, pch = 20) +
#   scale_color_manual(values = cols, drop = FALSE) +
#   theme_void() +
#   theme(legend.position = "bottom")
# 
# # Extract the legend from the dummy plot
# legend <- get_legend(dummy_plot)
# 
# # Combine all CN plots into one facet grid plot
# combined_cn <- ggarrange(plotlist = cn_plots, ncol = 3, nrow = 3)
# 
# # Combine all CS plots into one facet grid plot
# combined_cs <- ggarrange(plotlist = cs_plots, ncol = 3, nrow = 3)
# 
# # Add the legend to the combined plots
# final_cn <- plot_grid(combined_cn, legend, ncol = 1, rel_heights = c(10, 1))
# final_cs <- plot_grid(combined_cs, legend, ncol = 1, rel_heights = c(10, 1))
# # Save or display the final plots
# ggsave("figures/CS/biplots/combined_CN_with_legend.pdf", final_cn, units = "in", width = 15, height = 15)
# ggsave("figures/CS/biplots/combined_CS_with_legend.pdf", final_cs, units = "in", width = 15, height = 15)




sites <- c("SRS3","RB10", "SRS4", "SRS6",  "TS3", "TS7", "TS9", "TS10", "TS11")

site_order <- c("SRS Marsh", "Upper River", "Mid River", "Lower River", 
                "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")
# Placeholder for the path where your files are stored
path <- "data/Sources/"

# Initialize an empty dataframe
sources <- data.frame()

# Loop through each site to read the corresponding CSV file and bind them together
for (site in sites) {
  # Construct the file name
  file_name <- paste0(path, "sources", site, ".csv")
  
  # Read the file
  temp_df <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Optionally add a column to indicate the site (if not already included)
  temp_df$Site <- site
  
  # Bind the data
  sources <- rbind(sources, temp_df)
}

consumers = SIa %>% 
  filter(group %in% "Consumer")

consumers = consumers %>% 
  mutate(Site = site,
         reps = n) %>% 
  select(!site)

consumers = consumers %>% 
  mutate(zone = ifelse(Site == "RB10", "Upper River",
                       ifelse(Site == "SRS3", "SRS Marsh",
                              ifelse(Site == "SRS4", "Mid River",
                                     ifelse(Site == "SRS6", "Lower River",
                                            ifelse(Site == "TS3", "TS Marsh",
                                                   ifelse(Site == "TS7", "Mangrove Ecotone",
                                                          ifelse(Site == "TS9", "Inner Bay",
                                                                 ifelse(Site == "TS10", "Mid Bay",
                                                                        ifelse(Site == "TS11", "Outer Bay", site)
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
                       )
  ),
  zone = factor(zone, levels = c( "SRS Marsh","Upper River","Mid River","Lower River", "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")))

sources = sources %>% 
  mutate(zone = ifelse(Site == "RB10", "Upper River",
                       ifelse(Site == "SRS3", "SRS Marsh",
                              ifelse(Site == "SRS4", "Mid River",
                                     ifelse(Site == "SRS6", "Lower River",
                                            ifelse(Site == "TS3", "TS Marsh",
                                                   ifelse(Site == "TS7", "Mangrove Ecotone",
                                                          ifelse(Site == "TS9", "Inner Bay",
                                                                 ifelse(Site == "TS10", "Mid Bay",
                                                                        ifelse(Site == "TS11", "Outer Bay", site)
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
                       )
  ),
  zone = factor(zone, levels = c( "SRS Marsh","Upper River","Mid River","Lower River", "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")))

sources = sources %>% 
  mutate(plot_source = ifelse(source == "Sawgrass", "Sawgrass",
                              ifelse(source == "Mangrove", "Mang",
                                     ifelse(source == "Floc", "Floc",
                                            ifelse(source == "Red Macroalgae", "RMA",
                                                   ifelse(source == "Seagrass", "Seagrass",
                                                          ifelse(source == "Epiphytes", "EMA",
                                                                 ifelse(source == "Periphyton", "Peri",
                                                                        ifelse(source == "Filamentous Green Algae", "FGA",
                                                                               ifelse(source == "SPOM", "PMA",
                                                                                      ifelse(source == "Phytoplankton", "PMA", source)
                                                                               )
                                                                        )
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
  ),
  plot_source = factor(plot_source, levels = c( "Sawgrass","Mang","Floc","RMA", "Seagrass", "EMA", "Peri", "FGA", "PMA")))



sources = sources %>% 
  mutate(transect = case_when(
    Site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
    Site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"))

sources <- sources %>%
  mutate(
    Meand13C = Meand13C + 1.95,
    Meand15N = Meand15N + 5.1,
    Meand34S = Meand34S + 0.75
  )

# CN Biplots free scales
cn_plot <- ggplot(data = sources, aes(Meand13C, Meand15N)) +
  geom_point(size = 3, pch = 20) + 
  geom_errorbar(aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0,linetype = "dashed") + 
  geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0,linetype = "dashed") +
  geom_point(data = consumers, aes(x = md13C, y = md15N, color = common_name), size = 3, pch = 20) +
  scale_color_manual(values = cols, drop = FALSE) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(color = "Consumer Species") +
  facet_wrap(~zone, labeller = labeller(zone = site_names), scales = "free") +
  theme_classic() + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    strip.text = element_text(size = 28),
    strip.background = element_rect(fill = "gray"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  ) +
  geom_text(data = sources, aes(label = plot_source),fontface = "bold", size = 4, hjust = -0.04, vjust = -.5) + 
  guides(color = guide_legend(ncol = 1))

# CS Biplots free scales
cs_plot <- ggplot(data = sources, aes(Meand13C, Meand34S)) +
  geom_point(size = 3, pch = 20) + 
  geom_errorbar(aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0,linetype = "dashed") + 
  geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0, linetype = "dashed") +
  geom_point(data = consumers, aes(x = md13C, y = md34S, color = common_name), size = 4, pch = 20) +
  scale_color_manual(values = cols, drop = FALSE) +
  ylab(expression(paste(delta^{34}, "S (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(color = "Consumer Species") +
  facet_wrap(~zone, labeller = labeller(zone = site_names), scales = "free") +
  theme_classic() + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    strip.text = element_text(size = 28),
    strip.background = element_rect(fill = "gray"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  ) +
  geom_text(data = sources, aes(label = plot_source),fontface = "bold", size = 4, hjust = -0.04, vjust = -.5) + 
  guides(color = guide_legend(ncol = 1))
# Display and save plots
cn_plot
ggsave("figures/nested/CN_biplots.png", cn_plot, units = "in", width = 18, height = 15)

cs_plot
ggsave("figures/nested/CS_biplots.png", cs_plot, units = "in", width = 18, height = 15)

# CN Biplots fixed scales
cn_plot_standard_axes <- ggplot(data = sources, aes(Meand13C, Meand15N)) +
  geom_point(size = 3, pch = 20) + 
  geom_errorbar(aes(ymin = Meand15N - SDd15N, ymax = Meand15N + SDd15N), width = 0,linetype = "dashed") + 
  geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0,linetype = "dashed") +
  geom_point(data = consumers, aes(x = md13C, y = md15N, color = common_name), size = 3, pch = 20) +
  scale_color_manual(values = cols, drop = FALSE) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(color = "Consumer Species") +
  facet_wrap(~zone, labeller = labeller(zone = site_names)) +
  theme_classic() + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    strip.text = element_text(size = 28),
    strip.background = element_rect(fill = "gray"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  ) +
  geom_text(data = sources, aes(label = plot_source),fontface = "bold", size = 4, hjust = -0.04, vjust = -.5) + 
  guides(color = guide_legend(ncol = 1))



#CS biplots fixed scales
cs_plot_standard_axes <- ggplot(data = sources, aes(Meand13C, Meand34S)) +
  geom_point(size = 3, pch = 20) + 
  geom_errorbar(aes(ymin = Meand34S - SDd34S, ymax = Meand34S + SDd34S), width = 0,linetype = "dashed") + 
  geom_errorbarh(aes(xmin = Meand13C - SDd13C, xmax =  Meand13C + SDd13C), height = 0, linetype = "dashed") +
  geom_point(data = consumers, aes(x = md13C, y = md34S, color = common_name), size = 4, pch = 20) +
  scale_color_manual(values = cols, drop = FALSE) +
  ylab(expression(paste(delta^{34}, "S (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(color = "Consumer Species") +
  facet_wrap(~zone, labeller = labeller(zone = site_names)) +
  theme_classic() + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    strip.text = element_text(size = 28),
    strip.background = element_rect(fill = "gray"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  ) +
  geom_text(data = sources, aes(label = plot_source),fontface = "bold", size = 4, hjust = -0.04, vjust = -.5) + 
  guides(color = guide_legend(ncol = 1))



cn_plot_standard_axes
ggsave("figures/nested/CN_biplots_fixed.png", cn_plot_standard_axes, units = "in", width = 18, height = 15)

cs_plot_standard_axes
ggsave("figures/nested/CS_biplots_fixed.png", cs_plot_standard_axes, units = "in", width = 18, height = 15)
