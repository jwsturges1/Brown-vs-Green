# Authors: James Sturges, Ryan Rezek, Ryan James
# Last Updated 20 August 2023
# Seasonal differences in brown/green energy pathways and source specific contributions across multiple FCE food webs
# DATASETUP ----

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


# not including seasonality as a factor
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
# RB10 Mixing Models ----

# RB10 with Floc
RB10mix<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")

write.csv(RB10mix,"data/Consumers/RB10mix.csv",row.names = F) 



mix <- load_mix_data(filename="data/Consumers/RB10mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','hydroseason'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/Sources/sourcesRB10.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_RB10.csv"), mix)

plot_data(filename="figures/isospace/RB10_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "data/Consumers/RB10_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.RB10 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=F)


jags.RB10 <- run_model(run="normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=F)

  output_jags.RB10  <- list(summary_save = T,
                          summary_name = "data/JAGS_Output/RB10/FCERB10_sumstats",
                          sup_post = F,
                          plot_post_save_pdf = T,
                          plot_post_name = "data/JAGS_Output/RB10/FCERB10_plot",
                          sup_pairs = F,
                          plot_pairs_save_pdf = T,
                          plot_pairs_name = "data/JAGS_Output/RB10/FCERB10_pairs",
                          sup_xy = T,
                          plot_xy_save_pdf = T,
                          plot_xy_name = "data/JAGS_Output/RB10/FCERB10_plot",
                          gelman = T,
                          heidel = F,
                          geweke = T,
                          diag_save = T,
                          diag_name = "data/JAGS_Output/RB10/FCERB10_Diagnostic",
                          indiv_effect = F,
                          plot_post_save_png = F,
                          plot_pairs_save_png = F,
                          plot_xy_save_png = F)

 
output_JAGS(jags.RB10, mix, source, output_jags.RB10)

mixtable_RB10 = mixTable("data/JAGS_Output/RB10/FCERB10_sumstats.txt",type = "RB10", nest = T)

write.csv(mixtable_RB10, "data/Mix_Quants/MT_RB10.csv", row.names = FALSE)


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

write.csv(SRS3mix, "data/Consumers/SRS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/SRS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS3.csv"), mix)

plot_data(filename = "figures/isospace/SRS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/SRS3_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS3 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS3 <- run_model(run = "normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.SRS3 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS3/FCESRS3_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS3/FCESRS3_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS3/FCESRS3_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "data/JAGS_Output/SRS3/FCESRS3_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS3/FCESRS3_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.SRS3, mix, source, output_jags.SRS3)

mixtable_SRS3 = mixTable("data/JAGS_Output/SRS3/FCESRS3_sumstats.txt", type = "SRS3", nest = TRUE)

write.csv(mixtable_SRS3, "data/Mix_Quants/MT_SRS3.csv", row.names = FALSE)

# combinedSRS3 <- combine_sources(jags.SRS3, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton','Floc'), brown=c('Sawgrass, Periphyton')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS3$post, 2, median)
# # summary_stat(combinedSRS3, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS3_combined_sumstats" )


# SRS4 


SRS4mix <- SIa %>% filter(site == 'SRS4', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS4mix, "data/Consumers/SRS4mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/SRS4mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS4.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS4.csv"), mix)

plot_data(filename = "figures/isospace/SRS4_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/SRS4_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS4 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS4 <- run_model(run = "normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.SRS4 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS4/FCESRS4_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS4/FCESRS4_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS4/FCESRS4_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/SRS4/FCESRS4_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS4/FCESRS4_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.SRS4, mix, source, output_jags.SRS4)

mixtable_SRS4 = mixTable("data/JAGS_Output/SRS4/FCESRS4_sumstats.txt", type = "SRS4", nest = TRUE)

write.csv(mixtable_SRS4, "data/Mix_Quants/MT_SRS4.csv", row.names = FALSE)

# combinedSRS4 <- combine_sources(jags.SRS4, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton',"Epiphytes"), brown=c('Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS4$post, 2, median)
# # summary_stat(combinedSRS4, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS4_combined_sumstats" )

# SRS 6
SRS6mix <- SIa %>% filter(site == 'SRS6', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS6mix, "data/Consumers/SRS6mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/SRS6mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesSRS6.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_SRS6.csv"), mix)

plot_data(filename = "figures/isospace/SRS6_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/SRS6_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS6 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS6 <- run_model(run = "normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.SRS6 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS6/FCESRS6_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/SRS6/FCESRS6_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/SRS6/FCESRS6_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/SRS6/FCESRS6_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS6/FCESRS6_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.SRS6, mix, source, output_jags.SRS6)

mixtable_SRS6 = mixTable("data/JAGS_Output/SRS6/FCESRS6_sumstats.txt", type = "SRS6", nest = TRUE)

write.csv(mixtable_SRS6, "data/Mix_Quants/MT_SRS6.csv", row.names = FALSE)
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

write.csv(TS3mix, "data/Consumers/TS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/TS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS3.csv"), mix)

plot_data(filename = "figures/isospace/TS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS3_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS3 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS3 <- run_model(run = "normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.TS3 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS3/FCETS3_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS3/FCETS3_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS3/FCETS3_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS3/FCETS3_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS3/FCETS3_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

output_JAGS(jags.TS3, mix, source, output_jags.TS3)

mixtable_TS3 = mixTable("data/JAGS_Output/TS3/FCETS3_sumstats.txt", type = "TS3", nest = TRUE)

write.csv(mixtable_TS3, "data/Mix_Quants/MT_TS3.csv", row.names = FALSE)
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

write.csv(TS7mix, "data/Consumers/TS7mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/TS7mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS7.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS7.csv"), mix)

plot_data(filename = "figures/isospace/TS7_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS7_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS7 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS7 <- run_model(run = "normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS7 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS7/FCETS7_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS7/FCETS7_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS7/FCETS7_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS7/FCETS7_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS7/FCETS7_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

output_JAGS(jags.TS7, mix, source, output_jags.TS7)

mixtable_TS7 = mixTable("data/JAGS_Output/TS7/FCETS7_sumstats.txt", type = "TS7", nest = TRUE)

write.csv(mixtable_TS7, "data/Mix_Quants/MT_TS7.csv", row.names = FALSE)

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

write.csv(TS9mix, "data/Consumers/TS9mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/TS9mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS9.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS9.csv"), mix)

plot_data(filename = "figures/isospace/TS9_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS9_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS9 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS9 <- run_model(run = "normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.TS9 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS9/FCETS9_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS9/FCETS9_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS9/FCETS9_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS9/FCETS9_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS9/FCETS9_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

output_JAGS(jags.TS9, mix, source, output_jags.TS9)

mixtable_TS9 = mixTable("data/JAGS_Output/TS9/FCETS9_sumstats.txt", type = "TS9", nest = TRUE)

write.csv(mixtable_TS9, "data/Mix_Quants/MT_TS9.csv", row.names = FALSE)
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

write.csv(TS10mix, "data/Consumers/TS10mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/TS10mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS10.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS10.csv"), mix)

plot_data(filename = "figures/isospace/TS10_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS10_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS10 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS10 <- run_model(run = "normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.TS10 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS10/FCETS10_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = T,
                        plot_post_name = "data/JAGS_Output/TS10/FCETS10_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = T,
                        plot_pairs_name = "data/JAGS_Output/TS10/FCETS10_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = T,
                        plot_xy_name = "data/JAGS_Output/TS10/FCETS10_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS10/FCETS10_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)

output_JAGS(jags.TS10, mix, source, output_jags.TS10)

mixtable_TS10 = mixTable("data/JAGS_Output/TS10/FCETS10_sumstats.txt", type = "TS10", nest = TRUE)

write.csv(mixtable_TS10, "data/Mix_Quants/MT_TS10.csv", row.names = FALSE)
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

write.csv(TS11mix, "data/Consumers/TS11mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/Consumers/TS11mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/Sources/sourcesTS11.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/TEF/FCE_TEF_TS11.csv"), mix)

plot_data(filename = "figures/isospace/TS11_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/Consumers/TS11_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS11 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS11 <- run_model(run = "normal", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

output_jags.TS11 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/TS11/FCETS11_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "data/JAGS_Output/TS11/FCETS11_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = T,
                         plot_pairs_name = "data/JAGS_Output/TS11/FCETS11_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = T,
                         plot_xy_name = "data/JAGS_Output/TS11/FCETS11_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/TS11/FCETS11_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.TS11, mix, source, output_jags.TS11)

mixtable_TS11 = mixTable("data/JAGS_Output/TS11/FCETS11_sumstats.txt", type = "TS11", nest = TRUE)

write.csv(mixtable_TS11, "data/Mix_Quants/MT_TS11.csv", row.names = FALSE)
# ##combine posterior ground into brown/green
# combinedTS11 <- combine_sources(jags.TS11, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('SPOM','Epiphytes'), brown=c('Seagrass', 'Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedTS11$post, 2, median)
# summary_stat(combinedTS11, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
#              filename = "TS11_combined_sumstats" )


# Aggregating by Energy Channel ----
MixOut_RB10 = read.csv('data/Mix_Quants/MT_RB10.csv')
MixOut_SRS3 = read.csv('data/Mix_Quants/MT_SRS3.csv')
MixOut_SRS4 = read.csv('data/Mix_Quants/MT_SRS4.csv')
MixOut_SRS6 = read.csv('data/Mix_Quants/MT_SRS6.csv')
MixOut_TS3 = read.csv('data/Mix_Quants/MT_TS3.csv')
MixOut_TS7 = read.csv('data/Mix_Quants/MT_TS7.csv')
MixOut_TS9 = read.csv('data/Mix_Quants/MT_TS9.csv')
MixOut_TS10 = read.csv('data/Mix_Quants/MT_TS10.csv')
MixOut_TS11 = read.csv('data/Mix_Quants/MT_TS11.csv')

SRSMixout_gb <- rbind(MixOut_RB10, MixOut_SRS3, MixOut_SRS4, MixOut_SRS6)

SRSMixout_gb = SRSMixout_gb %>% 
  rename(site = type, season = code) %>%
  mutate(path = case_when(
    source %in% c("Epiphytes", "Phytoplankton", "Filamentous Green Algae", "Periphyton") ~ "green",
    source %in% c("Mangrove", "Sawgrass", "Red Macroalgae",'Floc') ~ "brown",
    TRUE ~ NA_character_  # For other cases, you can assign NA or something else if needed
  ))

SRSMixout_gb = SRSMixout_gb %>% 
  group_by(name, site, season, path) %>% 
  summarize(value = sum(mean)) %>% 
  pivot_wider(names_from = path, values_from = value)

TSMixout_gb <- rbind(MixOut_TS3, MixOut_TS7, MixOut_TS9, MixOut_TS10, MixOut_TS11)

TSMixout_gb = TSMixout_gb %>% 
  rename(site = type, season = code) %>%
  mutate(path = case_when(
    source %in% c("Epiphytes", "Phytoplankton", "Filamentous Green Algae", "Periphyton", "Epiphytic microalgae", 'SPOM') ~ "green",
    source %in% c("Mangrove", "Sawgrass", "Red Macroalgae", "Seagrass", 'Floc') ~ "brown",
    TRUE ~ NA_character_  # For other cases, you can assign NA or something else if needed
  ))

TSMixout_gb = TSMixout_gb %>% 
  group_by(name, site, season, path) %>% 
  summarize(value = sum(mean)) %>% 
  pivot_wider(names_from = path, values_from = value)


y_label_formatter <- function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}


combined_df = bind_rows(SRSMixout_gb,TSMixout_gb) %>%  
  mutate(transect = case_when(
    site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
    site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"))

unique_names <- combined_df %>%
  group_by(site, season) %>%
  summarize(unique_names = n_distinct(name))


# Combined Brown vs Green boxplot----
combined_df <- combined_df %>%
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

combined_df = combined_df %>% 
  group_by(site, season) %>% 
  mutate(fill = mean(green),
         site = factor(site, levels = c( "Outer Bay","Mid Bay","Inner Bay","Mangrove Ecotone" ,"TS Marsh", "Lower River", "Mid River","Upper River" ,"SRS Marsh" )))

# generates a table of average contributions for brown vs green at each site during each season
bvg_table = combined_df %>% 
  group_by(site,season) %>% 
  summarise(green = mean(green))


bvg_table_ft <- flextable(bvg_table,
                           col_keys = c("site", "season","green")) %>%
  add_header_row(colwidths = c(1,1,1), values = c("Site", "Season", "Mean Green Pathway Contribution")) %>%
  set_header_labels(site = "Site",
                    season = "Season",
                    green = "Mean Green Pathway Contribution") %>%
  colformat_double(digits = 2) %>%
  theme_box() %>%
  align(align = "center") %>%
  align(part = "header", align = "center") %>% 
  # compose(j = "Genus_spp",
  #         value = as_paragraph(as_i(Genus_spp))) %>%
  merge_v(part = "header")

bvg_table_ft

save_as_docx(bvg_table_ft, path = "tables/bvg_table_flex.docx")


mixoutput_bxplt_gb_combined <-ggplot(combined_df,aes(x=site, y = green, fill=fill, width=0.8))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_gradient2(low = "#663300",
                       high = "#92D050",
                       mid = 'white',
                       midpoint = 0.5,
                       limits = c(0,1),
                       na.value = "grey50") +
  facet_grid(transect~season, scales = "free_y") +
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black"), 
        plot.title = element_text(size = 18, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0,1),
    labels = y_label_formatter
  ) +
  labs(y="Green Pathway Source Contribution",
       x=NULL,
       fill = "Green Pathway\nSource Contribution\n ")+
  coord_flip()

mixoutput_bxplt_gb_combined

ggsave("figures/mixoutput_bxplt_gb.png", width = 11, height = 6, dpi = 600)

# Source Contribution Plots ----

# combined dataset from mixing model outputs
cont_df = rbind(MixOut_RB10, MixOut_SRS3, MixOut_SRS4, MixOut_SRS6, MixOut_TS3, MixOut_TS7, MixOut_TS9, MixOut_TS10, MixOut_TS11)

# renaming for plotting
cont_df <- cont_df %>%
  mutate(source = ifelse(source == "Phytoplankton", "Phyto.",
                         ifelse(source == "Mangrove", "Mang.",
                                ifelse(source == "Sawgrass", "Sawgrass",
                                       ifelse(source == "Floc", "Floc",
                                              ifelse(source == "Red Macroalgae", "RMA",
                                                     ifelse(source == "Seagrass", "Seagrass",
                                                            ifelse(source == "SPOM", "POM",
                                                                   ifelse(source == "Epiphytes", "Epi.",
                                                                          ifelse(source == "Periphyton", "Peri.",
                                                                                 ifelse(source == "Filamentous Green Algae", "FGA", source)
                                                                          )
                                                                   )
                                                            )
                                                     )
                                              )
                                       )
                                )
                         )))

cont_df = cont_df %>% 
  rename(site = type, season = code)

cont_df = cont_df %>% 
  group_by(site, season, source, name) %>% 
  summarize(value = sum(mean)) 

cont_df = cont_df %>% 
mutate(transect = case_when(
  site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
  site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"),
  source = factor(source, levels = c( "Mang.","Sawgrass","Floc","RMA", "Seagrass", "Epi.","Peri." ,"Phyto.",  "FGA", "POM")))


text_colors <- c("Mang." = "saddlebrown", "Sawgrass" = "saddlebrown", "Floc" = "saddlebrown","RMA" = "saddlebrown","Seagrass" = "saddlebrown",
                 "Epi." = "forestgreen", "Peri." = "forestgreen", "Phyto." = "forestgreen",
                 "FGA" = "forestgreen", "POM"= "forestgreen" )

cont_df <- cont_df %>%
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

source_cont_plot = ggplot(cont_df, aes(x = source, y = value, fill = season, width = 0.2)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site, scales = "free_x") +
  theme(
    axis.title = element_text(size = 32), 
    axis.text.y = element_text(size = 20, colour = "black", face = "bold", family = "arial"), 
    axis.text.x = element_text(
      size = 18,
      colour = c("black", "black", "black", "black")), 
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'top',
    legend.title = element_text(size = 32, face = "bold", family = "arial"),
    strip.text.x = element_text(size = 26, face = "bold", family = "arial"),
    strip.text.y = element_text(size = 18, face = "bold", family = "arial"),
    legend.text = element_text(size = 24, face = "bold", family = "arial")) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season")

source_cont_plot

ggsave("figures/source_cont_plot.png", width = 15, height = 11, dpi = 600)

# creates a table of mean values for source specific contributions. 
cont_table = cont_df %>% 
  group_by(site,season, source) %>% 
  summarise(value = mean(value))


cont_table_ft <- flextable(cont_table,
                       col_keys = c("site", "season","source",
                                    "value")) %>%
  add_header_row(colwidths = c(1,1,1,1), values = c("Site", "Season", "Source","Total Energy Contribution")) %>%
  set_header_labels(site = "Site",
                    season = "Season",
                    source = "Source",
                    value   = "Total Energy Contribution") %>%
  colformat_double(digits = 2) %>%
  theme_box() %>%
  align(align = "center") %>%
  align(part = "header", align = "center") %>% 
  # compose(j = "Genus_spp",
  #         value = as_paragraph(as_i(Genus_spp))) %>%
  merge_v(part = "header")

cont_table_ft

save_as_docx(cont_table_ft, path = "tables/cont_table_flex.docx")
write.csv(cont_table,"data/Sources/cont_table.csv",row.names = F) 


# SRS only source contribution plot
cont_SRS = cont_df %>% 
  filter(transect == "Shark River Slough")
  # group_by(site,season, source) %>% 
  # summarise(value = mean(value))

source_plot_SRS = ggplot(cont_SRS, aes(x = source, y = value, fill = season, width = 0.2)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site, scales = "free_x", nrow = 1) +
  theme(
    axis.ticks.length.x = unit(0.1, "inch"),
    axis.title = element_text(size = 32), 
    axis.text.y = element_text(size = 20, colour = "black", face = "bold", family = "arial"), 
    axis.text.x = element_text(
      size = 22,
      colour = c("black", "black", "black", "black", face = "bold", family = "arial") 
    ), 
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'top',
    legend.title = element_text(size = 32, face = "bold", family = "arial"),
    strip.text.x = element_text(size = 32, face = "bold", family = "arial"),
    strip.text.y = element_text(size = 18, face = "bold", family = "arial"),
    legend.text = element_text(size = 24, face = "bold", family = "arial")
  ) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter
  ) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season"
  ) 


source_plot_SRS

ggsave("figures/source_plot_SRS.png", width = 22, height = 10, dpi = 600)

# TS only source contribution plot

cont_TS = cont_df %>% 
  filter(transect == "Taylor Slough")

cont_TS = cont_TS %>% 
  mutate(site = factor(site, levels = c( "TS Marsh","Mangrove Ecotone", "Inner Bay","Mid Bay","Outer Bay")))


source_plot_TS <- ggplot(cont_TS, aes(x = source, y = value, fill = season, width = 0.2)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4) +
  geom_rect(aes(xmin=source, xmax=source, ymin=-Inf, ymax=Inf, fill = season), alpha=0.5, stat="identity") +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ site, scales = "free_x", nrow = 1) +
  theme(
    axis.title = element_text(size = 32), 
    axis.text.y = element_text(size = 20, colour = "black", face = "bold", family = "arial"), 
    axis.text.x = element_text(
      size = 18,
      colour = c("black", "black", "black", "black")
    ),
    plot.title = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'top',
    legend.title = element_text(size = 32, face = "bold", family = "arial"),
    strip.text.x = element_text(size = 28, face = "bold", family = "arial"),
    strip.text.y = element_text(size = 18, face = "bold", family = "arial"),
    legend.text = element_text(size = 24, face = "bold", family = "arial")
  ) +
  scale_fill_manual(values = c("Dry" = "black", "Wet" = "white")) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter
  ) +
  labs(
    y = "Basal Resource Energy Contribution (%)",
    x = NULL,
    fill = "Season"
  )

source_plot_TS


ggsave("figures/source_plot_TS.png", width = 22, height = 10, dpi = 600)

# FCE proposal plots ----
combined_df_dry = combined_df %>% 
  filter(season == "Dry")

dry_plot <-ggplot(combined_df_dry,aes(x=site, y = green, fill=fill, width=0.8))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_gradient2(low = "saddlebrown",
                       high = "forestgreen",
                       mid = 'white',
                       midpoint = 0.5,
                       limits = c(0,1),
                       na.value = "grey50") +
  facet_grid(transect~season, scales = "free_y") +
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black"), 
        plot.title = element_text(size = 18, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0,1),
    labels = y_label_formatter
  ) +
  labs(y="Green Pathway Source Contribution",
       x=NULL,
       fill = "Green Pathway\nSource Contribution\n ")+
  coord_flip()

dry_plot
ggsave("figures/dry_plot.png", width = 8, height = 6, dpi = 600)




no_season <-ggplot(combined_df,aes(x=site, y = green, fill=fill, width=0.8))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_gradient2(low = "saddlebrown",
                       high = "forestgreen",
                       mid = 'white',
                       midpoint = 0.5,
                       limits = c(0,1),
                       na.value = "grey50") +
  facet_grid(transect~season, scales = "free_y") +
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black"), 
        plot.title = element_text(size = 18, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0,1),
    labels = y_label_formatter
  ) +
  labs(y="Green Pathway Source Contribution",
       x=NULL,
       fill = "Green Pathway\nSource Contribution\n ")+
  coord_flip()

no_season

ggsave("figures/no_season.png", width = 8, height = 6, dpi = 600)


allmix = SRSMixout_gb <- rbind(MixOut_RB10, MixOut_SRS3, MixOut_SRS4, MixOut_SRS6,MixOut_TS3, MixOut_TS7, MixOut_TS9, MixOut_TS10, MixOut_TS11)


write.csv(allmix,"data/Consumers/allmix.csv",row.names = F) 
write.csv(combined_df,"data/Consumers/combined.csv",row.names = F)


  





  





stacked_source_cont = ggplot(cont_table, aes(x = season, y = value, fill = source)) +
  geom_col(colour = "black", position = "stack", width = 0.5) +  # Set width to adjust spacing
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  facet_wrap(facets = c("site"), nrow = 3, ncol = 3) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 24),
    strip.text.x = element_text(size = 24),
    strip.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)) +
  scale_x_discrete(expand = c(0,1)) +
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter
  ) +
   scale_fill_manual(values = c("Mang." = "wheat4", "Sawgrass" = "wheat4", "Floc" = "wheat1","RMA" = "wheat1","Seagrass" = "wheat1",
                               "Epi." = "darkolivegreen4", "Peri." = "darkolivegreen4", "Phyto." = "darkolivegreen1",
                               "FGA" = "darkolivegreen1", "POM"= "darkolivegreen1" )) +
  labs(
    y = "Proportional Energy Contribution (%)",
    x = NULL,
    fill = "Basal Resource"
  )

stacked_source_cont

ggsave("figures/stacked_source_cont.png", width = 12, height = 11, dpi = 600)

# END
