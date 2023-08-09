# Authors: James Sturges, Ryan Rezek, Ryan James
# Last Updated 19 July 2023
# Seasonal differences in brown/green energy pathways across 9 FCE sites
# DATASETUP ----

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
SI <- read.csv("data/FCE_SI_data_xls_08_August_2023_Sturges_edits.csv",na.strings=c("","NA"))%>%filter(is.na(outlier))

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


# RB10 Mixing Model ----
RB10mix<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")

write.csv(RB10mix,"data/RB10mix.csv",row.names = F) 



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

discr <- load_discr_data(file("data/FCE_TEF_RB10.csv"), mix)

plot_data(filename="figures/isospace/RB10_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "data/RB10_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.RB10 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=F)


jags.RB10 <- run_model(run="very long", mix, source, discr, model_filename,
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

# # combines sources into energy channel groups (brown or green pathway)
# combinedRB10 <- combine_sources(jags.RB10, mix, source, alpha.prior=1,
#                                 groups=list(green=c('Phytoplankton','Epiphytes'),
#                                             brown=c('Wet Mangrove', 'Dry Mangrove')))
# 
# 
# # get posterior medians for new source groupings
# apply(combinedRB10$post, 2, median)
# # summary_stat(combinedRB10, meanSD=T, quantiles=c(c(0.025, 0.25, 0.5, 0.75, 0.975)), savetxt=T,
# #              filename = "RB10_combined_sumstats" )



# Shark River Slough Mixing Models 3,4,6 ----

# SRS3 

SRS3mix <- SIa %>% filter(site == 'SRS3', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS3mix, "data/SRS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/SRS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesSRS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS3.csv"), mix)

plot_data(filename = "figures/isospace/SRS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/SRS3_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS3 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS3 <- run_model(run = "very long", mix, source, discr, model_filename,
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

write.csv(mixtable_SRS3, "data/Mix_Quants/SRS3/MT_SRS3.csv", row.names = FALSE)

# combinedSRS3 <- combine_sources(jags.SRS3, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton','Floc'), brown=c('Sawgrass, Periphyton')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS3$post, 2, median)
# # summary_stat(combinedSRS3, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS3_combined_sumstats" )


# SRS4 


SRS4mix <- SIa %>% filter(site == 'SRS4', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS4mix, "data/SRS4mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/SRS4mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesSRS4.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS4.csv"), mix)

plot_data(filename = "figures/isospace/SRS4_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/SRS4_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS4 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS4 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.SRS4 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS4/FCERSRS4_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "data/JAGS_Output/SRS4/FCERSRS4_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = FALSE,
                         plot_pairs_name = "data/JAGS_Output/SRS4/FCERSRS4_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "data/JAGS_Output/SRS4/FCERSRS4_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS4/FCERSRS4_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)
output_JAGS(jags.SRS4, mix, source, output_jags.SRS4)

mixtable_SRS4 = mixTable("data/JAGS_Output/SRS4/FCERSRS4_sumstats.txt", type = "SRS4", nest = TRUE)

# combinedSRS4 <- combine_sources(jags.SRS4, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton',"Epiphytes"), brown=c('Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS4$post, 2, median)
# # summary_stat(combinedSRS4, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS4_combined_sumstats" )

# SRS 6
SRS6mix <- SIa %>% filter(site == 'SRS6', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(SRS6mix, "data/SRS6mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/SRS6mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesSRS6.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS6.csv"), mix)

plot_data(filename = "figures/isospace/SRS6_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/SRS6_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.SRS6 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.SRS6 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.SRS6 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/SRS6/FCERSRS6_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "data/JAGS_Output/SRS6/FCERSRS6_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = FALSE,
                         plot_pairs_name = "data/JAGS_Output/SRS6/FCERSRS6_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "data/JAGS_Output/SRS6/FCERSRS6_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/SRS6/FCERSRS6_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)
output_JAGS(jags.SRS6, mix, source, output_jags.SRS6)

mixtable_SRS6 = mixTable("data/JAGS_Output/SRS6/FCERSRS6_sumstats.txt", type = "SRS6", nest = TRUE)


# combinedSRS6 <- combine_sources(jags.SRS6, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('Phytoplankton','Filamentous Green Algae' ), brown=c('Mangrove', 'Red Macroalgae')))
# 
# # get posterior medians for new source groupings
# apply(combinedSRS6$post, 2, median)
# # summary_stat(combinedSRS6, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
# #              filename = "SRS6_combined_sumstats" )

# Taylor Slough Transect Mixing Models ----

# TS3

TS3mix <- SIa %>% filter(site == 'TS3', common_name != "Egyptian paspalidium", group == 'Consumer') %>% rename('d13C' = 'md13C', 'd15N' = 'md15N', 'd34S' = 'md34S')

write.csv(TS3mix, "data/TS3mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/TS3mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesTS3.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS3.csv"), mix)

plot_data(filename = "figures/isospace/TS3_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/TS3_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS3 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS3 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS3 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS3/FCERTS3_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = FALSE,
                        plot_post_name = "data/JAGS_Output/TS3/FCERTS3_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = FALSE,
                        plot_pairs_name = "data/JAGS_Output/TS3/FCERTS3_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = FALSE,
                        plot_xy_name = "data/JAGS_Output/TS3/FCERTS3_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS3/FCERTS3_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)
output_JAGS(jags.TS3, mix, source, output_jags.TS3)
mixtable_TS3 = mixTable("data/JAGS_Output/TS3/FCERTS3_sumstats.txt", type = "TS3", nest = TRUE)

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

write.csv(TS7mix, "data/TS7mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/TS7mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesTS7.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS7.csv"), mix)

plot_data(filename = "figures/isospace/TS7_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/TS7_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS7 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS7 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS7 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS7/FCERTS7_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = FALSE,
                        plot_post_name = "data/JAGS_Output/TS7/FCERTS7_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = FALSE,
                        plot_pairs_name = "data/JAGS_Output/TS7/FCERTS7_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = FALSE,
                        plot_xy_name = "data/JAGS_Output/TS7/FCERTS7_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS7/FCERTS7_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)
output_JAGS(jags.TS7, mix, source, output_jags.TS7)
mixtable_TS7 = mixTable("data/JAGS_Output/TS7/FCERTS7_sumstats.txt", type = "TS7", nest = TRUE)

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

write.csv(TS9mix, "data/TS9mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/TS9mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesTS9.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS9.csv"), mix)

plot_data(filename = "figures/isospace/TS9_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/TS9_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS9 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS9 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS9 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS9/FCERTS9_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = FALSE,
                        plot_post_name = "data/JAGS_Output/TS9/FCERTS9_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = FALSE,
                        plot_pairs_name = "data/JAGS_Output/TS9/FCERTS9_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = FALSE,
                        plot_xy_name = "data/JAGS_Output/TS9/FCERTS9_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS9/FCERTS9_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)
output_JAGS(jags.TS9, mix, source, output_jags.TS9)
mixtable_TS9 = mixTable("data/JAGS_Output/TS9/FCERTS9_sumstats.txt", type = "TS9", nest = TRUE)

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

write.csv(TS10mix, "data/TS10mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/TS10mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesTS10.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS10.csv"), mix)

plot_data(filename = "figures/isospace/TS10_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/TS10_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS10 <- run_model(run = "test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS10 <- run_model(run = "very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS10 <- list(summary_save = TRUE,
                        summary_name = "data/JAGS_Output/TS10/FCERTS10_sumstats",
                        sup_post = FALSE,
                        plot_post_save_pdf = FALSE,
                        plot_post_name = "data/JAGS_Output/TS10/FCERTS10_plot",
                        sup_pairs = FALSE,
                        plot_pairs_save_pdf = FALSE,
                        plot_pairs_name = "data/JAGS_Output/TS10/FCERTS10_pairs",
                        sup_xy = TRUE,
                        plot_xy_save_pdf = FALSE,
                        plot_xy_name = "data/JAGS_Output/TS10/FCERTS10_plot",
                        gelman = TRUE,
                        heidel = FALSE,
                        geweke = TRUE,
                        diag_save = TRUE,
                        diag_name = "data/JAGS_Output/TS10/FCERTS10_Diagnostic",
                        indiv_effect = FALSE,
                        plot_post_save_png = FALSE,
                        plot_pairs_save_png = FALSE,
                        plot_xy_save_png = FALSE)
output_JAGS(jags.TS10, mix, source, output_jags.TS10)

mixtable_TS10 = mixTable("data/JAGS_Output/TS10/FCERTS10_sumstats.txt", type = "TS10", nest = TRUE)


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

write.csv(TS11mix, "data/TS11mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "data/TS11mix.csv",
                     iso_names = c("d13C", "d15N", "d34S"),
                     factors = c('common_name', 'hydroseason'),
                     fac_random = c(FALSE, FALSE),
                     fac_nested = c(FALSE, FALSE),
                     cont_effects = NULL)

source <- load_source_data(filename = "data/sourcesTS11.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS11.csv"), mix)

plot_data(filename = "figures/isospace/TS11_isospace_plot",
          plot_save_pdf = TRUE,
          plot_save_png = TRUE,
          mix, source, discr)

model_filename <- "data/TS11_mix.txt"
write_JAGS_model(model_filename, resid_err = FALSE, process_err = TRUE, mix, source)

jags.TS11 <- run_model(run = "test", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)

jags.TS11 <- run_model(run = "very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err = FALSE, process_err = FALSE)
output_jags.TS11 <- list(summary_save = TRUE,
                         summary_name = "data/JAGS_Output/TS11/FCERTS11_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "data/JAGS_Output/TS11/FCERTS11_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = FALSE,
                         plot_pairs_name = "data/JAGS_Output/TS11/FCERTS11_pairs",
                         sup_xy = TRUE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "data/JAGS_Output/TS11/FCERTS11_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "data/JAGS_Output/TS11/FCERTS11_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)
output_JAGS(jags.TS11, mix, source, output_jags.TS11)

mixtable_TS11 = mixTable("data/JAGS_Output/TS11/FCERTS11_sumstats.txt", type = "TS11", nest = TRUE)

# ##combine posterior ground into brown/green
# combinedTS11 <- combine_sources(jags.TS11, mix, source, alpha.prior=1, 
#                                 groups=list(green=c('SPOM','Epiphytes'), brown=c('Seagrass', 'Mangrove')))
# 
# # get posterior medians for new source groupings
# apply(combinedTS11$post, 2, median)
# summary_stat(combinedTS11, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T, 
#              filename = "TS11_combined_sumstats" )


# SRS boxplot ----

mixoutput_bxplt_gb_SRS<-ggplot(SRS_sumstats_gb,aes(x=season, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X105., middle = X50., ymin = X2.50., ymax = X910.50.), stat="identity")+
  theme_bw()+
  scale_fill_manual(values=c("saddlebrown",'limegreen'))+ 
  coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "top")+
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.105, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()

mixoutput_bxplt_gb_SRS

ggsave("figures/mixoutput_bxplt_gb_SRS.png", width = 10, height = 8, dpi = 300)


SRS_sumstats_gb$season <- gsub("wet19", "Wet 2019", SRS_sumstats_gb$season, ignore.case = TRUE)
SRS_sumstats_gb$season <- gsub("dry19", "Dry 2019", SRS_sumstats_gb$season, ignore.case = TRUE)

# TS boxplot ----

TS_sumstats_gb<-read.csv('data/TS_sumstats_gb.csv')
TS_sumstats_gb$site<-fct_relevel(TS_sumstats_gb$site, "TS3","TS10","TS9","TS10","TS11")

mixoutput_bxplt_gb_TS<-ggplot(TS_sumstats_gb,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X105., middle = X50., ymin = X2.50., ymax = X910.50.), stat="identity")+
  theme_bw()+
  scale_fill_manual(values=c("saddlebrown", "limegreen"))+ 
  coord_cartesian(ylim = c( 0,1))+facet_grid(site~season)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 12))+
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.105, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()


mixoutput_bxplt_gb_TS
ggsave("figures/mixoutput_bxplt_gb_TS.png", width = 10, height = 8, dpi = 300)



# Combined boxplot ----



mixoutput_bxplt_gb_combined <-ggplot(combine_df,aes(x=site, fill=X50., width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X105., middle = X50., ymin = X2.50., ymax = X910.50.), stat="identity")+
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
    breaks = c(0.0, 0.25, 0.5, 0.105, 1.0),
    limits = c(0,1),
    labels = y_label_formatter
  ) +
  labs(y="Green Pathway Source Contribution",
       x=NULL,
       fill = "Green Pathway\nSource Contribution\n ")+
  coord_flip()

mixoutput_bxplt_gb_combined

ggsave("figures/mixoutput_bxplt_gb_combined.png", width = 9, height = 6, dpi = 600)


