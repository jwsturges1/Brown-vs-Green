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
SRS_sumstats_gb<-read.csv('data/SRSMixout_gb.csv') %>% 
  mutate(transect = "Shark River Slough", 
         season = if_else(season == "wet19", "Wet Season", "Dry Season"),
         site = fct_relevel(site, "SRS6","SRS4","RB10", "SRS3"))

TS_sumstats_gb<-read.csv('data/TS_sumstats_gb.csv') %>% 
  mutate(transect = "Taylor Slough", 
         season = if_else(season == "wet19", "Wet Season", "Dry Season"),
         site = fct_relevel(site, "TS11","TS10","TS9","TS7","TS3"))


y_label_formatter <- function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

combine_df = SRS_sumstats_gb %>% 
  bind_rows(TS_sumstats_gb) %>% 
  filter(source == "green")


#### SRS ----

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
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "top")+
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()

mixoutput_bxplt_gb_SRS

ggsave("figures/mixoutput_bxplt_gb_SRS.png", width = 10, height = 8, dpi = 300)


SRS_sumstats_gb$season <- gsub("wet19", "Wet 2019", SRS_sumstats_gb$season, ignore.case = TRUE)
SRS_sumstats_gb$season <- gsub("dry19", "Dry 2019", SRS_sumstats_gb$season, ignore.case = TRUE)

#### TS ----

TS_sumstats_gb<-read.csv('data/TS_sumstats_gb.csv')
TS_sumstats_gb$site<-fct_relevel(TS_sumstats_gb$site, "TS3","TS7","TS9","TS10","TS11")

mixoutput_bxplt_gb_TS<-ggplot(TS_sumstats_gb,aes(x=source, fill=source, width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
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
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = y_label_formatter,
    expand = c(0.01, 0)
  ) +
  labs(y="Proportional Dietary Contribution")+
  coord_flip()


mixoutput_bxplt_gb_TS
ggsave("figures/mixoutput_bxplt_gb_TS.png", width = 10, height = 8, dpi = 300)



#### combined plot ----



mixoutput_bxplt_gb_combined <-ggplot(combine_df,aes(x=site, fill=X50., width=0.8))+
  geom_boxplot(aes(lower = X25., upper = X75., middle = X50., ymin = X2.50., ymax = X97.50.), stat="identity")+
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

mixoutput_bxplt_gb_combined

ggsave("figures/mixoutput_bxplt_gb_combined.png", width = 9, height = 6, dpi = 600)

TS3mix<-SIa %>% filter(site == 'TS3', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


#write.csv(TS3mix,"TS3mix.csv",row.names = F) 

mix <- load_mix_data(filename="data/TS3mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/sourcesTS3.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TS3.csv"), mix)

plot_data(filename="TS3_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "TS3_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.TS3 <- run_model(run="test", mix, source, discr, model_filename, 
                      alpha.prior = 1, resid_err=F, process_err=T)


jags.TS3 <- run_model(run="normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err=F, process_err=T)

output_jags.TS3  <- list(summary_save = TRUE,
                         summary_name = "FCETS3_sumstats",
                         sup_post = FALSE,
                         plot_post_save_pdf = T,
                         plot_post_name = "FCETS3_plot",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = F,
                         plot_pairs_name = "FCETS3_pairs",
                         sup_xy = T,
                         plot_xy_save_pdf = F,
                         plot_xy_name = "FCETS3_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = T,
                         diag_save = T,
                         diag_name = "FCETS3_Diagnostic",
                         indiv_effect = FALSE,
                         plot_post_save_png = T,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = F)

output_JAGS(jags.TS3 , mix, source, output_jags.TS3)






mixTable("FCETS3_sumstats.txt",type = "TS3", nest = T)



##combine posterior ground into brown/green

combinedTS3 <- combine_sources(jags.TS3, mix, source, alpha.prior=1, 
                               groups=list(green=c('Periphyton'), brown=c('Plant', 'Floc' )))

# get posterior medians for new source groupings
apply(combinedTS3$post, 2, median)
summary_stat(combinedTS3, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "TS3_combined_sumstats" )


#### Species Level ----

mixTable = function(file,type,ind = F,nest = F, csv = F){
  require(tidyverse)
  cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '75%', '95%', '97.5%')
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
      x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
      
      df = tibble(x$name, x$type, x$source, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                  x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
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
      x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
      
      df = tibble(x$name, x$type, x$source, x$code, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                  x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
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

