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



#### DATASETUP ----

# Create mixTable function from Ryan James -
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

# processed stable isotope values used in these analyses
SI <- read.csv("data/SI_data_csv_8_30_21.csv",na.strings=c("","NA"))%>%filter(is.na(outlier))

summary(SI)

SIa<-SI %>% group_by(site,season,group, common_name, functional_grp) %>% 
  summarise(n=n(),md13C = mean(as.numeric(d13C), na.rm = T),d13Csd=sd(as.numeric(d13C),na.rm = T),
            md15N=mean(as.numeric(d15N),na.rm = T),d15Nsd=sd(as.numeric(d15N),na.rm = T),
            md34S=mean(as.numeric(d34S),na.rm = T),d34Ssd=sd(d34S,na.rm = T)) %>%
  filter(!is.na(md13C),!is.na(md15N),!is.na(md34S))


SIb<-SI %>% group_by(site,group, common_name, functional_grp) %>% 
  summarise(n=n(),md13C = mean(d13C,na.rm = T),d13Csd=sd(d13C,na.rm = T), md15N=mean(d15N,na.rm = T),d15Nsd=sd(d15N,na.rm = T),md34S=mean(d34S,na.rm = T),d34Ssd=sd(d34S,na.rm = T))  %>%
  filter(!is.na(md13C),!is.na(md15N),!is.na(md34S))

# SRS3 Bi-plots 
SRS3<-SIa %>% filter(site == 'SRS3', common_name!="Egyptian paspalidium")


SRS3CSbiplot<-ggplot(SRS3,aes(x = md13C,y = md34S )) +            #coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS3,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS3,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS3, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=3,check_overlap = T)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


SRS3CNbiplot<-ggplot(SRS3,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS3,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS3,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS3, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name), size =3, check_overlap = T)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)

#TS3 Bi-plots 
TS3<-SIa %>% filter(site == 'TS3', common_name!="Egyptian paspalidium")


TS3CSbiplot<-ggplot(TS3,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS3,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS3,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS3, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


TS3CNbiplot<-ggplot(TS3,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS3,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS3,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS3, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()



#RB10 Bi-plots 
RB10<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium")


RB10CSbiplot<-ggplot(RB10,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = RB10,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = RB10,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=RB10, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = T)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


RB10CNbiplot<-ggplot(RB10,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = RB10,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = RB10,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=RB10, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name),check_overlap = T)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#SRS4 Bi-plots 

SRS4<-SIa %>% filter(site == 'SRS4', common_name!="Egyptian paspalidium")


SRS4CSbiplot<-ggplot(SRS4,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS4,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS4,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS4, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


SRS4CNbiplot<-ggplot(SRS4,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS4,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS4,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS4, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#SRS6 Bi-plots 
SRS6<-SIa %>% filter(site == 'SRS6', common_name!="Egyptian paspalidium")


SRS6CSbiplot<-ggplot(SRS6,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS6,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS6,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS6, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


SRS6CNbiplot<-ggplot(SRS6,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = SRS6,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = SRS6,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=SRS6, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#TS7 Bi-plots 
TS7<-SIa %>% filter(site == 'TS7', common_name!="Egyptian paspalidium")


TS7CSbiplot<-ggplot(TS7,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS7,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS7,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS7, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


TS7CNbiplot<-ggplot(TS7,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS7,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS7,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS7, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#TS9 Bi-plots 
TS9<-SIa %>% filter(site == 'TS9', common_name!="Egyptian paspalidium")


TS9CSbiplot<-ggplot(TS9,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS9,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS9,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS9, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


TS9CNbiplot<-ggplot(TS9,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS9,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS9,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS9, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#TS10 Bi-plots 
TS10<-SIa %>% filter(site == 'TS10')


TS10CSbiplot<-ggplot(TS10,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS10,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS10,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS10, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = T)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


TS10CNbiplot<-ggplot(TS10,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS10,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS10,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS10, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()


#TS11 Bi-plots 
TS11<-SIa %>% filter(site == 'TS11')


TS11CSbiplot<-ggplot(TS11,aes(x = md13C,y = md34S )) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS11,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS11,aes(ymin = md34S - d34Ssd,ymax = md34S + d34Ssd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS11, aes(color=functional_grp),size=2)+ geom_text(aes(label = common_name),size=4,check_overlap = F)+#scale_shape_manual(values=c(1,0,15,16,2))+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + facet_wrap(~season,nrow=2)+scale_color_viridis_d()


TS11CNbiplot<-ggplot(TS11,aes(x = md13C,y = md15N)) +#coord_fixed(ratio = 1)+
  geom_errorbarh(data = TS11,aes(xmin = md13C - d13Csd ,xmax = md13C + d13Csd), height=0 ,color="#999999") + 
  geom_errorbar(data = TS11,aes(ymin = md15N - d15Nsd,ymax = md15N + d15Nsd),width=0,color="#999999")+ theme_bw()+
  geom_point(data=TS11, aes(color=functional_grp),size=2)+scale_shape_manual(values=c(1,0,15,16,2))+
  geom_text(aes(label = common_name))+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+theme(text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ facet_wrap(~season,nrow=2)+scale_color_viridis_d()

SIa$site<-fct_relevel(SIa$site, "SRS3","RB10", "SRS4","SRS6","TS3","TS7","TS9","TS10","TS11")


SIa1<-SIa%>% filter(group != 'Source')
transectCS<-ggplot(SIa1,aes(x = md13C,y = md34S ))+
  theme_bw()+
  geom_point(data=SIa1, aes(color=group),size=2)+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  facet_wrap(season~site,nrow=2)+scale_color_viridis_d()




transectCN<-ggplot(SIa1,aes(x = md13C,y = md15N ))+
  theme_bw()+
  geom_point(data=SIa1, aes(color=group),size=2)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(text = element_text(size=14))+scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  facet_wrap(season~site,nrow=2)+scale_color_viridis_d()


#Mixing Models (Early version for Brownbag)----


#TS3
TS3mix<-SIa %>% filter(site == 'TS3', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


#write.csv(TS3mix,"TS3mix.csv",row.names = F) 

# load consumer data
mix <- load_mix_data(filename="data/TS3mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

# load source data
source <- load_source_data(filename="data/sourcesTS3.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

# load discrimination factors
discr <- load_discr_data(file("data/FCE_TEF_TS3.csv"), mix)

#generates bioplots for CN, CS, and NS in isospace
plot_data(filename="TS3_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

#writing model from mix files, residual error is set to false even though we have single replicate species representatives most of these are composite samples composed of multiple individual organisms
model_filename <- "TS3_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.TS3 <- run_model(run="test", mix, source, discr, model_filename, 
                      alpha.prior = 1, resid_err=F, process_err=T)

#running all models on very long
jags.TS3 <- run_model(run="very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err=F, process_err=T)

#generate output summary stats, plots, and model diagnostics
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

#out of TS3 mixing model
output_JAGS(jags.TS3 , mix, source, output_jags.TS3)

# use mixTable function
mixTable("FCETS3_sumstats.txt",type = "TS3", nest = T)

##combine posterior ground into brown/green
combinedTS3 <- combine_sources(jags.TS3, mix, source, alpha.prior=1, 
                               groups=list(green=c('Periphyton'), brown=c('Plant', 'Floc' )))

# get posterior medians for new source groupings
apply(combinedTS3$post, 2, median)
summary_stat(combinedTS3, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "TS3_combined_sumstats" )


#### Taylor Slough (Florida Bay) Marine Transect TS 7, 9, 10, & 11----

TSFBmix<-SIa %>% filter(site == 'TS7'|site =='TS9'|site =='TS10'|site =='TS11', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


write.csv(TSFBmix,"TSFBmix.csv",row.names = F) 

mix <- load_mix_data(filename="TSFBmix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)


source <- load_source_data(filename="data/FCE_sources_TSFB.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_TSFB.csv"), mix)

plot_data(filename="TSFB_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)


model_filename <- "TSFB_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.TSFB <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.TSFB <- run_model(run="very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)

output_jags.TSFB  <- list(summary_save = TRUE,
                          summary_name = "FCETSFB_sumstats",
                          sup_post = FALSE,
                          plot_post_save_pdf = T,
                          plot_post_name = "FCETSFB_plot",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "FCETSFB_pairs",
                          sup_xy = T,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "FCETSFB_plot",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = T,
                          diag_save = F,
                          diag_name = "FCETSFB_Diagnostic",
                          indiv_effect = FALSE,
                          plot_post_save_png = T,
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = F)

output_JAGS(jags.TSFB , mix, source, output_jags.TSFB)

##combine posterior ground into brown/green

combinedTSFB <- combine_sources(jags.TSFB, mix, source, alpha.prior=1, 
                                groups=list(green=c('Epiphytic microalgae','SPOM' ), brown=c('Seagrass', 'Red Mangrove')))

# get posterior medians for new source groupings
apply(combinedTSFB$post, 2, median)
summary_stat(combinedTSFB, meanSD=FALSE, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "TSFB_combined_sumstats" )







#### Shark River Slough Transects SRS 3,4,6 ----

### SRS 6 ---
SRS6mix <- SIa %>% filter(site == 'SRS6', group == 'Consumer') %>% rename(d13C = md13C, d15N = md15N, d34S = md34S)

write.csv(SRS6mix, "SRS6mix.csv", row.names = FALSE)

mix <- load_mix_data(filename = "SRS6mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename = "data/sourcesSRS6_2.csv",
                           source_factors = NULL,
                           conc_dep = T,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS6_1.csv"), mix)

#plot_data(filename = "SRS6_isospace_plot", plot_save_pdf = TRUE, plot_save_png = TRUE, mix, source, discr)

model_filename <- "SRS6_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.SRS6 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.SRS6 <- run_model(run="very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)

output_jags.SRS6  <- list(summary_save = TRUE,
                          summary_name = "FCESRS6_sumstats",
                          sup_post = FALSE,
                          plot_post_save_pdf = T,
                          plot_post_name = "FCESRS6_plot",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "FCESRS6_pairs",
                          sup_xy = T,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "FCESRS6_plot",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = T,
                          diag_save = F,
                          diag_name = "FCESRS6_Diagnostic",
                          indiv_effect = FALSE,
                          plot_post_save_png = T,
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = F)

output_JAGS(jags.SRS6 , mix, source, output_jags.SRS6)

##combine posterior ground into brown/green

#p.dry19.Macroalgae    0.036 0.035 0.001 0.002 0.010 0.026 0.050 0.106 0.125
#p.wet19.Macroalgae    0.014 0.016 0.000 0.001 0.003 0.009 0.019 0.046 0.057
#p.dry19.Mangrove      0.467 0.051 0.371 0.388 0.433 0.466 0.499 0.554 0.573
#p.wet19.Mangrove      0.848 0.058 0.733 0.752 0.808 0.849 0.887 0.946 0.957
#p.dry19.Phytoplankton 0.497 0.050 0.390 0.409 0.465 0.499 0.531 0.574 0.588
#p.wet19.Phytoplankton 0.138 0.054 0.034 0.045 0.102 0.138 0.174 0.227 0.244


combinedSRS6 <- combine_sources(jags.SRS6, mix, source, alpha.prior=1, 
                                groups=list(green=c('Phytoplankton','Macroalgae' ), brown=c('Mangrove')))




# get posterior medians for new source groupings
apply(combinedSRS6$post, 2, median)
summary_stat(combinedSRS6, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "SRS6_combined_sumstats" )


#               Mean    SD  2.5%   25%   50%   75% 97.5%
#p.green.dry19 0.533 0.051 0.427 0.501 0.534 0.567 0.629
#p.brown.dry19 0.467 0.051 0.371 0.433 0.466 0.499 0.573
#p.green.wet19 0.152 0.058 0.043 0.113 0.151 0.192 0.267
#p.brown.wet19 0.848 0.058 0.733 0.808 0.849 0.887 0.957


##################################SRS4


SRS4mix<-SIa %>% filter(site == 'SRS4', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


write.csv(SRS4mix,"SRS4mix.csv",row.names = F) 

mix <- load_mix_data(filename="SRS4mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/sourcesSRS4_2.csv",
                           source_factors = NULL,
                           conc_dep = T,
                           data_type = "means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS4_2.csv"), mix)

plot_data(filename="SRS4_isospace_plot", plot_save_pdf=T, plot_save_png=T, mix,source,discr)

model_filename <- "SRS4_mix_test.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.SRS4 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.SRS4 <- run_model(run="very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)


output_jags.SRS4  <- list(summary_save = TRUE,
                          summary_name = "FCESRS4_sumstats_demo",
                          sup_post = FALSE,
                          plot_post_save_pdf = T,
                          plot_post_name = "FCESRS4_plot_demo",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "FCESRS4_pairs_demo",
                          sup_xy = T,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "FCESRS4_plot_demo",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = T,
                          diag_save = F,
                          diag_name = "FCESRS4_Diagnostic_demo",
                          indiv_effect = FALSE,
                          plot_post_save_png = T,
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = F)

output_JAGS(jags.SRS4 , mix, source, output_jags.SRS4)

##combine posterior ground into brown/green

#                             
#                             Mean    SD  2.5%    5%   25%   50%   75%   95% 97.5%
#p.dry19.Epiphytic microalgae 0.536 0.150 0.256 0.288 0.431 0.534 0.639 0.786 0.844
#p.wet19.Epiphytic microalgae 0.474 0.181 0.186 0.216 0.334 0.455 0.592 0.817 0.884
#p.dry19.Mangrove             0.281 0.094 0.076 0.108 0.220 0.287 0.350 0.420 0.443
#p.wet19.Mangrove             0.299 0.098 0.045 0.096 0.247 0.316 0.371 0.423 0.442
#p.dry19.Phytoplankton        0.183 0.092 0.026 0.038 0.113 0.179 0.245 0.341 0.370
#p.wet19.Phytoplankton        0.228 0.112 0.014 0.031 0.146 0.236 0.313 0.399 0.420

combinedSRS4 <- combine_sources(jags.SRS4, mix, source, alpha.prior=1, 
                                groups=list(green=c('Phytoplankton',"Epiphytic microalgae"), brown=c('Mangrove')))

# get posterior medians for new source groupings
apply(combinedSRS4$post, 2, median)
summary_stat(combinedSRS4, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "SRS4_combined_sumstats_demo" )


#              Mean    SD  2.5%   25%   50%   75% 97.5%
#p.green.dry19 0.722 0.097 0.558 0.652 0.711 0.783 0.939
#p.brown.dry19 0.278 0.097 0.061 0.217 0.289 0.348 0.442
#p.green.wet19 0.706 0.100 0.560 0.634 0.685 0.763 0.955
#p.brown.wet19 0.294 0.100 0.045 0.237 0.315 0.366 0.440
#


#### SRS3 ----
SRS3mix<-SI %>% filter(site == 'SRS3', common_name!="Egyptian paspalidium",group=='Consumer')
# %>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


write.csv(SRS3mix,"SRS3mix.csv",row.names = F) 

mix <- load_mix_data(filename="data/SRS3mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/sourcesSRS3_1.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_SRS3_1.csv"), mix)

plot_data(filename="SRS3_isospace_plot", plot_save_pdf=FALSE, plot_save_png=F, mix,source,discr)

model_filename <- "SRS3_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.SRS3 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.SRS3 <- run_model(run="very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)

output_jags.SRS3  <- list(summary_save = TRUE,
                          summary_name = "FCESRS3_sumstats",
                          sup_post = FALSE,
                          plot_post_save_pdf = F,
                          plot_post_name = "FCESRS3_plot",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "FCESRS3_pairs",
                          sup_xy = F,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "FCESRS3_plot",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = F,
                          diag_save = F,
                          diag_name = "FCESRS3_Diagnostic",
                          indiv_effect = FALSE,
                          plot_post_save_png = F,
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = F)

output_JAGS(jags.SRS3 , mix, source, output_jags.SRS3)



#############RB10


RB10mix<-SIa %>% filter(site == 'RB10', common_name!="Egyptian paspalidium",group=='Consumer')%>% rename('d13C'='md13C',"d15N"= "md15N","d34S"="md34S")


write.csv(RB10mix,"RB10mix.csv",row.names = F) 

mix <- load_mix_data(filename="data/RB10mix.csv",
                     iso_names=c("d13C","d15N","d34S"),
                     factors=c('common_name','season'),
                     fac_random=c(F,F),
                     fac_nested=c(F,F),
                     cont_effects=NULL)

source <- load_source_data(filename="data/sources_RB10_1.csv",
                           source_factors=NULL,
                           conc_dep=T,
                           data_type="means",
                           mix)

discr <- load_discr_data(file("data/FCE_TEF_RB10_1.csv"), mix)

plot_data(filename="RB10_isospace_plot", plot_save_pdf=FALSE, plot_save_png=F, mix,source,discr)

model_filename <- "RB10_mix.txt"
write_JAGS_model(model_filename, resid_err=F, process_err=T, mix, source)


#run a test model to make sure it works
jags.RB10 <- run_model(run="test", mix, source, discr, model_filename, 
                       alpha.prior = 1, resid_err=F, process_err=T)


jags.RB10 <- run_model(run="very long", mix, source, discr, model_filename,
                       alpha.prior = 1, resid_err=F, process_err=T)

output_jags.RB10  <- list(summary_save = TRUE,
                          summary_name = "FCERB10_sumstats",
                          sup_post = FALSE,
                          plot_post_save_pdf = T,
                          plot_post_name = "FCERB10_plot",
                          sup_pairs = FALSE,
                          plot_pairs_save_pdf = F,
                          plot_pairs_name = "FCERB10_pairs",
                          sup_xy = T,
                          plot_xy_save_pdf = F,
                          plot_xy_name = "FCERB10_plot",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = T,
                          diag_save = F,
                          diag_name = "FCERB10_Diagnostic",
                          indiv_effect = FALSE,
                          plot_post_save_png = T,
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = F)

output_JAGS(jags.RB10 , mix, source, output_jags.RB10)

#                         Mean    SD  2.5%    5%   25%   50%   75%   95% 97.5%
#p.dry19.Epiphytic microalgae 0.845 0.066 0.703 0.730 0.801 0.851 0.894 0.946 0.957
#p.wet19.Epiphytic microalgae 0.694 0.071 0.564 0.582 0.647 0.690 0.739 0.816 0.842
#p.dry19.Mangrove             0.018 0.017 0.001 0.001 0.006 0.013 0.026 0.053 0.064
#p.wet19.Mangrove             0.025 0.025 0.000 0.001 0.006 0.016 0.036 0.077 0.092
#p.dry19.Phytoplankton        0.137 0.065 0.032 0.040 0.089 0.131 0.179 0.254 0.277
#p.wet19.Phytoplankton        0.281 0.071 0.134 0.162 0.236 0.284 0.331 0.390 0.412






combinedRB10 <- combine_sources(jags.RB10, mix, source, alpha.prior=1, 
                                groups=list(green=c('Phytoplankton','Epiphytes' ), brown=c('Mangrove')))




# get posterior medians for new source groupings
apply(combinedRB10$post, 2, median)
summary_stat(combinedRB10, meanSD=T, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), savetxt=T,
             filename = "RB10_combined_sumstats" )


#             Mean    SD  2.5%   25%   50%   75% 97.5%
#p.green.dry19 0.982 0.017 0.936 0.974 0.987 0.994 0.999
#p.brown.dry19 0.018 0.017 0.001 0.006 0.013 0.026 0.064
#p.green.wet19 0.975 0.025 0.908 0.964 0.984 0.994 1.000
#p.brown.wet19 0.025 0.025 0.000 0.006 0.016 0.036 0.092







#### Food Web Level Box Plots ----
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


