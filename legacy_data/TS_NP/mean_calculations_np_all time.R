getwd()
setwd("/Users/mack/Desktop/RESEARCH/FCE FOOD WEBS/FOOD WEBS R/")

################ FOR TAYLOR SLOUGH/FLORIDA BAY SITES ONLY ###################

datALL = read.csv("TS_NP_ALL.csv")
preWILMA = read.csv("TS_NP_PREWILMA.csv")
postWILMA = read.csv("TS_NP_POSTWILMA.csv")
postIRMA = read.csv("postirma_goods.csv")

##### CALCULATE MEAN FOR EACH POR

site_means = datALL %>%
  group_by(SITENAME) %>%
  summarise_at(vars(NP), list(name = mean))
View(site_means)

# TS 3 = 226
# TS 7 = 162
# TS 9 = 103
# TS 10 = 100
# TS 11 = 31

# PRE WILMA TS DATA

PREW_means = preWILMA %>%
  group_by(SITENAME) %>%
  summarise_at(vars(NP), list(name = mean))
View(PREW_means)

# TS 3 = 310
# TS 7 = 207
# TS 9 = 107
# TS 10 = 108
# TS 11 = 31

# POST WILMA TS DATA

POSTW_means = postWILMA %>%
  group_by(SITENAME) %>%
  summarise_at(vars(NP), list(name = mean))
View(POSTW_means)

# TS 3 = 218
# TS 7 = 125
# TS 9 = 103
# TS 10 = 96
# TS 11 = 33

# POST IRMA TS DATA

head(postIRMA)

POSTI_means = postIRMA %>%
  group_by(SITENAME) %>%
  summarise_at(vars(NP), list(name = mean))
View(POSTI_means)

# TS 3 = 223
# TS 7 = 79
# TS 9 = 81
# TS 10 = 78
# TS 11 = 30


###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


########### CLEAR OUT THE GLOBAL ENVIRONMENT AND START OVER WITH #######
################ SHARK RIVER SLOUGH DATA - CODE BELOW ##################

setwd("/Users/mack/Desktop/R")
datALL <- read.csv("SRS_NP_ALL.csv")
View(datALL)
# OCTOBER 2000 - DECEMBER 2020

preWILMA <- read.csv("SRS_NP_PREWILMA.csv")
View(preWILMA)
# OCTOBER 2000 - SEPT 2005

postWILMA <- read.csv("SRS_NP_POSTWILMA.csv")
View(postWILMA)
# OCTOBER 2005 - DEC 2020

postIRMA <- read.csv("SRS_NP_POSTIRMA.csv")
View(postIRMA)
# AUGUST 2017 - DEC 2020

# INITIALLY FOCUS ON POST WILMA (2005 THROUGH 2020)
# AS SUCH, SUBSET BY IND SITE AND FIGURE OUT MEAN AND OTHER IMPORTANT INFO

SRS1 = postWILMA %>% filter(postWILMA$Site == "SRS1")

SRS2 = postWILMA %>% filter(postWILMA$Site == "SRS2")

SRS3 = postWILMA %>% filter(postWILMA$Site == "SRS3")

SRS4 = postWILMA %>% filter(postWILMA$Site == "SRS4")

SRS5 = postWILMA %>% filter(postWILMA$Site == "SRS5")

SRS6 = postWILMA %>% filter(postWILMA$Site == "SRS6")


# FIND MEAN FOR EACH SITE ACROSS ENTIRE DATA SET FOR MASS

site_means = datALL %>%
  group_by(Site) %>%
  summarise_at(vars(NP_mass), list(name = mean, name = sd))
View(site_means)

#	SRS1 = 813.9556
# SRS2 = 928.6240
# SRS3 = 1363.5060
# SRS4 = 468.1940
# SRS5 = 420.1753
# SRS6 = 313.2003

# FIND MEAN FOR EACH SITE ACROSS POST WILMA DATA SET FOR MASS

site_means_WILMA = postWILMA %>%
  group_by(Site) %>%
  summarise_at(vars(NP_mass), list(name = mean, name = sd))
View(site_means_WILMA)

# SRS1 = 563.3792
# SRS2 = 997.1553
# SRS3 = 1094.9784
# SRS4 = 271.2009
# SRS5 = 207.9369
# SRS6 = 279.7601

# FIND MEAN FOR EACH SITE ACROSS POST WILMA DATA SET FOR MASS

site_means_PREWILMA = preWILMA %>%
  group_by(Site) %>%
  summarise_at(vars(NP_mass), list(name = mean)) 
View(site_means_PREWILMA)

# SRS1 = 1253.5142
# SRS2 = 800.2787
# SRS3 = 1806.5314
# SRS4 = 1113.2771
# SRS5 = 1155.0661
# SRS6 = 411.0686

# FIND MEAN FOR EACH SITE ACROSS POST IRMA DATA SET FOR MASS

site_means_POSTIRMA = postIRMA %>%
  group_by(Site) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(site_means_POSTIRMA)

# SRS1 = 81.51961
# SRS2 = 70.63400
# SRS3 = 3562.41976
# SRS4 = 34.27556
# SRS5 = 101.02142
# SRS6 = 14.53505

# FOR POST IRMA WITHOUT 

site_means_POSTIRMA2 = postIRMA2 %>%
  group_by(Site) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(site_means_POSTIRMA2)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

postIRMA %>%
  group_by(Site) %>%
  mutate(outlier = ifelse(is_outlier(NP_mass), NP_mass, as.numeric(NA))) %>%
  ggplot(., aes(x = factor(Site), y = NP_mass)) +
  geom_boxplot() +
  geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3)
          
# filter data out that may be considered an outlier

outliers = boxplot(postIRMA$NP_mass, plot = FALSE)$out

postIRMA_no = postIRMA

postIRMA_no = postIRMA_no[-which(postIRMA_no$NP_mass %in% outliers)]

# FIND MEAN FOR EACH SITE, MONTHLY ACROSS POST WILMA DATA SET

View(SRS1)

SRS1_month_means = SRS1 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS1_month_means)

SRS2_month_means = SRS2 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS2_month_means)

SRS3_month_means = SRS3 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS3_month_means)

SRS4_month_means = SRS4 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS4_month_means)

SRS5_month_means = SRS5 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS5_month_means)

SRS6_month_means = SRS6 %>%
  group_by(Month) %>%
  summarise_at(vars(NP_mass), list(name = mean))
View(SRS6_month_means)
