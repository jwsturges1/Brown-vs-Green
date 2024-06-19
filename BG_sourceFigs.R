#' """ Make basal resource figures for each transect with shaded background
#'     @author: Ryan James, James Sturges
#'     date: 9/22/23"""

library(tidyverse)

# load data ----
MixOut_RB10 = read.csv('data/Mix_Quants/MT_RB10.csv')
MixOut_SRS3 = read.csv('data/Mix_Quants/MT_SRS3.csv')
MixOut_SRS4 = read.csv('data/Mix_Quants/MT_SRS4.csv')
MixOut_SRS6 = read.csv('data/Mix_Quants/MT_SRS6.csv')
MixOut_TS3 = read.csv('data/Mix_Quants/MT_TS3.csv')
MixOut_TS7 = read.csv('data/Mix_Quants/MT_TS7.csv')
MixOut_TS9 = read.csv('data/Mix_Quants/MT_TS9.csv')
MixOut_TS10 = read.csv('data/Mix_Quants/MT_TS10.csv')
MixOut_TS11 = read.csv('data/Mix_Quants/MT_TS11.csv')

# SRSMixout_gb <- bind_rows(MixOut_RB10, MixOut_SRS3, MixOut_SRS4, MixOut_SRS6)
# 
# SRSMixout_gb = SRSMixout_gb %>% 
#   rename(site = type, season = code) %>%
#   mutate(path = case_when(
#     source %in% c("Epiphytes", "Phytoplankton", "Filamentous Green Algae", "Periphyton") ~ "green",
#     source %in% c("Mangrove", "Sawgrass", "Red Macroalgae",'Floc') ~ "brown",
#     TRUE ~ NA_character_  # For other cases, you can assign NA or something else if needed
#   ))
# 
# SRSMixout_gb = SRSMixout_gb %>% 
#   group_by(name, site, season, path) %>% 
#   summarize(value = sum(mean)) %>% 
#   pivot_wider(names_from = path, values_from = value)
# 
# TSMixout_gb <- rbind(MixOut_TS3, MixOut_TS7, MixOut_TS9, MixOut_TS10, MixOut_TS11)
# 
# TSMixout_gb = TSMixout_gb %>% 
#   rename(site = type, season = code) %>%
#   mutate(path = case_when(
#     source %in% c("Epiphytes", "Phytoplankton", "Filamentous Green Algae", "Periphyton", "Epiphytic microalgae", 'SPOM') ~ "green",
#     source %in% c("Mangrove", "Sawgrass", "Red Macroalgae", "Seagrass", 'Floc') ~ "brown",
#     TRUE ~ NA_character_  # For other cases, you can assign NA or something else if needed
#   ))
# 
# TSMixout_gb = TSMixout_gb %>% 
#   group_by(name, site, season, path) %>% 
#   summarize(value = sum(mean)) %>% 
#   pivot_wider(names_from = path, values_from = value)
# 
# 
# y_label_formatter <- function(x) {
#   ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
# }
# 
# 
# combined_df = bind_rows(SRSMixout_gb,TSMixout_gb) %>%  
#   mutate(transect = case_when(
#     site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
#     site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"))
# 
# unique_names <- combined_df %>%
#   group_by(site, season) %>%
#   summarize(unique_names = n_distinct(name))
# 
# 
# # Combined Brown vs Green boxplot----
# combined_df <- combined_df %>%
#   mutate(site = ifelse(site == "RB10", "Upper River",
#                        ifelse(site == "SRS3", "SRS Marsh",
#                               ifelse(site == "SRS4", "Mid River",
#                                      ifelse(site == "SRS6", "Lower River",
#                                             ifelse(site == "TS3", "TS Marsh",
#                                                    ifelse(site == "TS7", "Mangrove Ecotone",
#                                                           ifelse(site == "TS9", "Inner Bay",
#                                                                  ifelse(site == "TS10", "Mid Bay",
#                                                                         ifelse(site == "TS11", "Outer Bay", site)
#                                                                  )
#                                                           )
#                                                    )
#                                             )
#                                      )
#                               )
#                        )
#   ),
#   site = factor(site, levels = c( "SRS Marsh","Upper River","Mid River","Lower River", "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")))
# 
# combined_df = combined_df %>% 
#   group_by(site, season) %>% 
#   mutate(fill = mean(green),
#          site = factor(site, levels = c( "Outer Bay","Mid Bay","Inner Bay","Mangrove Ecotone" ,"TS Marsh", "Lower River", "Mid River","Upper River" ,"SRS Marsh" )))
# 
# 
# combined dataset from mixing model outputs
s_df= tibble(source = c("Sawgrass", "Mangrove",
                         "Floc","Red Macroalgae","Seagrass",
                         "Epiphytes","Periphyton","Filamentous Green Algae",
                         'SPOM',"Phytoplankton"),
             lab = c('Sawgrass', 'Man',
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
  mutate(lab = factor(lab, levels = c('Sawgrass', 'Man',
                                      "Floc","RMA","Seagrass",
                                      'EMA','Peri','FGA',
                                      'POM', 'Phyto')),
         lb = factor(lab, levels = c( "Man","Sawgrass","Floc","RMA",
                                      "Seagrass", "EMA","Peri" ,
                                      "Phyto",  "FGA", "POM")))

names = tibble(site = c("RB10","SRS3", "SRS4","SRS6",
                        "TS3", "TS7", "TS9", "TS10", "TS11"),
               slab = c("Upper River", "SRS Marsh", "Mid River", "Lower River",
                         "TS Marsh","Mangrove Ecotone","Inner Bay", "Mid Bay", "Outer Bay"))

cont_df = bind_rows(MixOut_RB10, MixOut_SRS3, MixOut_SRS4, MixOut_SRS6, MixOut_TS3, MixOut_TS7, MixOut_TS9, MixOut_TS10, MixOut_TS11) |> 
  as_tibble() |> 
  rename(site = type, season = code) |> 
  left_join(s_df, by = 'source') |> 
  left_join(names, by = 'site') |> 
  mutate(transect = case_when(
    site %in% c("SRS3", "SRS4","SRS6", "RB10") ~ "Shark River Slough",
    site %in% c("TS3", "TS7", "TS9", "TS10", "TS11") ~ "Taylor Slough"),
    slab = factor(slab, levels = c( "SRS Marsh","Upper River","Mid River","Lower River", "TS Marsh", "Mangrove Ecotone", "Inner Bay", "Mid Bay", "Outer Bay")))

#663300 for brown and #92D050 for green
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
mutate(lab = factor(lab, levels = c('Sawgrass', 'Man',
                                    "Floc","RMA","Seagrass",
                                    'EMA','Peri','FGA',
                                    'PMA')))

ggplot(data = cont_df, aes(x = lab, y = mean, fill = season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = cont_df |> group_by(slab) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = cont_df |> group_by(slab) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~slab, scales = "free_x") +
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
  # theme(
  #   axis.title = element_text(size = 32), 
  #   axis.text.y = element_text(size = 20, colour = "black", face = "bold"), 
  #   axis.text.x = element_text(
  #     size = 18,
  #     colour = 'black'), 
  #   plot.title = element_text(size = 18, hjust = 0.5),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   legend.position = 'top',
  #   legend.title = element_text(size = 32, face = "bold"),
  #   strip.text.x = element_text(size = 26, face = "bold"),
  #   strip.text.y = element_text(size = 18, face = "bold"),
  #   legend.text = element_text(size = 24, face = "bold"))

ggsave("figures/source_cont_plot_GBcol.png", width = 9.5, height = 7.5, dpi = 600)

# SRS only 
srs = cont_df |> filter(transect == "Shark River Slough")
# SRS only 
ggplot(data = srs, aes(x = lab, y = mean, fill = season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = srs |> group_by(slab) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = srs |> group_by(slab) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~slab, scales = "free_x", drop = T, ncol = 1) +
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
  # theme(
  #   axis.title = element_text(size = 14), 
  #   axis.text.y = element_text(size = 14, colour = "black", face = "bold"), 
  #   axis.text.x = element_text(
  #     size = 12,
  #     colour = 'black'), 
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   legend.position = 'top',
  #   legend.title = element_text(size = 14, face = "bold"),
  #   strip.text = element_text(size = 14, face = "bold"),
  #   legend.text = element_text(size = 12, face = "bold"))

ggsave("figures/Source_SRS_tall.png", width = 4, height = 8, dpi = 600)

# SRS only 
ggplot(data = srs, aes(x = lab, y = mean, fill = season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = srs |> group_by(slab) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = srs |> group_by(slab) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~slab, scales = "free_x", drop = T, nrow = 1) +
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

ggsave("figures/Source_SRS_wide.png", width = 12, height = 4, dpi = 600)

# TS only 
ts = cont_df |> filter(transect == "Taylor Slough")
# TS only 
ggplot(data = ts, aes(x = lab, y = mean, fill = season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = ts |> group_by(slab) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = ts |> group_by(slab) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~slab, scales = "free_x", drop = T, ncol = 1) +
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

ggsave("figures/Source_TS_tall.png", width = 4, height = 10, dpi = 600)

# TS only 
ggplot(data = ts, aes(x = lab, y = mean, fill = season)) +
  geom_boxplot() + # gets drawn over but sets order of axes not sure why it changes order when start with geom_rect
  geom_rect(data = ts |> group_by(slab) |> slice(1),
            aes(xmax = Inf, xmin = 2.5, ymax = Inf, ymin =-Inf), fill = '#92D050', alpha = 0.4)+
  geom_rect(data = ts |> group_by(slab) |> slice(1),
            aes(xmin = -Inf, xmax = 2.5, ymax = Inf, ymin =-Inf), fill = '#663300', alpha = 0.4)+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 'dashed', alpha = 0.4)  +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~slab, scales = "free_x", drop = T, nrow = 1) +
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

ggsave("figures/Source_TS_wide.png", width = 15, height = 4, dpi = 600)

