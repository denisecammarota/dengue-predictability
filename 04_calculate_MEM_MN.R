# Calculating activity thresholds for MUNS ###################
# Using the moving epidemic method (MEM) ####################

library(tidyverse)
library(sf)
library(lubridate)
library(tidyverse)
library(dplyr)
library(readr)
library(zoo)
library(geobr)
library(ggthemes)
library(patchwork)
library(viridis)
library(ggplot2)
library(reshape2)
library(mem)

# Loading the names of municipalities ###########################
map_br <- read_municipality()
map_br <- map_br %>% mutate(code_muni_6 = substr(code_muni, 1, 6))
unique_muns <- unique(map_br$code_muni_6)

# Loading incidence data ########################################
df_res <- read.csv('timeseries_dengue_lambda_2010_2023.csv')

# Shifting seasons for doing MEM ########################
df_res <- df_res %>% 
  mutate(season = ifelse(epiweek <= 40, 
                         paste0(anoepi - 1, '/', anoepi), 
                         paste0(anoepi, '/', anoepi + 1)))
df_res <- df_res %>% filter(!season %in% c('2009/2010', '2023/2024'))
df_res <- df_res %>% filter(epiweek != 53)
df_res

df_res <- df_res %>% mutate(epiweek_aux = ifelse(epiweek >= 40, epiweek - 40 + 1, epiweek + 52 - 40 + 1))

# Selecting municipalities that are relevant ###########

df_res <- df_res %>% filter(ID_MN_RESI %in% unique_muns)
a <- df_res %>% mutate(bin = ifelse(inc == 0, 1, 0)) %>%
  group_by(ID_MN_RESI) %>% summarise(n_bin = sum(bin))
a <- a %>% filter(n_bin < 500)
sel_muns <- unique(a$ID_MN_RESI)

df_res <- df_res %>% filter(ID_MN_RESI %in% sel_muns)
df_res <- df_res %>% select(season, epiweek_aux, ID_MN_RESI, inc)

# Eliminating seasons with extreme incidence ########################

# Calculating peaks and their geom mean and distance
df_sum <- df_res %>% group_by(ID_MN_RESI, season) %>% summarise(inc_max = max(inc))
df_mean_dist <- df_sum %>% group_by(ID_MN_RESI) %>% summarise(geom_mean = exp(mean(log(inc_max + 1))) - 1)
df_sum <- df_sum %>% left_join(df_mean_dist, by = join_by(ID_MN_RESI))
rm(df_mean_dist)
df_sum <- df_sum %>% mutate(geom_dist = inc_max/geom_mean)
df_sum

# Calculating the distribution of geom mean peaks
perc_10_90 <- quantile(df_sum$geom_dist, probs = c(0.10, 0.90))

# Filtering by the percentile criterion
df_sum <- df_sum %>% mutate(keep = case_when(
  geom_dist >= perc_10_90[1] & geom_dist <= perc_10_90[2] ~ 1,
  T ~ 0
))

# Checking the number that remain for each state - ones with < 5 seasons
df_tmp_seasons <- df_sum %>% group_by(ID_MN_RESI) %>% 
  summarise(seasons_keep = sum(keep)) %>%
  filter(seasons_keep < 5) %>%
  mutate(seasons_need = 5 - seasons_keep)
df_tmp_seasons$ID_MN_RESI
df_tmp_seasons$seasons_need

# Adding seasons that are needed 
for(id_mun in unique(df_tmp_seasons$ID_MN_RESI)){
  n_seasons_need <- df_tmp_seasons %>% filter(ID_MN_RESI == id_mun)
  n_seasons_need <- n_seasons_need$seasons_need[1]
  df_tmp_sg <- df_sum %>% filter(ID_MN_RESI == id_mun) %>% filter(keep == 0)
  df_tmp_sg <- df_tmp_sg %>% arrange(abs(geom_dist - perc_10_90[1]), abs(geom_dist - perc_10_90[2]))
  df_tmp_sg <- df_tmp_sg %>% slice_head(n = n_seasons_need) %>% mutate(keep = 1)
  df_sum <- rbind(df_sum, df_tmp_sg)
}

rm(n_seasons_need, df_tmp_sg, df_tmp_seasons)
gc()

df_sum <- df_sum %>% group_by(ID_MN_RESI, season) %>% summarise(keep = max(keep))
df_mem <- df_res %>% left_join(df_sum, by = join_by(ID_MN_RESI, season)) 
df_mem <- df_mem %>% filter(keep == 1)
df_mem <- df_mem %>% select(season, epiweek_aux, 
                            ID_MN_RESI, inc)

list_pre <- list()
list_pos <- list()
list_start <- list()
list_duration <- list()
i = 1
for(id_mun in unique(df_mem$ID_MN_RESI)){
  print(id_mun)
  df_sp <- df_mem %>% filter(ID_MN_RESI == id_mun) %>% 
    select(season, epiweek_aux, inc) %>%
    pivot_wider(names_from = season, values_from = inc)
  list_epiweek <- df_sp$epiweek_aux
  df_sp <- df_sp %>% select(!epiweek_aux)
  dengue.memmodel <- memmodel(df_sp, i.season = 30, i.method = 2)
  plot(dengue.memmodel)
  pre_aux <- as.numeric(dengue.memmodel$pre.post.intervals["pre.i", 3])
  pos_aux <- as.numeric(dengue.memmodel$pre.post.intervals["post.i", 3])
  start_aux <- as.numeric(dengue.memmodel[["mean.start"]])
  dur_aux <- as.numeric(dengue.memmodel[["mean.length"]])
  list_pre[[i]] <- pre_aux[1]
  list_pos[[i]] <- pos_aux[1]
  list_start[[i]] <- start_aux[1]
  list_duration[[i]] <- dur_aux[1]
  i = i + 1
}

list_pre <- unlist(list_pre)
list_pos <- unlist(list_pos)
list_start <- unlist(list_start)
list_duration <- unlist(list_duration)

df_summary <- data.frame(ID_MN_RESI = unique(df_mem$ID_MN_RESI), 
                         pre_thr = list_pre,
                         pos_thr = list_pos,
                         start_time = list_start,
                         dur_time = list_duration)
df_summary <- df_summary %>% mutate(start_time_real = ifelse(start_time <= 13, start_time + 39, start_time - 13))

# Making maps of quantities ###########################

map_br <- map_br %>% mutate(code_muni_6 = as.numeric(code_muni_6))
map_br <- map_br %>% left_join(df_summary, by = join_by(code_muni_6 == ID_MN_RESI))

ggplot(map_br, aes(fill = pre_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Pre-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$pre_thr, na.rm = TRUE),  # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = pos_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Post-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$pos_thr, na.rm = TRUE),         # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = start_time_real)) + geom_sf() + theme_map() +
  labs(fill = 'Start Week') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$start_time_real, na.rm = TRUE),   # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = dur_time)) + geom_sf() + theme_map() +
  labs(fill = 'Duration') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$dur_time, na.rm = TRUE), # set your desired middle point
    guide = "colorbar"
  )





