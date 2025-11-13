# Put together a dataset for ML, using epidemiological and climate data #####################
# using the previously selected municipalities through the MEM method ######################

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

# Loading MEM municipalities data ######################################
df_muns <- read.csv('mem_municipalities.csv')
list_muns <- unique(df_muns$ID_MN_RESI)

# Loading epidemiological data ########################################
df_res <- read.csv('timeseries_dengue_lambda_2010_2023.csv')
df_res <- df_res %>% 
  mutate(season = ifelse(epiweek <= 40, 
                         paste0(anoepi - 1, '/', anoepi), 
                         paste0(anoepi, '/', anoepi + 1)))
df_res <- df_res %>% filter(!season %in% c('2009/2010', '2023/2024'))
df_res <- df_res %>% filter(epiweek != 53)
df_res <- df_res %>% mutate(epiweek_aux = ifelse(epiweek >= 40, epiweek - 40 + 1, epiweek + 52 - 40 + 1))
df_res <- df_res %>% filter(ID_MN_RESI %in% list_muns)
df_muns <- df_muns %>% select(ID_MN_RESI, pre_thr)
df_res <- df_res %>% left_join(df_muns, by = join_by(ID_MN_RESI))
colnames(df_res)
df_res <- df_res %>% mutate(inc_bin = ifelse(inc >= pre_thr, 1, 0))
df_res <- df_res %>% mutate(inc_lambda_bin = ifelse(inc_lambda >= 1, 1, 0),
                            smooth_lambda_bin = ifelse(smooth_lambda >= 1, 1, 0))
a <- df_res %>% filter(inc_bin == 1) %>% 
  group_by(ID_MN_RESI, season) %>% 
  summarise(week_aux_first = min(epiweek_aux),
            week_aux_last = max(epiweek_aux))

b <- df_res %>% filter(ID_MN_RESI == 110011) %>% filter(season == '2021/2022')
#b <- b %>% mutate(tplot_aux = epiweek_aux + anoepi)

# Has a few problems 
# - some years start and then end at high volumes
# - some years have no epidemics really, how do we compute lambda then?

# Loading climate data ###############################################
df_clim <- read.csv('data/climate_mosqlimate/climate.csv')
colnames(df_clim)
df_clim <- df_clim %>% filter(id_mun %in% list_muns)

# Joining all data together #########################################
df_total


# Grouping and summarizing ##########################################
df_aux <- df_total %>% filter(inter_epi = ifelse(epiweek_ >= 40))
df_total <- df_total %>% group_by()









