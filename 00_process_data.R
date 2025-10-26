library(tidyverse)
library(dplyr)
library(read.dbc)
library(ggplot2)
library(lubridate)
  

# Dengue ##################################################################

file_list <- list.files('data/dengue/', pattern = "\\.dbc$", full.names = TRUE)
file_list <- file_list[4:length(file_list)]
#file_list <- file_list[-1]
df <- data.frame()

for(file in file_list){
  print(file)
  df_aux <- read.dbc(file)
  df_aux <- df_aux %>% filter(CLASSI_FIN %in% c(1, 2, 3, 4, 10, 11, 12))
  df_aux <- df_aux %>% mutate(anoepi = epiyear(DT_SIN_PRI))
  df_aux <- df_aux %>% mutate(epiweek = epiweek(DT_SIN_PRI))
  df_aux <- df_aux %>% group_by(anoepi, epiweek, ID_MN_RESI) %>% summarise(n = n())
  df_aux <- df_aux %>% ungroup()
  df <- rbind(df, df_aux)
}

df <- df %>% filter(anoepi >= 2010) %>% filter(anoepi <= 2023)

df <- df %>% group_by(anoepi, epiweek, ID_MN_RESI) %>% summarise(n = sum(n))
save(df, file = 'timeseries_dengue_2010_2023.RData')

