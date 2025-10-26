# Calculating the index proposed by Rypdal and Sugihara, 2019 ##############

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

# Carregando dados de series temporais #############################
load('timeseries_dengue_2010_2023.RData')
df <- df %>% ungroup()
df <- df %>% complete(anoepi, epiweek, ID_MN_RESI)
valid_53_week_years <- c(2015, 2020)
df <- df %>% filter(epiweek != 53 | anoepi %in% valid_53_week_years)
df <- df %>% filter(anoepi >= 2010) %>% filter(anoepi <= 2023)
df <- df %>% mutate(n = replace_na(n, 0))
df <- df %>% mutate(tplot = anoepi + (epiweek/52))
df <- df %>% select(anoepi, epiweek, tplot, ID_MN_RESI, n)
df_cases <- df
rm(df)

# Carregando dados populacionais ###################################
load('pop_muns.RData')
df_pop <- df
df_pop <- df_pop %>% select('Ano', 'Mun', 'Total')
df_pop <- df_pop %>% mutate(Ano = as.numeric(Ano))
#df_pop <- df_pop %>% filter(Ano == 2016)
#df_pop <- df_pop %>% select(!Ano)
rm(df)
gc()

# Calculando incidencias semanais ##################################
df_cases <- df_cases %>% left_join(df_pop, by = join_by(ID_MN_RESI == Mun,
                                                        anoepi == Ano))
df_cases <- df_cases %>% mutate(inc = (10**5)*(n/Total))
df_cases <- df_cases %>% select(anoepi, epiweek, tplot, ID_MN_RESI, inc)
df_cases <- df_cases %>% drop_na()

# Calculando o susceptibility proxy (Rypdal and Sugihara, 2019) ###########

## Defining lambda function #####################################
lambda <- function(x){
  y <- lead(x, 1)
  x <- x[1:(length(x)-1)]
  y <- y[!(is.na(y))]
  coeflm <- as.numeric(lm(y ~ x - 1)$coefficients[1])
  return(coeflm)
}

## Calculando para cada serie temporal ###########################
# Pular se ja foi executada essa parte ###########################

ws1 <- 12 # window size
df_res <- data.frame()
muns_sel <- unique(df_cases$ID_MN_RESI)
#muns_sel <- c(310620)

i = 0

for(mun in muns_sel){
  print(i)
  df_aux <- df_cases %>% filter(ID_MN_RESI == mun)
  df_aux['inc_lambda'] <- rollapply(df_aux['inc'], width = ws1, FUN = function(w) lambda(w), align = "right", by.column = FALSE, fill = NA)
  df_aux['smooth_lambda'] <- rollmean(df_aux['inc_lambda'], k = 10, align = "right", fill = NA)
  df_res <- rbind(df_res, df_aux)
  i = i + 1
}

write_csv(df_res, 'timeseries_dengue_lambda_2010_2023.csv')

# Carregando o resultado #####################################
#df_res <- read.csv('timeseries_dengue_lambda_2010_2023.csv')

# Ate agora:
# df_pop dados pop
# df_cases incidencias semanais
# df_res valores de lambda

muns <- c(355030, 354850, 330455, 310620, 431490)
muns_name <- c('SP', 'Santos', 'RJ', 'BH', 'PA')

i = 1

list_plots <- list()

for(mun in muns){
  
  df_mun <- df_res %>% filter(ID_MN_RESI == mun)
  p <- ggplot(df_mun, aes(x = tplot)) +
    geom_line(aes(y = inc), color = "blue") +
    geom_line(aes(y = smooth_lambda*100, color = "red")) +  # Transform y2 to y1 scale
    scale_y_continuous(
      name = "inc",
      sec.axis = sec_axis(~ . /100, name = "lambda_s")  # Transform back
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(muns_name[i]) +
    geom_hline(yintercept = 100)
  
  print(p)
  
  list_plots[[i]] <- p
  
  i = i + 1
}

wrap_plots(list_plots, ncol = 2)
  


