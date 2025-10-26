# Calculating the index proposed by Rypdal and Sugihara, 2019, for UFS ##############

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
load('timeseries_dengue_2010_2023_SG.RData')
df <- df %>% ungroup()
df <- df %>% complete(anoepi, epiweek, SG_UF)
valid_53_week_years <- c(2015, 2020)
df <- df %>% filter(epiweek != 53 | anoepi %in% valid_53_week_years)
df <- df %>% filter(anoepi >= 2010) %>% filter(anoepi <= 2023)
df <- df %>% mutate(n = replace_na(n, 0))
df <- df %>% mutate(tplot = anoepi + (epiweek/52))
df <- df %>% select(anoepi, epiweek, tplot, SG_UF, n)
df_cases <- df
df_cases <- df_cases %>% drop_na()
rm(df)

# Carregando dados populacionais ###################################
load('pop_SG.RData')
df_pop <- df
df_pop <- df_pop %>% select('Ano', 'SG_UF', 'Total')
df_pop <- df_pop %>% mutate(Ano = as.numeric(Ano))
#df_pop <- df_pop %>% filter(Ano == 2016)
#df_pop <- df_pop %>% select(!Ano)
rm(df)
gc()

# Calculando incidencias semanais ##################################
df_cases <- df_cases %>% mutate(SG_UF = as.numeric(as.character(SG_UF)),
                                anoepi = as.numeric(as.character(anoepi)))
df_cases <- df_cases %>% left_join(df_pop, by = join_by(SG_UF == SG_UF,
                                                        anoepi == Ano))
df_cases <- df_cases %>% mutate(inc = (10**5)*(n/Total))
df_cases <- df_cases %>% select(anoepi, epiweek, tplot, SG_UF, inc)
df_cases <- df_cases %>% drop_na()
df_cases

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
sg_sel <- unique(df_cases$SG_UF)
#muns_sel <- c(310620)

i = 0

for(sg in sg_sel){
  print(i)
  df_aux <- df_cases %>% filter(SG_UF == sg)
  df_aux['inc_lambda'] <- rollapply(df_aux['inc'], width = ws1, FUN = function(w) lambda(w), align = "right", by.column = FALSE, fill = NA)
  df_aux['smooth_lambda'] <- rollmean(df_aux['inc_lambda'], k = 10, align = "right", fill = NA)
  df_res <- rbind(df_res, df_aux)
  i = i + 1
}

write_csv(df_res, 'timeseries_dengue_lambda_2010_2023_SG.csv')

# Plotting lambda and timeseries #######################################
# we are going to plot by region of the country

# lambda without smoothing #############################################


rgs <- c('NN', 'NE', 'SE', 'SL', 'CO')

for(name_rg in rgs){
  
  print(name_rg)
  if(name_rg == 'NN'){
    sg_sel <- c(11, 12, 13, 14, 15, 16, 17)
    sgs_name <- c('RO', 'AC', 'AM', 'RR', 'PA', 'AP', 'TO')
  }else if(name_rg == 'NE'){
    sg_sel <- c(21, 22, 23, 24, 25, 26, 27, 28, 29) 
    sgs_name <- c('MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA')
  }else if(name_rg == 'SE'){
    sg_sel <- c(31, 32, 33, 35)
    sgs_name <- c('MG', 'ES', 'RJ', 'SP')
  }else if(name_rg == 'SL'){
    sg_sel <- c(41, 42, 43)
    sgs_name <- c('PR', 'SC', 'RS')
  }else{
    sg_sel <- c(50, 51, 52, 53)
    sgs_name <- c('MS', 'MT', 'GO', 'DF')
  }
  
  i = 1
  
  list_plots <- list()
  
  for(sg in sg_sel){
    print(sg)
    df_mun <- df_res %>% filter(SG_UF == sg)
    p <- ggplot(df_mun, aes(x = tplot)) +
      geom_line(aes(y = inc), color = "blue") +
      geom_line(aes(y = inc_lambda*50, color = "red")) +  # Transform y2 to y1 scale
      scale_y_continuous(
        name = "inc",
        sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
      ) +
      theme_minimal() +
      theme(legend.position = "none") +
      ggtitle(sgs_name[i]) +
      geom_hline(yintercept = 50)
    
    #print(p)
    
    list_plots[[i]] <- p
    
    i = i + 1
  }
  
  print(wrap_plots(list_plots, ncol = 2))
  
}


# smoothed lambda ######################################################

rgs <- c('NN', 'NE', 'SE', 'SL', 'CO')

for(name_rg in rgs){
  
  print(name_rg)
  if(name_rg == 'NN'){
    sg_sel <- c(11, 12, 13, 14, 15, 16, 17)
    sgs_name <- c('RO', 'AC', 'AM', 'RR', 'PA', 'AP', 'TO')
  }else if(name_rg == 'NE'){
    sg_sel <- c(21, 22, 23, 24, 25, 26, 27, 28, 29) 
    sgs_name <- c('MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA')
  }else if(name_rg == 'SE'){
    sg_sel <- c(31, 32, 33, 35)
    sgs_name <- c('MG', 'ES', 'RJ', 'SP')
  }else if(name_rg == 'SL'){
    sg_sel <- c(41, 42, 43)
    sgs_name <- c('PR', 'SC', 'RS')
  }else{
    sg_sel <- c(50, 51, 52, 53)
    sgs_name <- c('MS', 'MT', 'GO', 'DF')
  }
  
  i = 1
  
  list_plots <- list()
  
  for(sg in sg_sel){
    print(sg)
    df_mun <- df_res %>% filter(SG_UF == sg)
    p <- ggplot(df_mun, aes(x = tplot)) +
      geom_line(aes(y = inc), color = "blue") +
      geom_line(aes(y = smooth_lambda*50, color = "red")) +  # Transform y2 to y1 scale
      scale_y_continuous(
        name = "inc",
        sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
      ) +
      theme_minimal() +
      theme(legend.position = "none") +
      ggtitle(sgs_name[i]) +
      geom_hline(yintercept = 50)
    
    #print(p)
    
    list_plots[[i]] <- p
    
    i = i + 1
  }
  
  print(wrap_plots(list_plots, ncol = 2))
  
}

















