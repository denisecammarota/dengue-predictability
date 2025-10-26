# Calculating activity thresholds for muns ###################
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

# Using the selected municipalities ###########################
muns_sel <- c(310620, 520870, 530010, 230440, 
              330455, 520140, 350950, 354980, 
              355030, 354340)
muns_names <- c(
  "Belo Horizonte",
  "Goiânia",
  "Brasília",
  "Fortaleza",
  "Rio de Janeiro",
  "Anápolis",
  "Campinas",
  "São José do Rio Preto",
  "São Paulo",
  "Ribeirão Preto"
)

# Loading incidence data ####################################
df_res <- read.csv('timeseries_dengue_lambda_2010_2023.csv')
df_res <- df_res %>% filter(ID_MN_RESI %in% muns_sel)

## Looking at seasonal incidence #########################
i = 1
list_plots <- list()
for(mun in muns_sel){
  df_plot <- df_res %>% filter(ID_MN_RESI == mun)
  df_plot <- df_plot %>% mutate(anoepi = factor(anoepi))
  p <- ggplot(df_plot, aes(x = epiweek, y = inc, color = anoepi)) + 
    geom_line() +
    geom_point() +
    xlab('Sem Epidemiológica') +
    ylab('Incidencia') +
    theme_minimal() +
    ggtitle(muns_names[i])
  print(p)
  list_plots[[i]] <- p
  i = i + 1
}

wrap_plots(list_plots, ncol = 2)
wrap_plots(list_plots) + plot_layout(guides = "collect")

# Selecting seasons ######################################
df_res <- df_res %>% 
  mutate(season = ifelse(epiweek <= 40, 
                         paste0(anoepi - 1, '/', anoepi), 
                         paste0(anoepi, '/', anoepi + 1)))
df_res <- df_res %>% filter(!season %in% c('2009/2010', '2023/2024'))
df_res <- df_res %>% filter(epiweek != 53)

df_sel <-  df_res
df_sel <- df_sel %>%
  mutate(keep = 1) %>%
  mutate(keep = case_when(
    ID_MN_RESI == 530010 & season == '2021/2022' ~ 0, 
    ID_MN_RESI == 230440 & season == '2011/2012' ~ 0,
    ID_MN_RESI == 330455 & season == '2010/2011' ~ 0,
    ID_MN_RESI == 350950 & season == '2013/2014' ~ 0,
    ID_MN_RESI == 350950 & season == '2014/2015' ~ 0,
    ID_MN_RESI == 354980 & season == '2018/2019' ~ 0,
    ID_MN_RESI == 355030 & season == '2014/2015' ~ 0,
    ID_MN_RESI == 355030 & season == '2013/2014' ~ 0,
    ID_MN_RESI == 354340 & season == '2015/2016' ~ 0,
    T ~ keep
  ))

df_sel <- df_sel %>% filter(keep == 1)
df_sel <- df_sel %>% select(!keep)

# Calculate MEM #######################################
list_pre <- list()
i = 1
for(mun in muns_sel){
  print(mun)
  df_sp <- df_sel %>% filter(ID_MN_RESI == mun) %>% 
    select(season, epiweek, inc) %>%
    pivot_wider(names_from = season, values_from = inc)
  list_epiweek <- df_sp$epiweek
  df_sp <- df_sp %>% select(!epiweek)
  dengue.memmodel <- memmodel(df_sp, i.season = 30, i.method = 1)
  plot(dengue.memmodel)
  pre_aux <- as.numeric(dengue.memmodel$pre.post.intervals["pre.i", 3])
  list_pre[[i]] <- pre_aux[1]
  i = i + 1
}
list_pre
list_pre <- unlist(list_pre)
list_pre

# Previous lambda tells us about the pre-epidemic threshold #########
df_pre <- data.frame(ID_MN_RESI = muns_sel,
                     pre_thr = list_pre)
df_res <- df_res %>% left_join(df_pre, by = join_by(ID_MN_RESI == ID_MN_RESI))
df_res <- df_res %>% mutate(inc_bin = ifelse(inc > pre_thr, 1, 0)) %>%
  mutate(lambda_bin = ifelse(inc_lambda > 1, 1, 0)) %>%
  mutate(smooth_bin = ifelse(smooth_lambda > 1, 1, 0))

## Calculando accuracy com regular lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>%
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP + TN
    den <- TP + TN + FN + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Accuracy') +
  labs(color = 'Municipality') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")


## Calculando accuracy com smoothed lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(smooth_bin = lead(smooth_bin, -lag)) %>%
      group_by(inc_bin, smooth_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP + TN
    den <- TP + TN + FN + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Accuracy') +
  labs(color = 'Municipality (smooth)') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")



## Calculando sensitivity com regular lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP
    den <- TP + FN
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Sensitivity') +
  labs(color = 'Municipality') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")
  


## Calculando sensitivity com smoothed lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1
#muns_sel <- muns_sel[1]

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(smooth_bin = lead(smooth_bin, -lag)) %>% 
      group_by(inc_bin, smooth_bin) %>% 
      summarise(n = n())
    print(df_aux) 
    TP <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 0)
    #print(TP)
    #print(FN)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP
    den <- TP + FN 
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Sensitivity') +
  labs(color = 'Municipality (smooth)') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")

## Calculando precision com regular lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP
    den <- TP + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Precision') +
  labs(color = 'Municipality') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")



## Calculando precision com smoothed lambda #################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(muns_sel), ncol = length(lags))
i = 1
#muns_sel <- muns_sel[1]

for(mun in muns_sel){
  print(mun)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(ID_MN_RESI == mun) %>%
      mutate(smooth_bin = lead(smooth_bin, -lag)) %>% 
      group_by(inc_bin, smooth_bin) %>% 
      summarise(n = n())
    print(df_aux) 
    TP <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & smooth_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & smooth_bin == 0)
    #print(TP)
    #print(FN)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP
    den <- TP + FP 
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- muns_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  xlab('Lag') +
  ylab('Precision') +
  labs(color = 'Municipality (smooth)') +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black")


## Final plots ##############################################
muns <- muns_sel

i = 1

list_plots <- list()

for(mun in muns){
  
  df_mun <- df_res %>% filter(ID_MN_RESI == mun)
  thr <- unique(df_mun$pre_thr)
  p <- ggplot(df_mun, aes(x = tplot)) +
    geom_line(aes(y = inc), color = "blue") +
    geom_line(aes(y = inc_lambda*500, color = "red")) +  # Transform y2 to y1 scale
    scale_y_continuous(
      name = "inc",
      sec.axis = sec_axis(~ . /500, name = "lambda")  # Transform back
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(mun) +
    geom_hline(yintercept = 500, color = 'red') +
    geom_hline(yintercept = thr, color = 'blue')
  
  print(p)
  
  list_plots[[i]] <- p
  
  i = i + 1
}

wrap_plots(list_plots, ncol = 2)

muns <- muns_sel

i = 1

list_plots <- list()

for(mun in muns){
  
  df_mun <- df_res %>% filter(ID_MN_RESI == mun)
  thr <- unique(df_mun$pre_thr)
  p <- ggplot(df_mun, aes(x = tplot)) +
    geom_line(aes(y = inc), color = "blue") +
    geom_line(aes(y = smooth_lambda*500, color = "red")) +  # Transform y2 to y1 scale
    scale_y_continuous(
      name = "inc",
      sec.axis = sec_axis(~ . /500, name = "lambda")  # Transform back
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(mun) +
    geom_hline(yintercept = 500, color = 'red') +
    geom_hline(yintercept = thr, color = 'blue')
  
  print(p)
  
  list_plots[[i]] <- p
  
  i = i + 1
}

wrap_plots(list_plots, ncol = 2)

## Calculate per epidemiological week #######################

### Accuracy ################################################

## Regular Lambda

lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1
mun <- 355030

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP + TN
    den <- TP + TN + FN + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Accuracy",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Accuracy')
  


## Smoothed Lambda
lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(smooth_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP + TN
    den <- TP + TN + FN + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Accuracy (s)",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Accuracy  (s)')

### Sensitivity ##############################################

## Regular Lambda

lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP 
    den <- TP + FN 
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Sensitivity",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Sensitivity')


## Smoothed Lambda
lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(smooth_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP 
    den <- TP + FN 
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Sensitivity (s)",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Sensitivity (s)')

### Precision ##############################################

## Regular Lambda

lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(lambda_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP 
    den <- TP + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Precision",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Sensitivity')


## Smoothed Lambda
lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      #filter(ID_MN_RESI == mun) %>%
      mutate(lambda_bin = lead(smooth_bin, -lag)) %>%
      filter(epiweek == week) %>% 
      group_by(inc_bin, lambda_bin) %>% 
      summarise(n = n())
    TP <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 1)
    FP <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 1)
    TN <- df_aux %>% filter(inc_bin == 0 & lambda_bin == 0)
    FN <- df_aux %>% filter(inc_bin == 1 & lambda_bin == 0)
    if(TP %>% nrow() == 0){
      TP <- 0
    }else{
      TP <- TP$n
    }
    if(FP %>% nrow() == 0){
      FP <- 0
    }else{
      FP <- FP$n
    }
    if(TN %>% nrow() == 0){
      TN <- 0
    }else{
      TN <- TN$n
    }
    if(FN %>% nrow() == 0){
      FN <- 0
    }else{
      FN <- FN$n
    }
    #TP <- TP$n
    #FP <- FP$n
    #TN <- TN$n
    #FN <- FN$n
    num <- TP 
    den <- TP + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- weeks
rownames(lag_mat) <- lags
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 0.5,         # set your desired middle point
    name = "Precision (s)",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Sensitivity (s)')
