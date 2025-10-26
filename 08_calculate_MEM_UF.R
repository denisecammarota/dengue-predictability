# Calculating activity thresholds for SG_UFs ###################
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

# Loading incidence data ########################################
df_res <- read.csv('timeseries_dengue_lambda_2010_2023_SG.csv')

# Extracting manually extreme seasons ###########################

## Visualizing seasons ##########################################
sg_sel <- unique(df_res$SG_UF)
sgs_names <- c('RO', 'AC', 'AM', 'RR', 'PA', 'AP', 'TO',
               'MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA',
               'MG', 'ES', 'RJ', 'SP',
               'PR', 'SC', 'RS',
               'MS', 'MT', 'GO', 'DF')

i = 1
list_plots <- list()
for(sg in sg_sel){
  df_plot <- df_res %>% filter(SG_UF == sg)
  df_plot <- df_plot %>% mutate(anoepi = factor(anoepi))
  p <- ggplot(df_plot, aes(x = epiweek, y = inc, color = anoepi)) + 
    geom_line(size = 1.5) +
    geom_point(size = 1.5) +
    xlab('Sem Epidemiológica') +
    ylab('Incidencia') +
    theme_minimal() +
    ggtitle(sgs_names[i])
  print(p)
  list_plots[[i]] <- p
  i = i + 1
}

## Shifting seasons for doing MEM ########################
df_res <- df_res %>% 
  mutate(season = ifelse(epiweek <= 40, 
                         paste0(anoepi - 1, '/', anoepi), 
                         paste0(anoepi, '/', anoepi + 1)))
df_res <- df_res %>% filter(!season %in% c('2009/2010', '2023/2024'))
df_res <- df_res %>% filter(epiweek != 53)

## Eliminating seasons with extreme incidence ##############
df_sum <- df_res %>% group_by(SG_UF, season) %>% summarise(inc_max = max(inc))
n_seasons <- length(unique(df_res$season))
df_sum <- df_sum %>% group_by(SG_UF) %>% summarise(log_geom_mean = mean(log(inc_max)),
                                                   log_geom_upper = mean(log(inc_max)) + ((qt(0.975, n_seasons - 1)*sd(log(inc_max)))/(sqrt(n_seasons))),
                                                   log_geom_lower = mean(log(inc_max)) - ((qt(0.975, n_seasons - 1)*sd(log(inc_max)))/(sqrt(n_seasons))))
df_sum <- df_sum %>% mutate(geom_mean = exp(log_geom_mean), 
                               geom_upper = exp(log_geom_upper),
                               geom_lower = exp(log_geom_lower))
df_sum <- df_sum %>% select(!log_geom_mean) %>% select(!log_geom_upper) %>% select(!log_geom_lower)
df_sum

df_seasons <- df_res %>% group_by(SG_UF, season) %>% summarise(inc_max = max(inc))
df_seasons <- df_seasons %>% left_join(df_sum, by = join_by(SG_UF))
df_seasons <- df_seasons %>% mutate(keep = 1) %>%
  mutate(keep = case_when(
    inc_max < geom_lower ~ 0,
    inc_max > geom_upper ~ 0,
    T ~ 1
  )) %>%
  mutate(keep = case_when(
    SG_UF == 51 & season %in% c('2013/2014', '2010/2011') ~ 1,
    SG_UF == 42 & season %in% c('2016/2017') ~ 1,
    SG_UF == 31 & season %in% c('2020/2021', '2016/2017') ~ 1,
    SG_UF == 26 & season %in% c('2022/2023') ~ 1,
    SG_UF == 24 & season %in% c('2020/2021') ~ 1,
    SG_UF == 17 & season %in% c('2017/2018') ~ 1,
    SG_UF == 15 & season %in% c('2018/2019') ~ 1,
    SG_UF == 12 & season %in% c('2016/2017') ~ 1,
    T ~ keep
  ))

df_seasons <- df_seasons %>% select(SG_UF, season, keep)

rm(df_sum)
gc()

df_res <- df_res %>% left_join(df_seasons, by = join_by(SG_UF, season))
df_res <- df_res %>% filter(keep == 1)


i = 1
list_plots <- list()
for(sg in sg_sel){
  df_plot <- df_res %>% filter(SG_UF == sg)
  df_plot <- df_plot %>% mutate(anoepi = factor(anoepi))
  p <- ggplot(df_plot, aes(x = epiweek, y = inc, color = factor(season))) + 
    geom_line(size = 1.5) +
    geom_point(size = 1.5) +
    xlab('Sem Epidemiológica') +
    ylab('Incidencia') +
    theme_minimal() +
    ggtitle(sgs_names[i])
  print(p)
  list_plots[[i]] <- p
  i = i + 1
}


# Calculate MEM ###################################################
list_pre <- list()
list_pos <- list()
list_start <- list()
list_duration <- list()
i = 1
for(sg in sg_sel){
  print(sg)
  df_sp <- df_res %>% filter(SG_UF == sg) %>% 
    select(season, epiweek, inc) %>%
    pivot_wider(names_from = season, values_from = inc)
  list_epiweek <- df_sp$epiweek
  df_sp <- df_sp %>% select(!epiweek)
  dengue.memmodel <- memmodel(df_sp, i.season = 30, i.method = 1)
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

df_summary <- data.frame(sg_uf = sg_sel, 
                         pre_thr = list_pre,
                         pos_thr = list_pos,
                         start_time = list_start,
                         dur_time = list_duration)

# Making maps of quantities ###########################

map_br <- read_state()
map_br <- map_br %>% left_join(df_summary, by = join_by(code_state == sg_uf))

ggplot(map_br, aes(fill = pre_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Pre-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 5,         # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = pos_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Pre-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 5,         # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = start_time)) + geom_sf() + theme_map() +
  labs(fill = 'Start Week') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 15,         # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = dur_time)) + geom_sf() + theme_map() +
  labs(fill = 'Start Week') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = 24,         # set your desired middle point
    guide = "colorbar"
  )

# Measuring accuracy and precision of lambda as EW ###########

df_res <- df_res %>% select(!keep)
df_res <- df_res %>% left_join(df_summary, by = join_by(SG_UF == sg_uf))
df_res <- df_res %>% select(!start_time) %>% select(!dur_time) %>%
  select(!pos_thr)
df_res <- df_res %>% mutate(inc_bin = ifelse(inc > pre_thr, 1, 0)) %>%
  mutate(lambda_bin = ifelse(inc_lambda > 1, 1, 0)) %>%
  mutate(smooth_bin = ifelse(smooth_lambda > 1, 1, 0))
df_res

## Overall accuracy ############################################

### Regular lambda #############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF") +
  xlab("Lag") +
  ylab("Accuracy") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

### Smoothed lambda ############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF (Smooth)") +
  xlab("Lag") +
  ylab("Accuracy") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

## Overall precision ###########################################


### Regular lambda #############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF") +
  xlab("Lag") +
  ylab("Precision") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

### Smoothed lambda ############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
    num <- TP
    den <- TP + FP
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF (Smooth)") +
  xlab("Lag") +
  ylab("Precision") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

## Overall recall ###########################################


### Regular lambda #############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF") +
  xlab("Lag") +
  ylab("Recall") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

### Smoothed lambda ############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
    num <- TP
    den <- TP + FN
    print(paste0(i,j))
    lag_mat[i, j] <- num/den
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF (Smooth)") +
  xlab("Lag") +
  ylab("Recall") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )


## Overall F1-score ###########################################


### Regular lambda #############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
    prec <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    print(paste0(i,j))
    lag_mat[i, j] <- 2*((prec**-1 + recall**-1)**-1)
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF") +
  xlab("Lag") +
  ylab("F1-Score") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

### Smoothed lambda ############################################
lags <- seq(-26, 0, 1)
lag_mat <- matrix(0, nrow = length(sg_sel), ncol = length(lags))
i = 1

for(sg in sg_sel){
  print(sg)
  j = 1
  for(lag in lags){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
    prec <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    print(paste0(i,j))
    lag_mat[i, j] <- 2*((prec**-1 + recall**-1)**-1)
    j = j + 1
  }
  i = i + 1
}

colnames(lag_mat) <- lags
rownames(lag_mat) <- sg_sel
df_plot <- melt(lag_mat)
df_plot

ggplot(df_plot, aes(x = Var2, y = value, color = factor(Var1), group = Var1)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "turbo", name = "UF (Smooth)") +
  xlab("Lag") +
  ylab("F1-Score") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )







