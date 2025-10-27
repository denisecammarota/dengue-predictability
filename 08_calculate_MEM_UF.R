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

rm(p, list_plots, df_plot)

## Shifting seasons for doing MEM ########################
df_res <- df_res %>% 
  mutate(season = ifelse(epiweek <= 40, 
                         paste0(anoepi - 1, '/', anoepi), 
                         paste0(anoepi, '/', anoepi + 1)))
df_res <- df_res %>% filter(!season %in% c('2009/2010', '2023/2024'))
df_res <- df_res %>% filter(epiweek != 53)
df_res

df_res <- df_res %>% mutate(epiweek_aux = ifelse(epiweek >= 40, epiweek - 40 + 1, epiweek + 52 - 40 + 1))

i = 1
list_plots <- list()
for(sg in sg_sel){
  df_plot <- df_res %>% filter(SG_UF == sg)
  df_plot <- df_plot %>% mutate(anoepi = factor(anoepi))
  p <- ggplot(df_plot, aes(x = epiweek_aux, y = inc, color = factor(season))) + 
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

rm(p, list_plots, df_plot)


## Eliminating seasons with extreme incidence ########################

# Calculating peaks and their geom mean and distance
df_sum <- df_res %>% group_by(SG_UF, season) %>% summarise(inc_max = max(inc))
df_mean_dist <- df_sum %>% group_by(SG_UF) %>% summarise(geom_mean = exp(mean(log(inc_max))))
df_sum <- df_sum %>% left_join(df_mean_dist, by = join_by(SG_UF))
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
df_tmp_seasons <- df_sum %>% group_by(SG_UF) %>% 
  summarise(seasons_keep = sum(keep)) %>%
  filter(seasons_keep < 5) %>%
  mutate(seasons_need = 5 - seasons_keep)
df_tmp_seasons$SG_UF
df_tmp_seasons$seasons_need

# Adding seasons that are needed 
for(sg in unique(df_tmp_seasons$SG_UF)){
  n_seasons_need <- df_tmp_seasons %>% filter(SG_UF == sg)
  n_seasons_need <- n_seasons_need$seasons_need[1]
  df_tmp_sg <- df_sum %>% filter(SG_UF == sg) %>% filter(keep == 0)
  df_tmp_sg <- df_tmp_sg %>% arrange(abs(geom_dist - perc_10_90[1]), abs(geom_dist - perc_10_90[2]))
  df_tmp_sg <- df_tmp_sg %>% slice_head(n = n_seasons_need) %>% mutate(keep = 1)
  df_sum <- rbind(df_sum, df_tmp_sg)
}

rm(n_seasons_need, df_tmp_sg, df_tmp_seasons)
gc()

df_sum <- df_sum %>% group_by(SG_UF, season) %>% summarise(keep = max(keep))
df_mem <- df_res %>% left_join(df_sum, by = join_by(SG_UF, season)) 
df_mem <- df_mem %>% filter(keep == 1)
df_mem <- df_mem %>% select(anoepi, epiweek, tplot, season, epiweek_aux, 
                            SG_UF, inc, inc_lambda, smooth_lambda
                            )

i = 1
list_plots <- list()
for(sg in sg_sel){
  df_plot <- df_mem %>% filter(SG_UF == sg)
  df_plot <- df_plot %>% mutate(anoepi = factor(anoepi))
  p <- ggplot(df_plot, aes(x = epiweek_aux, y = inc, color = factor(season))) + 
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

rm(p, list_plots, df_plot)


# Calculate MEM ###################################################
list_pre <- list()
list_pos <- list()
list_start <- list()
list_duration <- list()
i = 1
for(sg in sg_sel){
  print(sg)
  df_sp <- df_mem %>% filter(SG_UF == sg) %>% 
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

df_summary <- data.frame(sg_uf = sg_sel, 
                         pre_thr = list_pre,
                         pos_thr = list_pos,
                         start_time = list_start,
                         dur_time = list_duration)
df_summary <- df_summary %>% mutate(start_time_real = ifelse(start_time <= 13, start_time + 39, start_time - 13))

# Making maps of quantities ###########################

map_br <- read_state()
map_br <- map_br %>% left_join(df_summary, by = join_by(code_state == sg_uf))

ggplot(map_br, aes(fill = pre_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Pre-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$pre_thr),  # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = pos_thr)) + geom_sf() + theme_map() +
  labs(fill = 'Post-epidemic Threshold') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$pos_thr),         # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = start_time_real)) + geom_sf() + theme_map() +
  labs(fill = 'Start Week') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$start_time_real),   # set your desired middle point
    guide = "colorbar"
  )

ggplot(map_br, aes(fill = dur_time)) + geom_sf() + theme_map() +
  labs(fill = 'Duration') +
  scale_fill_gradient2(
    low = "#313695",      # dark blue for low
    mid = "#FFFFBF",      # light yellow at midpoint
    high = "#A50026",     # dark red for high
    midpoint = mean(map_br$dur_time), # set your desired middle point
    guide = "colorbar"
  )

# Measuring accuracy and precision of lambda as EW ###########

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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Accuracy (mean)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Accuracy (mean, smooth)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Precision (mean)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Precision (mean, smooth)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Recall (mean)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("Recall (mean, smooth)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("F1-Score (mean)") +
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

df_plot <- df_plot %>% group_by(Var2) %>% summarise(
  value = mean(value)
)

ggplot(df_plot, aes(x = Var2, y = value)) +
  geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  xlab("Lag") +
  ylab("F1-Score (mean, smooth)") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm")
  )

## Calculate per epidemiological week #######################

### Accuracy ################################################

## Regular Lambda

lags <- seq(-26, 0, 1)
weeks <- sort(unique(df_res$epiweek))
lag_mat <- matrix(0, nrow = length(lags), ncol = length(weeks))
i = 1
sg <- 35

for(lag in lags){
  print(lag)
  j = 1
  for(week in weeks){
    df_aux <- df_res %>% 
      filter(SG_UF == sg) %>%
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
      filter(SG_UF == 35) %>%
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
      filter(SG_UF == 35) %>%
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
  labs(fill = 'Precision')


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
      filter(SG_UF == 35) %>%
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
  labs(fill = 'Precision (s)')


### Recall ##############################################

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
      filter(SG_UF == 35) %>%
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
    name = "Recall",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Recall')


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
      filter(SG_UF == 35) %>%
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
    name = "Recall (s)",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'Recall (s)')

### F1-Score ##############################################

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
      filter(SG_UF == 35) %>%
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
    prec <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    print(paste0(i,j))
    lag_mat[i, j] <- 2*((prec**-1 + recall**-1)**-1)
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
    name = "F1-Score",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'F1-Score')


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
      filter(SG_UF == 35) %>%
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
    prec <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    print(paste0(i,j))
    lag_mat[i, j] <- 2*((prec**-1 + recall**-1)**-1)
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
    name = "F1-Score (s)",
    guide = "colorbar"
  ) +
  theme_bw() +
  xlab('Epidemiological Week') +
  ylab('Lags') +
  labs(fill = 'F1-Score (s)')




