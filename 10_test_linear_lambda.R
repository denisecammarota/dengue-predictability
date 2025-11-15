# Simple test for RJ and SP (Rypdal and Sugihara, 2019, fig 1) #######
# SP = 355030
# RJ = 330455

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

# Doing linear test of dependence in SP and RJ ###################################

# Loading epidemiological data ########################################
df_res <- read.csv('timeseries_dengue_lambda_2010_2023.csv')
df_res <- df_res %>% filter(ID_MN_RESI %in% c(330455, 355030))
df_res <- df_res %>% mutate(inc_lambda_bin = ifelse(inc_lambda >= 1, 1, 0),
                            smooth_lambda_bin = ifelse(smooth_lambda >= 1, 1, 0))
df_sp <- df_res %>% filter(ID_MN_RESI == 355030) %>% filter(anoepi >= 2011)
df_rj <- df_res %>% filter(ID_MN_RESI == 330455) %>% filter(anoepi >= 2011)

# Loading MEM data #################################################
df_mem <- read.csv('mem_municipalities.csv')
df_mem <- df_mem %>% filter(ID_MN_RESI %in% c(330455, 355030))
pre_sp <- 0.5*(4.033144 + 3.320181)
pre_rj <- 0.5*(7.628334 + 6.376419)
df_sp <- df_sp %>% mutate(inc_bin = ifelse(inc >= pre_sp, 1, 0))
df_rj <- df_rj %>% mutate(inc_bin = ifelse(inc >= pre_rj, 1, 0))

# Plotting for SP and RJ ################################################

## SP ####################################################################

ggplot(df_sp, aes(x = tplot)) +
  geom_line(aes(y = inc), color = "blue") +
  geom_line(aes(y = inc_lambda*50, color = "red")) +  # Transform y2 to y1 scale
  scale_y_continuous(
    name = "inc",
    sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle('São Paulo') +
  geom_hline(yintercept = 50)

ggplot(df_sp, aes(x = tplot)) +
  geom_line(aes(y = inc), color = "blue") +
  geom_line(aes(y = smooth_lambda*50, color = "red")) +  # Transform y2 to y1 scale
  scale_y_continuous(
    name = "inc",
    sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle('São Paulo') +
  geom_hline(yintercept = 50)


## RJ #######################################################################

ggplot(df_rj, aes(x = tplot)) +
  geom_line(aes(y = inc), color = "blue") +
  geom_line(aes(y = smooth_lambda*50, color = "red")) +  # Transform y2 to y1 scale
  scale_y_continuous(
    name = "inc",
    sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle('Rio de Janeiro') +
  geom_hline(yintercept = 50)

# Taking SP and doing the full code #########################################
p <- ggplot(df_sp, aes(x = tplot)) +
  geom_line(aes(y = inc), color = "blue") +
  geom_line(aes(y = smooth_lambda*50, color = "red")) +  # Transform y2 to y1 scale
  scale_y_continuous(
    name = "inc",
    sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle('São Paulo') +
  geom_hline(yintercept = 50) +
  geom_hline(yintercept = pre_sp)
p

colnames(df_sp)
df_sp <- df_sp %>% select(tplot, inc, inc_bin, smooth_lambda, smooth_lambda_bin)
mean(df_sp$smooth_lambda[df_sp$smooth_lambda < 1], na.rm = TRUE)


df_periods <- df_sp %>%
  arrange(tplot) %>%
  mutate(flag = smooth_lambda < 1,
         # start a new group every time flag switches from FALSE → TRUE
         grp = cumsum(flag & !lag(flag, default = FALSE))) %>%
  filter(flag) %>%  # keep only periods where lambda < 1
  group_by(grp) %>%
  summarise(
    start = first(tplot),
    end   = last(tplot),
    n     = n(),
    mean_lambda = mean(smooth_lambda, na.rm = TRUE)
  )
df_epi_periods <- df_sp %>%
  arrange(tplot) %>%
  mutate(flag = smooth_lambda >= 1,
         # start a new group every time flag switches from FALSE → TRUE
         grp = cumsum(flag & !lag(flag, default = FALSE))) %>%
  filter(flag) %>%  # keep only periods where lambda < 1
  group_by(grp) %>%
  summarise(
    start = first(tplot),
    end   = last(tplot),
    n     = n(),
    sum_inc = sum(inc, na.rm = TRUE)
  )

df_periods <- df_periods %>% filter(grp != 1) %>% filter(grp != 13)
df_epi_periods <- df_epi_periods %>% filter(grp != 1) %>% filter(grp != 13)
df_total <- df_periods %>% left_join(df_epi_periods, by = join_by(grp))
model <- lm(sum_inc ~ mean_lambda, data = df_total)
coef_vals <- coef(model)

ggplot(df_total, aes(x = mean_lambda, y = sum_inc)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text",
           x = Inf, y = -Inf, label = eq_text,
           hjust = 1.1, vjust = -1.1, size = 4) +
  theme_bw() +
  xlab("Total Incidence") +
  ylab("Mean Lambda (lambda < 1)")


colnames(df_rj)
df_rj <- df_rj %>% select(tplot, inc, inc_bin, smooth_lambda, smooth_lambda_bin)

p <- ggplot(df_rj, aes(x = tplot)) +
  geom_line(aes(y = inc), color = "blue") +
  geom_line(aes(y = smooth_lambda*50, color = "red")) +  # Transform y2 to y1 scale
  scale_y_continuous(
    name = "inc",
    sec.axis = sec_axis(~ . /50, name = "lambda_s")  # Transform back
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle('São Paulo') +
  geom_hline(yintercept = 50) +
  geom_hline(yintercept = pre_sp)
p


df_periods <- df_rj %>%
  arrange(tplot) %>%
  mutate(flag = smooth_lambda < 1,
         # start a new group every time flag switches from FALSE → TRUE
         grp = cumsum(flag & !lag(flag, default = FALSE))) %>%
  filter(flag) %>%  # keep only periods where lambda < 1
  group_by(grp) %>%
  summarise(
    start = first(tplot),
    end   = last(tplot),
    n     = n(),
    mean_lambda = mean(smooth_lambda, na.rm = TRUE)
  )
df_epi_periods <- df_rj %>%
  arrange(tplot) %>%
  mutate(flag = smooth_lambda >= 1,
         # start a new group every time flag switches from FALSE → TRUE
         grp = cumsum(flag & !lag(flag, default = FALSE))) %>%
  filter(flag) %>%  # keep only periods where lambda < 1
  group_by(grp) %>%
  summarise(
    start = first(tplot),
    end   = last(tplot),
    n     = n(),
    sum_inc = sum(inc, na.rm = TRUE)
  )

df_periods <- df_periods %>% filter(grp != 1) %>% filter(grp != 13)
df_epi_periods <- df_epi_periods %>% filter(grp != 1) %>% filter(grp != 13)
df_total <- df_periods %>% left_join(df_epi_periods, by = join_by(grp))
model <- lm(sum_inc ~ mean_lambda, data = df_total)
coef_vals <- coef(model)

eq_text <- paste0(
  "y = ", round(coef_vals[1], 3),
  " + ", round(coef_vals[2], 3), "x"
)

ggplot(df_total, aes(x = mean_lambda, y = sum_inc)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  #annotate("text",
  #         x = Inf, y = -Inf, label = eq_text,
  #         hjust = 1.1, vjust = -1.1, size = 4) +
  theme_bw() +
  xlab("Total Incidence") +
  ylab("Mean Lambda (lambda < 1)")






