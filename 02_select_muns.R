library(tidyverse)
library(geobr)

# Select the municipalities with the most cases ###################
load('timeseries_dengue_2010_2023.RData')
df_grouped <- df %>% group_by(ID_MN_RESI) %>% summarise(cases = sum(n))
df_grouped
rm(df)
gc()

## Loading population data #########################################

load('pop_muns.RData')
df_pop <- df
df_pop <- df_pop %>% select('Ano', 'Mun', 'Total')
df_pop <- df_pop %>% filter(Ano == 2016)
df_pop <- df_pop %>% select(!Ano)
rm(df)
gc()


## Calculating incidence ############################################
df <- df_grouped
df <- df %>% left_join(df_pop, by = join_by(ID_MN_RESI == Mun))
df <- df %>% mutate(inc_100k = (10**5)*(cases/Total))
df

## Selecting municipalities ########################################
df <- df %>% arrange(desc(cases))
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

rm(df, df_grouped)

# Plotting municipalities ###############################################
load('timeseries_dengue_2010_2023.RData')
df <- df %>% ungroup()
df <- df %>% complete(anoepi, epiweek, ID_MN_RESI)
df <- df %>% mutate(n = replace_na(n, 0))
df <- df %>% mutate(tplot = anoepi + (epiweek/52))
df <- df %>% filter(ID_MN_RESI %in% muns_sel)
df <- df %>% select(ID_MN_RESI, tplot, n)

ggplot(df, aes(x = tplot, y = n, color = ID_MN_RESI)) + 
  geom_line() +
  theme_bw()


