library(tidyverse)
library(car)
library(RColorBrewer)


file_name = 'cf_all_clean.csv'
file_name = 'md_all_clean.csv'
df_main = data.table::fread(file_name, 
                            colClasses = c('mbr_zip_5_cd' = 'character',
                                           'age_band_cd' = 'character',
                                           'ndc' = 'character',
                                           'z_patid' = 'character'))

df_cost = df_main[, c('z_patid', 'mbr_state','sex', 
                      'age_band_cd', 'inp_cost', 'op_cost', 'ph_cost', 'pharmacy_cost')]
df_cost[is.na(df_cost)] = 0
df_cost$total_cost = rowSums(df_cost[, c('inp_cost', 'op_cost', 'ph_cost', 'pharmacy_cost')])

df_inp = df_main[, c('z_patid', 'mbr_state','sex', 
                     'inp_cost', 'inp_freq', 'inp_avg_lag')]
df_inp = df_inp %>% drop_na() %>% filter(inp_cost>0)

kmdat = df_inp[,c('inp_cost', 'inp_freq', 'inp_avg_lag')]
kmdat$inp_cost = log(kmdat$inp_cost)
kmdat = kmdat %>% data.frame()
for(i in 1:ncol(kmdat)){
  kmdat[,i] = (kmdat[,i] - min(kmdat[,i]))/max(kmdat[,i])
}

km = kmeans(kmdat, 5, iter.max = 300)

df_inp$clusters = km$cluster
colors = brewer.pal(n = 5, 'Set1')

scatter3d(
  x = log(df_inp$inp_avg_lag+1),
  y = log(df_inp$inp_cost), 
  z = log(df_inp$inp_freq+1),
  xlab = 'Lag between Claims', ylab = 'Cost', zlab = 'Frequency of Claims',
  groups = factor(df_inp$clusters),
  ellipsoid = TRUE, surface.col = colors,
  surface = FALSE, revolutions =  3
)

df_summary = df_inp %>% group_by(clusters) %>% 
summarise(members = n(),
          cost = median(inp_cost),
          freq = median(inp_freq),
          avg_lag = median(inp_avg_lag))



