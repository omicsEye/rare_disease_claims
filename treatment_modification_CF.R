library(tidyverse)
library(RODBC)
library(scales)
library(usmap)

report_directory = paste(getwd(), "/", 'reports/treatment_modification/', sep = "")
cbp <- c('#e69f00', '#0072b2','#d55e00', '#009e73', 
         '#f0e442', '#56b4e9')
patient_list = read.csv('reports/prevalence/overall/patient_list.csv')

ch = odbcConnect('vertica',
                 uid = 'XXXXX',
                 pwd = 'XXX')

## Cystic Fibrosis
patients_raw = patient_list %>% 
  filter(disease =='Cystic_Fibrosis') %>% 
  select(z_patid) %>% 
  unlist()

patients_raw_query = paste("'",patients_raw,"'", sep="")


first_dates = sqlQuery(ch, paste("SELECT rs.z_patid, min(rs.fill_dt) as first_date
from HCCI_2.RX_SDDV2 rs
where rs.yr>2015 AND rs.ndc like '511670%'
and rs.ndc not in ('51167010001' ,'51167010003') AND rs.z_patid in (", 
                                 paste(patients_raw_query, collapse = ','), ")
GROUP by 1"))

patients = paste("'",first_dates %>% select(z_patid) %>% unlist(),"'", sep="")

inp = sqlQuery(ch, paste("SELECT inp.z_patid, inp.fst_dt, 
SUM(inp.calc_allwd) as cost, 'Inpatient' as claim_type
from HCCI_2.INP_SDDV2 inp
where inp.yr > 2015 AND inp.z_patid in (", 
                         paste(patients, collapse = ','), ")
                           group by 1,2"))

op = sqlQuery(ch, paste("SELECT op.z_patid, op.fst_dt,
SUM(op.calc_allwd) as cost, 'Outpatient' as claim_type
from HCCI_2.OP_SDDV2 op
where op.yr>2015 AND op.z_patid in (", 
                        paste(patients, collapse = ','), ")
group by 1,2"))

ph = sqlQuery(ch, paste("SELECT ph.z_patid, ph.fst_dt,
SUM(ph.calc_allwd) as cost, 'Physician' as claim_type
from HCCI_2.PHYS_SDDV2 ph 
where ph.yr > 2015 AND ph.z_patid in (", 
                        paste(patients, collapse = ','), ")
group by 1,2"))

pharmacy = sqlQuery(ch, paste("SELECT rs.z_patid, rs.fill_dt as fst_dt,
SUM(rs.calc_allwd) as cost, 'Pharmacy' as claim_type
from HCCI_2.RX_SDDV2 rs
where rs.yr > 2015 AND rs.z_patid in (",
                              paste(patients, collapse = ','), ")",
                              "group by 1,2
"))


pt0 = inp %>% 
  bind_rows(op) %>% 
  bind_rows(ph) %>% 
  bind_rows(pharmacy) %>%
  left_join(first_dates) %>% 
  mutate('treatment' = ifelse(fst_dt<first_date, 'Before', 'After'),
         day_diff = as.numeric(fst_dt-first_date)) 

vc = c()
for (val in seq(-5, 5, 1)){
  if(val!=0){
    if(val<0){
      vc = append(vc, paste(abs(val), 'years before'))
      if(val == -1){
        vc[length(vc)] = 'last year'
      }
    }else{
      vc = append(vc, paste(abs(val), 'years after'))
      if(val == 1){
        vc[length(vc)] = '1st year'
      }
    }
  }
}

pt0$time_cat = cut(pt0$day_diff, 
                   breaks = seq(-365*5, 365*5, by = 365), 
                   right = FALSE,
                   labels = vc, ordered_result = TRUE)

df_main = pt0 %>% 
  filter(claim_type!='Pharmacy') %>% 
  group_by(z_patid, time_cat,treatment) %>% 
  summarise(cost = sum(cost))
df_info = df_main %>% 
  group_by(time_cat) %>% 
  summarise(average = log(median(cost+1),2),
            patients = n_distinct(z_patid),
            qu = max(log((cost+1),2))*1.1) %>% 
  mutate(qu = max(qu))

df_main %>% ggplot(mapping = aes(x = time_cat, 
                                 y = log(cost+1,2)))+
  geom_boxplot(aes(fill = treatment))+
  scale_fill_manual(values = cbp)+
  geom_point(data = df_info, 
             mapping = aes(x = time_cat, y = average, size = patients),
             color = 'black')+
  geom_line(data = df_info, 
            mapping = aes(x = time_cat, y = average, group = 1),
            color = 'black', size = 1)+
  geom_label(data = df_info, mapping = aes(x = time_cat, 
                                           y = qu, label = patients),
             fill = 'green')+
  ylab('Cost log()')+
  xlab('')+
  ggtitle("Costs (without pharmacy) in different time intervals",
          subtitle = '2016-2020')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12))

ggsave(paste(report_directory, 'time_intervals_cost_Cystic_Fibrosis.png', sep = ""),
       dpi = 600, width = 7.2, height = 5, units = 'in')  


pt = pt0%>% 
  group_by(z_patid, treatment, claim_type) %>% 
  summarise(cost = sum(cost),
            claim_count = n())

pt$treatment = factor(pt$treatment,
                      levels = c('Before', 'After'))

p1 = ggplot(pt,
            aes(x=treatment, y = log(cost+1), fill = treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape = 1, aes(color = treatment), 
             alpha = .15, position = position_jitter())+
  facet_wrap(~claim_type, nrow = 1)+
  ylab(label = 'Cost log()')+
  xlab(label = "")+
  theme_bw()+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b=-5, unit = 'in'))+
  scale_fill_manual(values = cbp[1:2])+
  scale_color_manual(values = cbp[1:2])

p2 = ggplot(pt,
            aes(x=treatment, y = log(claim_count), fill = treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape = 1, aes(color = treatment), 
             alpha = .15, position = position_jitter())+
  facet_wrap(~claim_type, nrow = 1)+
  ylab(label = 'Claim counts')+
  xlab(label = "")+
  theme_bw()+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(t = -5, unit = 'in'))+
  scale_fill_manual(values = cbp[1:2])+
  scale_color_manual(values = cbp[1:2])

p3 = cowplot::plot_grid(p1,
                        p2,nrow = 2,
                        align = 'hv',
                        labels = c('a', 'b')
)



write.csv(pt0, paste(report_directory,'Cystic_Fibrosis.csv', sep = ''), row.names = FALSE)
ggsave(plot = p3, filename = paste(report_directory, 'Cystic_Fibrosis.png', sep = ""),
       dpi = 600, width = 7.2, height = 5.5, units = 'in')

pvals = c()
mean_vec = c()
median_vec = c()
for (i in c('all','Pharmacy','Without Pharmacy')){
  df = pt0 %>% 
    group_by(z_patid) %>% 
    summarise(start_day = min(fst_dt),
              last_day = max(fst_dt),
              first_date = max(first_date)) %>% 
    mutate(days_before = as.numeric(first_date - start_day),
           days_after = as.numeric(last_day -first_date))
  if (i == 'all'){
    df2 = pt0 %>% 
      group_by(z_patid, treatment) %>% 
      summarise(cost = n()) %>% 
      spread(treatment, cost)
  }else{
    if(i == 'Pharmacy'){
      df2 = pt0 %>% 
        filter(claim_type=='Pharmacy') %>%
        group_by(z_patid, treatment) %>% 
        summarise(cost = n()) %>% 
        spread(treatment, cost)
    }else{
      df2 = pt0 %>% 
        filter(claim_type!='Pharmacy') %>%
        group_by(z_patid, treatment) %>% 
        summarise(cost = n()) %>% 
        spread(treatment, cost)
    }
  }
  
  df3 = df2 %>% 
    drop_na() %>% 
    left_join(df) %>% 
    filter(days_before>90, days_after >90)
  
  df3 = df3 %>% 
    mutate(daily_before = Before/days_before,
           daily_after = After/days_after,
           cost_diff = daily_after-daily_before)
  
  df3$claim = i
  wl = wilcox.test(df3$cost_diff)
  pvals = append(pvals,wl$p.value)
  mean_vec =append(mean_vec,mean(df3$cost_diff))
  median_vec = append(median_vec,median(df3$cost_diff))
  
  if(i =='all'){
    tmp = df3  
  }else{
    tmp = tmp %>% bind_rows(df3)
  }
  
}
tmp$cost_diff = tmp$cost_diff*365
mu = tmp %>% group_by(claim) %>% summarise(mn = mean(cost_diff))

p1 = ggplot(tmp %>% filter(between(cost_diff, -150,150)) ,
            aes(x = cost_diff, color = claim, fill = claim))+
  geom_histogram(aes(y=..density..), position = 'identity',
                 alpha = 0.2)+
  geom_density(alpha = .3)+
  geom_vline(data = mu, aes(xintercept = mn, color = claim), linetype = 'dashed')+
  scale_color_manual(values = cbp[1:3])+
  scale_fill_manual(values = cbp[1:3])+
  ggtitle('Difference in the number of claims')+
  ylab(label = 'Density')+
  xlab(label = 'Claims frequency difference (per year)')+
  theme_classic()+
  theme(legend.position = c(.2,.6))



ggsave(p1, filename = paste(report_directory, 'Cystic_Fibrosis_claims_plot_year.png', sep = ""),
       dpi = 600, width = 7.2, height = 5.5, units = 'in')  


## costs
for (i in c('all','Pharmacy','Without Pharmacy')){
  df = pt0 %>% 
    group_by(z_patid) %>% 
    summarise(start_day = min(fst_dt),
              last_day = max(fst_dt),
              first_date = max(first_date)) %>% 
    mutate(days_before = as.numeric(first_date - start_day),
           days_after = as.numeric(last_day -first_date))
  if (i == 'all'){
    df2 = pt0 %>% 
      group_by(z_patid, treatment) %>% 
      summarise(cost = sum(cost)) %>% 
      spread(treatment, cost)
  }else{
    if(i == 'Pharmacy'){
      df2 = pt0 %>% 
        filter(claim_type=='Pharmacy') %>%
        group_by(z_patid, treatment) %>% 
        summarise(cost = sum(cost)) %>% 
        spread(treatment, cost)
    }else{
      df2 = pt0 %>% 
        filter(claim_type!='Pharmacy') %>%
        group_by(z_patid, treatment) %>% 
        summarise(cost = sum(cost)) %>% 
        spread(treatment, cost)
    }
  }
  
  df3 = df2 %>% 
    drop_na() %>% 
    left_join(df) %>% 
    filter(days_before>90, days_after >90)
  
  df3 = df3 %>% 
    mutate(daily_before = Before/days_before,
           daily_after = After/days_after,
           cost_diff = daily_after-daily_before)
  
  df3$claim = i
  wl = wilcox.test(df3$cost_diff)
  pvals = append(pvals,wl$p.value)
  mean_vec =append(mean_vec,mean(df3$cost_diff))
  median_vec = append(median_vec,median(df3$cost_diff))
  
  if(i =='all'){
    tmp = df3  
  }else{
    tmp = tmp %>% bind_rows(df3)
  }
  
}

tmp$cost_diff = tmp$cost_diff*365/1000
mu = tmp %>% group_by(claim) %>% summarise(mn = mean(cost_diff))

p2 = ggplot(tmp %>% filter(between(cost_diff, -500, 800)),
            aes(x = cost_diff, color = claim, fill = claim))+
  geom_histogram(aes(y=..density..), position = 'identity',
                 alpha = 0.2)+
  geom_density(alpha = .3)+
  geom_vline(data = mu, aes(xintercept = mn, color = claim), linetype = 'dashed')+
  scale_color_manual(values = cbp[1:3])+
  scale_fill_manual(values = cbp[1:3])+
  ggtitle('Difference in the costs')+
  ylab(label = 'Density')+
  xlab(label = 'Costs difference (thousands of dollars per year)')+
  theme_classic()



ggsave(p2, filename = paste(report_directory, 'Cystic_Fibrosis_costs_plot_year.png', sep = ""),
       dpi = 600, width = 7.2, height = 5.5, units = 'in') 

p3 = cowplot::plot_grid(p1,byrow = F,
                        p2+theme(legend.position = 'none' ),nrow = 1,
                        labels = c('a', 'b')
)
ggsave(p3, filename = paste(report_directory, 'Cystic_Fibrosis_comb_plot_year.png', sep = ""),
       dpi = 600, width = 10, height = 3, units = 'in') 


p3_cf_dist =  cowplot::plot_grid(p1+theme(axis.title.x = element_blank()),
                                 p2+theme(legend.position = 'none',
                                          axis.title.x = element_blank(),
                                          axis.title.y = element_blank()),
                                 byrow = F,nrow = 1,
                                 labels = c('a', 'b')
                                 )
