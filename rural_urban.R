library(tidyverse)
library(RODBC)
library(scales)
library(usmap)
library(readxl)

ch = odbcConnect('vertica',
                 uid = 'XXXXX',
                 pwd = 'XXX')
cbp <- c('#e69f00', '#0072b2', '#009e73', 
         '#d55e00', '#f0e442', '#56b4e9')

rural_urban = read_excel('RUCA2010zipcode_1.xlsx', sheet = 'Data')
rural_urban_def = read_excel('RUCA2010zipcode_1.xlsx', sheet = 'description')
rural_urban = rural_urban %>% left_join(rural_urban_def)

key = read_excel('keys.xlsx', sheet = 'disease')

patient_list = read.csv('reports/prevalence/overall/patient_list.csv')


zip_members = sqlQuery(ch, paste("SELECT mbr_zip_5_cd, 
count(distinct z_patid) as members
from HCCI_2.MBR_SDDV2
where yr >'2015'
group by 1;"), stringsAsFactors = F,  as.is = c(T,T))
zip_members$members = as.numeric(zip_members$members)

zip_members$mbr_zip_5_cd = zip_members$mbr_zip_5_cd %>% sprintf(fmt = '%05d')
patient_list$mbr_zip_5_cd = patient_list$mbr_zip_5_cd %>% sprintf(fmt = '%05d')

zip_members = zip_members %>% 
  left_join(rural_urban, by = c('mbr_zip_5_cd' = 'ZIP_CODE'))
zip_members = zip_members %>% drop_na()

diseases = unique(patient_list$disease)

for (i in 1:length(diseases)) {
  if(i==1){
    prev_table = zip_members %>% mutate(disease = diseases[i])
  }else{
    prev_table = prev_table %>% bind_rows(
      zip_members %>% mutate(disease = diseases[i])
    )
  }
  
}

df_patient = patient_list %>% 
  group_by(disease, mbr_zip_5_cd) %>% 
  summarise(patients= n())

prev_table = prev_table %>% 
  left_join(df_patient)
prev_table$prevalence = prev_table$patients*100000/prev_table$members
prev_table$prevalence[is.na(prev_table$prevalence)] = 0
prev_table$patients[is.na(prev_table$patients)] = 0
prev_table = prev_table %>% left_join(key)



# do the same plot without dividing 
# us heat map based on the zip codes
prev_table$short_description = factor(prev_table$short_description, 
                                   levels = c('Metropolitan', 'Micropolitan',
                                              'Small town', 'Rural area'))

ggplot(prev_table %>% filter(short_description !='Unknown'),
       aes(x=short_description, y = log(prevalence+1), fill = short_description))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape = 1, aes(color = short_description), 
             alpha = .15, position = position_jitter())+
  facet_wrap(~short_disease, nrow = 2, ncol = 7)+
  ylab(label = 'Prevalence - 100,000 per members log()')+
  xlab(label = "")+
  theme_bw()+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = 'bold'),
        axis.text.x = element_text(angle = 90, vjust = 0.5,  hjust = 1))+
  scale_fill_manual(values = cbp[1:4])+
  scale_color_manual(values = cbp[1:4])


ggsave('./reports/prevalence/overall/rural_urban_all_areas.png',
       dpi = 600, width = 7.2, height = 5.5, units = 'in')

for (i in 1:length(diseases)) {
  
  tmp = prev_table %>% filter(disease==diseases[i])
  kt = kruskal.test(tmp$prevalence~tmp$short_description)
  wl = pairwise.wilcox.test(tmp$prevalence,tmp$short_description)
  
  tmp = tmp %>% 
    group_by(short_description) %>% 
    summarise(med = mean(prevalence)) %>% 
    drop_na() %>% 
    spread(short_description, med) %>% 
    mutate(p_value = kt$p.value,
           disease = diseases[i])
  tmp2 = wl$p.value %>% 
    data.frame() %>% 
    mutate(gr1 = rownames(.)) %>% 
    gather('gr2', 'pvalue', -'gr1') %>% 
    drop_na() %>% 
    mutate(diseas = diseases[i])
  
  if(i==1){
    ks_test = tmp
    wl_test = tmp2
  }else{
    ks_test = ks_test %>% 
      bind_rows(tmp)
    wl_test = wl_test %>% 
      bind_rows(tmp2)
  }
  
}

write.csv(ks_test, './reports/prevalence/overall/kruskall_all.csv', row.names = FALSE)
write.csv(wl_test, './reports/prevalence/overall/wil_pairwise_all.csv', row.names = FALSE)


ggplot(prev_table %>% 
         filter(short_description !='Unknown') %>% 
         filter(prevalence>0),
       aes(x=short_description, y = log(prevalence+1), fill = short_description))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape = 1, aes(color = short_description), 
             alpha = .15, position = position_jitter())+
  facet_wrap(~short_disease, nrow = 2, ncol = 7)+
  ylab(label = 'Prevalence - 100,000 per members log()')+
  xlab(label = "")+
  theme_bw()+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = 'bold'),
        axis.text.x = element_text(angle = 90, vjust = 0.5,  hjust = 1))+
  scale_fill_manual(values = cbp[1:4])+
  scale_color_manual(values = cbp[1:4])

ggsave('./reports/prevalence/overall/rural_urban_only_with_patient.png',
       dpi = 600, width = 7.2, height = 5.5, units = 'in')

for (i in 1:length(diseases)) {
  
  tmp = prev_table %>%
    filter(disease==diseases[i]) %>% 
    filter(prevalence>0)
  kt = kruskal.test(tmp$prevalence~tmp$short_description)
  tmp = tmp %>% 
    group_by(short_description) %>% 
    summarise(med = median(prevalence)) %>% 
    drop_na() %>% 
    spread(short_description, med) %>% 
    mutate(p_value = kt$p.value,
           disease = diseases[i])
  
  tmp2 = wl$p.value %>% 
    data.frame() %>% 
    mutate(gr1 = rownames(.)) %>% 
    gather('gr2', 'pvalue', -'gr1') %>% 
    drop_na() %>% 
    mutate(diseas = diseases[i])
  
  if(i==1){
    ks_test = tmp
    wl_test = tmp2
  }else{
    ks_test = ks_test %>% 
      bind_rows(tmp)
    wl_test = wl_test %>% 
      bind_rows(tmp2)
  }
  
}

write.csv(ks_test, './reports/prevalence/overall/kruskall_only_with_patients.csv', row.names = FALSE)
write.csv(wl_test, './reports/prevalence/overall/wil_pairwise_only_with_patients.csv', row.names = FALSE)
