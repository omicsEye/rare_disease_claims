library(tidyverse)
library(RODBC)
library(scales)
library(usmap)


plot_list = list()
sex_test = data.frame(disease=character(), sex = character(), patients=integer(),
                      total=integer(), non_patient=integer(), chisq_pvalue=double(), 
                      fisher_pvalue=double(), or= double(), lb=double(), ub=double())

age_test = data.frame(disease=character(), age = integer(), patients=integer(), 
                      total=integer(), non_patient=integer(), chisq_pvalue=double())

age_sex_test = data.frame(disease=character(), age = integer(), sex = character(),
                          patients=integer(), total=integer(), 
                          non_patient=integer(), chisq_pvalue=double(), 
                          fisher_pvalue=double(), or= double(), lb=double(), ub=double())

patient_list = data.frame(disease = character(), z_patid = integer(), sex = character(),
                          age_band_cd = integer(), mbr_zip_5_cd = integer(), mbr_state = character())

ch = odbcConnect('vertica',
                 uid = 'XXXXX',
                 pwd = 'XXX')
icd10  = list('Cystic_Fibrosis' = c("'E840'", "'E8411'", "'E8419'", "'E848'", "'E849'"), 
            'Sickle_Cell'= c("'D5700'", "'D5701'", "'D5702'", "'D571'", "'D5720'",
                             "'D57211'", "'D57212'", "'D57212'", "'D57219'",
                             "'D5740'", "'D57411'", "'D57412'", "'D57419'",
                             "'D5780'", "'D57811'", "'D57812'", "'D57819'"),
            'Muscular_Dystrophy' = c("'G710'", "'G7101'", "'G7100'", "'G7102'",
                                       "'G7109'", "'G7111'"),
            'Lennox_Gastaut_Syndrome' = c("'G40811'", "'G40813'", "'G40812'",
                                          "'G40814'"),
            'Urea_Cycle_Disorder' = c("'E7220'","'E7221'","'E7222'","'E7223'",
                                      "'E7229'", "'E724'"),
            "Takayasu's_Arteritis" = c("'M314'"),
            'Pheochromocytoma' = c("'C7410'", "'D3500'"),
            'Hereditary_Hemorrhagic_Telangiectasia' = c("'I780'"),
            'Osteogenesis_Imperfecta' = c("'Q780'"),
            'Eosinophilic_Esophagitis' = c("'K200'"),
            'Charcot_Marie_Tooth' = c("'G600'"),
            'Batten_Disease' = c("'E754'", "'G111'"),
            'Focal_and_Segmental_Glomerulosclerosis' = c("'N040'","'N041'", "'N042'",
                                                         "'N047'","'N031'"),
            'Mitochondrial_Neurogastrointestinal_Encephalopathy' = c("'E8849'")
            )
state_members = sqlQuery(ch, paste("SELECT mbr_state, 
count(distinct z_patid) as members
from HCCI_2.MBR_SDDV2
where yr >'2015'
group by 1;"))

write.csv(state_members %>% drop_na(), 'state_members.csv', row.names = FALSE)

state_members = read.csv('state_members.csv')
state_members = state_members %>% drop_na()

state_members_sex_age = read.csv('state_members_sex_age.csv')
state_members_sex_age = state_members_sex_age %>% drop_na()

state_members_yr = sqlQuery(ch, paste("SELECT mbr_state,yr, 
count(distinct z_patid) as members
from HCCI_2.MBR_SDDV2
where yr >'2015'
group by 1,2;"))
write.csv(state_members_yr %>% drop_na(), 'state_members_yr.csv', row.names = FALSE)

state_members_yr = read.csv('state_members_yr.csv')
state_members_yr = state_members_yr %>% drop_na()

diseases = names(icd10)
t = Sys.time()
for(i in 1:length(diseases)){
print(i)
print(diseases[i])
report_directory = paste(getwd(), "/", 'reports/prevalence/overall/', diseases[i], sep = "")
ifelse(!dir.exists(report_directory), dir.create(report_directory),FALSE)

condition = paste(
  paste('inp.icd10_cm',1:10, ' in (', 
        paste(icd10[[diseases[i]]], collapse = ","),
        ') or',sep = ''),collapse = ' ')
condition = substr(condition, 1, nchar(condition)-3)

inp = sqlQuery(ch, paste("SELECT inp.z_patid, inp.fst_dt,
SUM(inp.calc_allwd) as cost
from HCCI_2.INP_SDDV2 inp
where",
condition,
"
group by 1,2
having SUM(inp.calc_allwd) > 0;"))


condition = paste(
  paste('op.diag_icd10_cm',1:10, ' in (', 
        paste(icd10[[diseases[i]]], collapse = ","),
        ') or',sep = ''),collapse = ' ')
condition = substr(condition, 1, nchar(condition)-3)


op = sqlQuery(ch, paste("SELECT op.z_patid, op.fst_dt,
SUM(op.calc_allwd) as cost
from HCCI_2.OP_SDDV2 op
where",
condition,
"group by 1,2
having SUM(op.calc_allwd) > 0;"))

condition = paste(
  paste('ph.diag_icd10_cm',1:10, ' in (', 
        paste(icd10[[diseases[i]]], collapse = ","),
        ') or',sep = ''),collapse = ' ')
condition = substr(condition, 1, nchar(condition)-3)

ph = sqlQuery(ch, paste("SELECT ph.z_patid, ph.fst_dt,
SUM(ph.calc_allwd) as cost
from HCCI_2.PHYS_SDDV2 ph 
where",
condition,
"group by 1,2
having SUM(ph.calc_allwd) > 0;"))

pt = inp %>% 
  bind_rows(op) %>% 
  bind_rows(ph) %>% 
  group_by(z_patid,fst_dt) %>% 
  summarise(cost = sum(cost)) %>% 
  group_by(z_patid) %>% 
  mutate(rank = row_number(fst_dt),
         behind = lag(fst_dt),
         time_diff= as.numeric(fst_dt-behind)) %>% 
  drop_na() %>% 
  filter(time_diff < 91) %>% 
  group_by(z_patid) %>% 
  summarise(claims = n()+1)


patients_query = paste("'", pt$z_patid, "'", sep = "")

members = sqlQuery(ch, paste("select ss.z_patid, ss.sex, ss.age_band_cd, ss.mbr_zip_5_cd, ss.mbr_state
from(
SELECT *, 
ROW_NUMBER() over(PARTITION by z_patid ORDER by yr DESC , mnth DESC) as rn
FROM HCCI_2.MBR_SDDV2 ms 
where ms.z_patid in (", paste(patients_query, collapse = ','),
")) as ss
where ss.rn = 1;"))

 
patient_temp = members %>% mutate(disease = diseases[i])
patient_list = patient_list %>% bind_rows(patient_temp)

df_main = inp %>% 
  bind_rows(op) %>% 
  bind_rows(ph) %>% 
  mutate(yr = lubridate::year(fst_dt)) %>% 
  group_by(z_patid, yr) %>% 
  summarise(cost = sum(cost),
            claim_count = n()) %>% 
  filter(z_patid %in% pt$z_patid) %>% 
  full_join(members)


## tests
tmp_t = df_main %>% 
  group_by(sex) %>% 
  summarise(patients = n_distinct(z_patid)) %>% 
  left_join(state_members_sex_age %>% 
              group_by(sex) %>% 
              summarise(total = sum(members))) %>% 
  mutate(non_patient = total-patients) %>% 
  drop_na()

M = as.matrix(tmp_t[,c(F,T,F,T)])
ct = chisq.test(t(M))$p.value
ft = fisher.test(M)

sex_test_tmp = data.frame(disease=character(), sex = character(), patients=integer(), total=integer(), 
                          non_patient=integer(), chisq_pvalue=double(), 
                          fisher_pvalue=double(), or= double(), lb=double(), ub=double())
sex_test_tmp[1:nrow(tmp_t),1] =diseases[i]
sex_test_tmp[,c('sex',"patients", "total", "non_patient")] = tmp_t
sex_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value
sex_test_tmp$fisher_pvalue = ft$p.value
sex_test_tmp$or = ft$estimate
sex_test_tmp[,'lb'] = ft$conf.int[1]
sex_test_tmp[,'ub'] = ft$conf.int[2]

sex_test = sex_test %>% bind_rows(sex_test_tmp)


tmp_t = df_main %>% 
  group_by(age_band_cd) %>% 
  summarise(patients = n_distinct(z_patid)) %>% 
  left_join(state_members_sex_age %>% 
              group_by(age_band_cd) %>% 
              summarise(total = sum(members))) %>% 
  mutate(non_patient = total-patients) %>% 
  drop_na()


M = as.matrix(tmp_t[,c(F,T,F,T)])

age_test_tmp = data.frame(disease=character(), age = integer(), patients=integer(), total=integer(), 
                          non_patient=integer(), chisq_pvalue=double())
age_test_tmp[1:nrow(tmp_t),1] =diseases[i]
age_test_tmp[,c('age',"patients", "total", "non_patient")] = tmp_t
age_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value

age_test = age_test %>% bind_rows(age_test_tmp)

for(j in 1:7){
  tmp_t = df_main %>% 
    filter(age_band_cd==j) %>% 
    group_by(sex) %>% 
    summarise(patients = n_distinct(z_patid)) %>% 
    left_join(state_members_sex_age %>% 
                filter(age_band_cd==j) %>% 
                group_by(sex) %>% 
                summarise(total = sum(members))) %>% 
    mutate(non_patient = total-patients) %>% 
    drop_na()
  if(nrow(tmp_t)==2){
    M = as.matrix(tmp_t[,c(F,T,F,T)])
    ct = chisq.test(t(M))$p.value
    ft = fisher.test(M)
    ft$p.value
    c(ft$conf.int)
    as.numeric(ft$estimate)
    
    age_sex_test_tmp = data.frame(disease=character(),age=integer(), sex = character(),
                                  patients=integer(), total=integer(), 
                                  non_patient=integer(), chisq_pvalue=double(), 
                                  fisher_pvalue=double(), or= double(), lb=double(), ub=double())
    age_sex_test_tmp[1:nrow(tmp_t),1] =diseases[i]
    age_sex_test_tmp[1:nrow(tmp_t),'age'] =j
    age_sex_test_tmp[,c('sex',"patients", "total", "non_patient")] = tmp_t
    age_sex_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value
    age_sex_test_tmp$fisher_pvalue = ft$p.value
    age_sex_test_tmp$or = ft$estimate
    age_sex_test_tmp[,'lb'] = ft$conf.int[1]
    age_sex_test_tmp[,'ub'] = ft$conf.int[2]
    
    age_sex_test = age_sex_test %>% bind_rows(age_sex_test_tmp)
  }
}
## overall prevalence
tmp = df_main %>% 
  filter(yr>2015) %>% 
  group_by(mbr_state) %>% 
  summarise(count = n_distinct(z_patid))

data("statepop")
statepop = statepop %>% 
  left_join(tmp, by = c('abbr' = 'mbr_state')) %>% 
  replace(is.na(.), 0) %>% 
  left_join(state_members, by = c('abbr' = 'mbr_state')) %>% 
  mutate(count = replace(count, count>0 & count<10 , NA),
    patient_per_hundered_thousands= count*100000/pop_2015,
         patient_per_hundered_thousands_mem= count*100000/members)


p1 = plot_usmap(data = statepop,size = 0.1 ,
             values = 'patient_per_hundered_thousands') +
    labs(title = paste("Prevalence of ",  gsub("_", " ", diseases[i]) ),
         subtitle = '2016-2020')+
    scale_fill_continuous(low = 'white', high = 'red',
                          name = gsub('_', ' ', 'Patients per 100,000 people\n'), 
                          label = scales::comma)+
    theme(legend.position = c(0.9, 0.2),
          text = element_text(family = 'sans'),
          plot.title = element_text(face = "bold",size = 10),
          plot.subtitle = element_text(size = 8),
          legend.title = element_text(face="bold.italic", size = 6),
          legend.text = element_text(size = 4),
          legend.key.size = unit(.15, 'in'),
          legend.background = element_blank())

ggsave(plot = p1, paste(report_directory,'/',diseases[i], '.png', sep = ""),
       dpi = 600, width = 7.2, height = 3, units = 'in')
plot_list[[paste(diseases[i],'_total_prev_pop', sep = '')]] = p1

p2 = plot_usmap(data = statepop,size = 0.1 ,
                values = 'patient_per_hundered_thousands_mem') +
  labs(title = paste("Prevalence of ",  gsub("_", " ", diseases[i]) ),
       subtitle = '2016-2020')+
  scale_fill_continuous(low = 'white', high = 'red',
                        name = gsub('_', ' ', 'Patients per 100,000 members\n'), 
                        label = scales::comma)+
  theme(legend.position = c(0.9, 0.2),
        text = element_text(family = 'sans'),
        plot.title = element_text(face = "bold",size = 10),
        plot.subtitle = element_text(size = 8),
        legend.title = element_text(face="bold.italic", size = 6),
        legend.text = element_text(size = 4),
        legend.key.size = unit(.15, 'in'),
        legend.background = element_blank())

ggsave(plot = p2,paste(report_directory,'/',diseases[i], '_members.png', sep = ""),
       dpi = 600, width = 7.2, height = 3, units = 'in')
plot_list[[paste(diseases[i],'_total_prev_mem', sep = '')]] = p2

write.csv(statepop, paste(report_directory,'/',diseases[i], '.csv', sep = ""))

tmp = df_main %>% 
  group_by(yr, mbr_state) %>% 
  summarise(count = n_distinct(z_patid)) %>% 
  filter(yr>2015) %>% 
  spread(yr, count) %>% 
  replace(is.na(.),0) %>% 
  gather('yr', 'count', -mbr_state)

data("statepop")
statepop = statepop %>% 
  left_join(tmp,by = c('abbr' = 'mbr_state'))

if(sum(is.na(statepop$count))>0){
  df_na = statepop[is.na(statepop$count),]
  for(r in 1:nrow(df_na)){
    for(year in (2016:2020)){
    statepop[nrow(statepop)+1,] = c(df_na[r,1:4], as.character(year), 0)
    }
  }
}
statepop$yr = as.numeric(statepop$yr)
statepop = statepop %>% 
  drop_na() %>% 
  left_join(state_members_yr, by = c('abbr' = 'mbr_state', 'yr'='yr')) %>% 
  mutate(count = replace(count, count>0 & count<10 , NA),
         patient_per_hundered_thousands = count*100000/pop_2015,
         patient_per_hundered_thousands_mem = count*100000/members)

p3 = plot_usmap(data = statepop, size = 0.01,
           values = 'patient_per_hundered_thousands') +
  labs(title = paste("Prevalence of ",  gsub("_", " ", diseases[i])),
       subtitle = '2016-2020')+
  scale_fill_continuous(low = 'white', high = 'red',
                        name = gsub('_', ' ', 'Patients per 100,000 people\n'), 
                        label = scales::comma)+
  facet_wrap(~yr,nrow = 1)+
  theme(legend.position ='right',
        text = element_text(family = 'sans'),
        plot.title = element_text(face = "bold",size = 10),
        plot.subtitle = element_text(size = 8),
        legend.title = element_text(face="bold.italic", size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(.1, 'in'),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6))
ggsave(plot = p3, paste(report_directory,'/',diseases[i], '_years.png', sep = ""),
       dpi = 600, width = 7.2, height = 3, units = 'in')
plot_list[[paste(diseases[i],'_year_prev', sep = '')]] = p3

p4 = plot_usmap(data = statepop, size = 0.01,
                values = 'patient_per_hundered_thousands_mem') +
  labs(title = paste("Prevalence of ",  gsub("_", " ", diseases[i])),
       subtitle = '2016-2020')+
  scale_fill_continuous(low = 'white', high = 'red',
                        name = gsub('_', ' ', 'Patients per 100,000 members\n'), 
                        label = scales::comma)+
  facet_wrap(~yr,nrow = 1)+
  theme(legend.position ='right',
        text = element_text(family = 'sans'),
        plot.title = element_text(face = "bold",size = 10),
        plot.subtitle = element_text(size = 8),
        legend.title = element_text(face="bold.italic", size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(.1, 'in'),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6))
ggsave(plot = p4, paste(report_directory,'/',diseases[i], '_years_members.png', sep = ""),
       dpi = 600, width = 7.2, height = 3, units = 'in')
plot_list[[paste(diseases[i],'_year_prev_mem', sep = '')]] = p4

write.csv(statepop, paste(report_directory,'/',diseases[i], '_years.csv', sep = ""))

}
Sys.time()-t
report_directory = paste(getwd(), "/", 'reports/prevalence/overall/', sep = "")
write.csv(age_test, paste(report_directory,'age_test.csv', sep = ''), row.names = FALSE)
write.csv(sex_test, paste(report_directory,'sex_test.csv', sep = ''), row.names = FALSE)
write.csv(age_sex_test, paste(report_directory,'age_sex_test.csv', sep = ''), row.names = FALSE)
write.csv(patient_list, paste(report_directory,'patient_list.csv', sep = ''), row.names = FALSE)


report_directory = paste(getwd(), "/", 'reports/prevalence/overall/', sep = "")
saveRDS(object = plot_list, file = paste(report_directory, 'plot_list.RData',sep = ''))
plot_list = readRDS(file = 'reports/prevalence/overall/plot_list.RData')

p1 = plot_list$Cystic_Fibrosis_year_prev_mem+
  ggtitle(label = 'a: Cystic Fibrosis',subtitle = '')+
  theme(plot.margin = margin(0,0,0,0))

p2 = plot_list$Sickle_Cell_year_prev_mem+
  ggtitle(label = 'b: Sickle Cell',subtitle = element_blank())+
  theme(strip.text = element_blank(),
        plot.title = element_text(margin = margin(0,0,-1,0,unit = 'cm')))

cowplot::plot_grid(p1,
                   p2,nrow = 2,
                   align = 'hv'
                   )
ggsave(paste(report_directory,'main_figure.png', sep = ""),
       dpi = 600, width = 7.2, height = 3, units = 'in')
