library(tidyverse)
library(RODBC)
library(scales)
library(usmap)


plot_list = list()
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
state_members = sqlQuery(ch, paste("SELECT mbr_state, sex, age_band_cd,
count(distinct z_patid) as members
from HCCI_2.MBR_SDDV2
where yr >'2015'
group by 1,2,3;"))
write.csv(state_members %>% drop_na(), 'state_members_sex_age.csv', row.names = FALSE)
state_members = read.csv('state_members_sex_age.csv')
state_members = state_members %>% drop_na()

diseases = names(icd10)
sex_test = data.frame(disease=character(), sex = character(), patients=integer(),
                      total=integer(), non_patient=integer(), chisq_pvalue=double(), 
                      fisher_pvalue=double(), or= double(), lb=double(), ub=double())

age_test = data.frame(disease=character(), age = integer(), patients=integer(), 
                      total=integer(), non_patient=integer(), chisq_pvalue=double())

age_sex_test = data.frame(disease=character(), age = integer(), sex = character(),
                          patients=integer(), total=integer(), 
                      non_patient=integer(), chisq_pvalue=double(), 
                      fisher_pvalue=double(), or= double(), lb=double(), ub=double())


for(i in 1:length(diseases)){
  
  report_directory = paste(getwd(), "/", 'reports/prevalence/overall/', diseases[i], sep = "")
  ifelse(!dir.exists(report_directory), dir.create(report_directory),FALSE)
  
  inp = sqlQuery(ch, paste("SELECT inp.z_patid, inp.yr,
SUM(inp.calc_allwd) as inp_cost
from HCCI_2.INP_SDDV2 inp
where 
inp.icd10_cm1 in (", paste(icd10[[diseases[i]]], collapse = ",") ,
                           ") group by 1,2
having SUM(inp.calc_allwd) > 0;"))
  
  
  op = sqlQuery(ch, paste("SELECT op.z_patid, op.yr,
SUM(op.calc_allwd) as op_cost
from HCCI_2.OP_SDDV2 op
where 
op.diag_icd10_cm1 in (", paste(icd10[[diseases[i]]], collapse = ",") ,
                          ") group by 1,2
having SUM(op.calc_allwd) > 0;"))
  
  ph = sqlQuery(ch, paste("SELECT ph.z_patid, ph.yr,
SUM(ph.calc_allwd) as ph_cost
from HCCI_2.PHYS_SDDV2 ph 
where ph.diag_icd10_cm1 in (", paste(icd10[[diseases[i]]], collapse = ",") ,
                          ") group by 1,2
having SUM(ph.calc_allwd) > 0;"))
  
patients = unique(c(inp$z_patid, op$z_patid, ph$z_patid))
patients_query = paste("'", patients, "'", sep = "")
  
members = sqlQuery(ch, paste("select ss.z_patid, ss.sex, ss.age_band_cd, ss.mbr_zip_5_cd, ss.mbr_state
from(
SELECT *, 
ROW_NUMBER() over(PARTITION by z_patid ORDER by yr DESC , mnth DESC) as rn
FROM HCCI_2.MBR_SDDV2 ms 
where ms.z_patid in (", paste(patients_query, collapse = ','),
                               ")) as ss
where ss.rn = 1;"))
  
df_main = inp %>% 
    full_join(op) %>% 
    full_join(ph) %>% 
    full_join(members)

tmp = df_main %>% 
  group_by(sex) %>% 
  summarise(patients = n_distinct(z_patid)) %>% 
  left_join(state_members %>% 
              group_by(sex) %>% 
              summarise(total = sum(members))) %>% 
  mutate(non_patient = total-patients) %>% 
  drop_na()

M = as.matrix(tmp[,c(F,T,F,T)])
ct = chisq.test(t(M))$p.value
ft = fisher.test(M)
ft$p.value
c(ft$conf.int)
as.numeric(ft$estimate)

sex_test_tmp = data.frame(disease=character(), sex = character(), patients=integer(), total=integer(), 
                      non_patient=integer(), chisq_pvalue=double(), 
                      fisher_pvalue=double(), or= double(), lb=double(), ub=double())
sex_test_tmp[1:nrow(tmp),1] =diseases[i]
sex_test_tmp[,c('sex',"patients", "total", "non_patient")] = tmp
sex_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value
sex_test_tmp$fisher_pvalue = ft$p.value
sex_test_tmp$or = ft$estimate
sex_test_tmp[,'lb'] = ft$conf.int[1]
sex_test_tmp[,'ub'] = ft$conf.int[2]

sex_test = sex_test %>% bind_rows(sex_test_tmp)

tmp = df_main %>% 
  group_by(age_band_cd) %>% 
  summarise(patients = n_distinct(z_patid)) %>% 
  left_join(state_members %>% 
              group_by(age_band_cd) %>% 
              summarise(total = sum(members))) %>% 
  mutate(non_patient = total-patients) %>% 
  drop_na()


M = as.matrix(tmp[,c(F,T,F,T)])

age_test_tmp = data.frame(disease=character(), age = integer(), patients=integer(), total=integer(), 
                      non_patient=integer(), chisq_pvalue=double())
age_test_tmp[1:nrow(tmp),1] =diseases[i]
age_test_tmp[,c('age',"patients", "total", "non_patient")] = tmp
age_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value

age_test = age_test %>% bind_rows(age_test_tmp)

for(j in 1:7){
tmp = df_main %>% 
  filter(age_band_cd==j) %>% 
  group_by(sex) %>% 
  summarise(patients = n_distinct(z_patid)) %>% 
  left_join(state_members %>% 
              filter(age_band_cd==j) %>% 
              group_by(sex) %>% 
              summarise(total = sum(members))) %>% 
  mutate(non_patient = total-patients) %>% 
  drop_na()
if(nrow(tmp)==2){
M = as.matrix(tmp[,c(F,T,F,T)])
ct = chisq.test(t(M))$p.value
ft = fisher.test(M)
ft$p.value
c(ft$conf.int)
as.numeric(ft$estimate)

age_sex_test_tmp = data.frame(disease=character(),age=integer(), sex = character(),
                              patients=integer(), total=integer(), 
                          non_patient=integer(), chisq_pvalue=double(), 
                          fisher_pvalue=double(), or= double(), lb=double(), ub=double())
age_sex_test_tmp[1:nrow(tmp),1] =diseases[i]
age_sex_test_tmp[1:nrow(tmp),'age'] =j
age_sex_test_tmp[,c('sex',"patients", "total", "non_patient")] = tmp
age_sex_test_tmp$chisq_pvalue = chisq.test(t(M))$p.value
age_sex_test_tmp$fisher_pvalue = ft$p.value
age_sex_test_tmp$or = ft$estimate
age_sex_test_tmp[,'lb'] = ft$conf.int[1]
age_sex_test_tmp[,'ub'] = ft$conf.int[2]

age_sex_test = age_sex_test %>% bind_rows(age_sex_test_tmp)
}
 
}
}

odbcClose(ch)

report_directory = paste(getwd(), "/", 'reports/prevalence/overall/', sep = "")
write.csv(age_test, paste(report_directory,'age_test.csv', sep = ''), row.names = FALSE)
write.csv(sex_test, paste(report_directory,'sex_test.csv', sep = ''), row.names = FALSE)
write.csv(age_sex_test, paste(report_directory,'age_sex_test.csv', sep = ''), row.names = FALSE)
