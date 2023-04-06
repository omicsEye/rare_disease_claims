library(tidyverse)
library(RODBC)
library(scales)
library(usmap)


patient_list = read.csv('reports/prevalence/overall/patient_list.csv')

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

diseases = names(icd10)
t = Sys.time()
for(i in 1:length(diseases)){
  print(i)
  print(diseases[i])
  # report_directory = paste(getwd(), "/", 'reports/cost_cluster/', diseases[i], sep = "")
  # ifelse(!dir.exists(report_directory), dir.create(report_directory),FALSE)

  patients = patient_list$z_patid[patient_list$disease==diseases[i]]
  patients = paste("'", patients, "'", sep = "")
  condition = paste(
    paste('inp.icd10_cm',1:10, ' in (', 
          paste(icd10[[diseases[i]]], collapse = ","),
          ') or',sep = ''),collapse = ' ')
  condition = paste("(", substr(condition, 1, nchar(condition)-3), ")")
  
  inp = sqlQuery(ch, paste("SELECT inp.z_patid, inp.fst_dt, 
SUM(inp.calc_allwd) as cost, 'Inpatient' as claim_type
from HCCI_2.INP_SDDV2 inp
where inp.yr > 2015 AND inp.z_patid in (", 
                           paste(patients, collapse = ','), ")
                           group by 1,2
having SUM(inp.calc_allwd) > 0;"))
  
  
  condition = paste(
    paste('op.diag_icd10_cm',1:10, ' in (', 
          paste(icd10[[diseases[i]]], collapse = ","),
          ') or',sep = ''),collapse = ' ')
  condition = paste("(", substr(condition, 1, nchar(condition)-3), ")")
  
  
  op = sqlQuery(ch, paste("SELECT op.z_patid, op.fst_dt,
SUM(op.calc_allwd) as cost, 'Outpatient' as claim_type
from HCCI_2.OP_SDDV2 op
where op.yr>2015 AND op.z_patid in (", 
paste(patients, collapse = ','), ")
group by 1,2
having SUM(op.calc_allwd) > 0;"))
  
  condition = paste(
    paste('ph.diag_icd10_cm',1:10, ' in (', 
          paste(icd10[[diseases[i]]], collapse = ","),
          ') or',sep = ''),collapse = ' ')
  condition =  paste("(", substr(condition, 1, nchar(condition)-3), ")")
  
  ph = sqlQuery(ch, paste("SELECT ph.z_patid, ph.fst_dt,
SUM(ph.calc_allwd) as cost, 'Physician' as claim_type
from HCCI_2.PHYS_SDDV2 ph 
where ph.yr > 2015 AND ph.z_patid in (", 
paste(patients, collapse = ','), ")
group by 1,2
having SUM(ph.calc_allwd) > 0;"))

  pharmacy = sqlQuery(ch, paste("SELECT rs.z_patid, rs.fill_dt as fst_dt,
SUM(rs.calc_allwd) as cost, 'Pharmacy' as claim_type
from HCCI_2.RX_SDDV2 rs
where rs.yr > 2015 AND rs.z_patid in (",
                          paste(patients, collapse = ','), ")",
                         "group by 1,2
having SUM(rs.calc_allwd) > 0;"))

  
  pt = inp %>% 
    bind_rows(op) %>% 
    bind_rows(ph) %>% 
    bind_rows(pharmacy) %>%
    mutate(disease = diseases[i]) %>% 
    group_by(disease, claim_type) %>% 
    summarise(cost = sum(cost),
              claim_count = n(),
              distinct_patient = n_distinct(z_patid)) %>% 
    mutate(total_patients = length(patients),
           total_cost = sum(cost),
           pp = total_cost/total_patients,
           pppy = pp/5)
  
  
  df_main = inp %>% 
    bind_rows(op) %>% 
    bind_rows(ph) %>%
    bind_rows(pharmacy) %>%
    group_by(z_patid) %>% 
    summarise(cost = sum(cost),
              claim_count = n(),
              date_diff = as.numeric(max(fst_dt)-min(fst_dt))) %>% 
    mutate(average_lag = date_diff/claim_count)
  
  kmdat = df_main %>% 
    select(cost, claim_count, average_lag) %>% 
    mutate(cost = log(cost)/log(max(cost)),
           claim_count = claim_count/max(claim_count),
           average_lag = average_lag/max(average_lag))
  km = kmeans(kmdat, 5, iter.max = 300)
  
  df_main$cluster = km$cluster
  df_main$disease = diseases[i]
  
  if(i==1){
    report = pt
    cluster_report = df_main
  }else{
    report = report %>% bind_rows(pt)
    cluster_report = cluster_report %>% bind_rows(df_main)
  }
  
}
Sys.time()-t
report_directory = paste(getwd(), "/", 'reports/cost_cluster/', sep = "")
write.csv(report, paste(report_directory,'cost_report_withoutICD.csv', sep = ''), row.names = FALSE)
write.csv(cluster_report, paste(report_directory,'cluster_report_withoutICD.csv', sep = ''), row.names = FALSE)
