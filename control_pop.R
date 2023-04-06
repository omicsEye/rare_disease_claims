library(tidyverse)
library(RODBC)


patient_list = read.csv('reports/prevalence/overall/patient_list.csv')
grouped_patients = patient_list %>% 
  group_by(age_band_cd) %>% 
  summarise(patients = n()) %>% 
  mutate(weights = patients/sum(patients))

ch = odbcConnect('vertica',
                 uid = 'XXXXX',
                 pwd = 'XXX')

inp = sqlQuery(ch, paste("select distinct(inp.z_patid)
from HCCI_2.INP_SDDV2 inp
where inp.proc_cd  in ('90750', '90751', '90752', '90753', '90754')
and inp.yr > '2015'
group by 1"))

op = sqlQuery(ch, paste("select distinct(op.z_patid)
from HCCI_2.OP_SDDV2 op
where op.proc_cd  in ('90750', '90751', '90752',  '90753','90754')
and op.yr > '2015'
group by 1"))

ph = sqlQuery(ch, paste("select distinct(ph.z_patid)
from HCCI_2.PHYS_SDDV2 ph
where ph.proc_cd  in ('90750', '90751', '90752', '90753', '90754')
and ph.yr > '2015'
group by 1"))


control_list = unique(unlist(c(inp$z_patid, op$z_patid, ph$z_patid)))

control_list = setdiff(control_list, patient_list$z_patid)

cut_points = round(seq(1,length(control_list), length.out = 5))
for( i in 1:(length(cut_points)-1)){
  print(i)
control_query = paste("'", control_list[cut_points[i]:cut_points[i+1]], "'", sep = "")

if(i ==1){
control_members = sqlQuery(ch, paste("select ss.z_patid, ss.sex, ss.age_band_cd, ss.mbr_zip_5_cd, ss.mbr_state
from(
SELECT *, 
ROW_NUMBER() over(PARTITION by z_patid ORDER by yr DESC) as rn
FROM HCCI_2.MBR_SDDV2 ms 
where ms.z_patid in (", paste(control_query, collapse = ','),
                             ")) as ss
where ss.rn = 1;"))
}else{
  tmp = sqlQuery(ch, paste("select ss.z_patid, ss.sex, ss.age_band_cd, ss.mbr_zip_5_cd, ss.mbr_state
from(
SELECT *, 
ROW_NUMBER() over(PARTITION by z_patid ORDER by yr DESC) as rn
FROM HCCI_2.MBR_SDDV2 ms 
where ms.z_patid in (", paste(control_query, collapse = ','),
                           ")) as ss
where ss.rn = 1;"))
  control_members = control_members %>% bind_rows(tmp)
}
}


grouped_control = control_members %>% 
  group_by(age_band_cd) %>% 
  summarise(patients = n())


patients = control_members$z_patid[sample(1:nrow(control_members), nrow(patient_list)*4)]
patients = paste("'", patients, "'", sep = "")

inp = sqlQuery(ch, paste("SELECT inp.z_patid,
SUM(inp.calc_allwd) as cost, 'Inpatient' as claim_type
from HCCI_2.INP_SDDV2 inp
where inp.yr > 2015 AND inp.z_patid in (", 
                         paste(patients, collapse = ','), ")
group by 1
having SUM(inp.calc_allwd) > 0;"))

op = sqlQuery(ch, paste("SELECT op.z_patid,
SUM(op.calc_allwd) as cost, 'Outpatient' as claim_type
from HCCI_2.OP_SDDV2 op
where op.yr > 2015 AND op.z_patid in (", 
                         paste(patients, collapse = ','), ")
group by 1
having SUM(op.calc_allwd) > 0;"))

ph = sqlQuery(ch, paste("SELECT ph.z_patid,
SUM(ph.calc_allwd) as cost, 'Physician' as claim_type
from HCCI_2.PHYS_SDDV2 ph
where ph.yr > 2015 AND ph.z_patid in (", 
                        paste(patients, collapse = ','), ")
group by 1
having SUM(ph.calc_allwd) > 0;"))

pharmacy = sqlQuery(ch, paste("SELECT rs.z_patid,
SUM(rs.calc_allwd) as cost, 'Pharmacy' as claim_type
from HCCI_2.RX_SDDV2 rs
where rs.yr > 2015 AND rs.z_patid in (", 
                              paste(patients, collapse = ','), ")
group by 1
having SUM(rs.calc_allwd) > 0;"))


df_main = inp %>% 
  bind_rows(op) %>% 
  bind_rows(ph) %>%
  bind_rows(pharmacy) %>%
  group_by(z_patid, claim_type) %>% 
  summarise(cost_patient = sum(cost)) 
# %>% 
  # left_join(control_members)

report = df_main %>% 
  group_by(claim_type) %>% 
  summarise(cost = sum(cost_patient),
            patients = n())


tmp = df_main %>% 
  group_by(age_band_cd) %>% 
  summarise(cost = sum(cost_patient),
            med = median(cost_patient),
            count = n()) %>% 
  left_join(grouped_patients)

tmp$cost_per_patient = tmp$cost/tmp$count
tmp$cost_group_weighted = tmp$cost_per_patient*tmp$patients
tmp$pp = sum(tmp$cost_group_weighted)/sum(tmp$patients)
tmp$pppy = tmp$pp/5
