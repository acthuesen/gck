# ------------ #
# Dependencies #
# ------------ #
library(tidyverse)
library(readxl)
library(lubridate)
library(data.table)

# --------- #
# OGTT data #
# --------- #
ogtt$var <- as.factor(ogtt$var)
ogtt$var <- fct_relevel(ogtt$var, "glu", "cpep", "ins")

# wide
ogtt_wide <- ogtt %>%
  group_by(var, project_id) %>%
  pivot_wider(names_from = c(var,time), values_from = val) %>%
  mutate(fpg = mean(glu_0,`glu_-5`, na.rm=T),
         fpi = mean(ins_0,`ins_-5`, na.rm=T),
         fpc = mean(cpep_0,`cpep_-5`, na.rm=T),
         dpg = mean(glu_30,glu_60,glu_90,glu_120, na.rm=T),
         dpi = mean(ins_30,ins_60,ins_90,ins_120, na.rm=T),
         dpc = mean(cpep_30,cpep_60,cpep_90,cpep_120, na.rm=T),
         glu_incr_2hr = glu_120-fpg,
         isi = 10000/sqrt(((fpg*18.018)*(fpi/6.945))*((dpg*18.018)*(dpi/6.945))),
         igi = (ins_30-fpi)/(glu_30-fpg))

# maximal glucose measurement
ogttmax <- ogtt %>%
  group_by(var, project_id) %>%
  summarise(max_val = max(val, na.rm=T),
            mean_val = mean(val, na.rm=T)) %>%
  pivot_wider(names_from = var, values_from = max_val:mean_val)

# merge
ogttagg <- reduce(list(ogtt_wide, ogttmax), left_join, by='project_id') %>%
  select(project_id,fpg:mean_val_ins) %>%
  mutate(glu_incr_max = max_val_glu - fpg)

rm(ogttmax, ogtt_wide)

# ---------------- #
# Biochemical data #
# ---------------- #
bc <- bc %>%
  mutate(HbA1cpc = HbA1c/10.929+2.15,
         HbA1cpc_3mo = HbA1c_3mo/10.929+2.15)

# ------------------- #
# Anthropometric data #
# ------------------- #
# aggregate multiple measurements
anthro <- anthro %>%
  rowwise() %>%
  mutate(weight = mean(c(weight_1,weight_2,weight_3), na.rm = T)) %>%
  mutate(height = mean(c(height_1,height_2,height_3), na.rm = T)) %>%
  mutate(wcmid = mean(c(wc_mid_1,wc_mid_2,wc_mid_3), na.rm = T)) %>%
  mutate(wcic = mean(c(wc_ic_1,wc_ic_2,wc_ic_2), na.rm = T)) %>%
  mutate(hip = mean(c(hip_1,hip_2,hip_3), na.rm = T)) %>%
  mutate(sysbp = mean(c(sys_bp_1,sys_bp_2,sys_bp_3), na.rm = T)) %>%
  mutate(diabp = mean(c(dia_bp_1,dia_bp_2,dia_bp_3), na.rm = T)) %>%
  mutate(pulse = mean(c(pulse_1,pulse_2,pulse_3), na.rm = T)) %>%
  select(-c(height_1:pulse_3)) %>%
  mutate(weight_3mo = mean(c(weight_1_3mo,weight_2_3mo,weight_3_3mo), na.rm = T)) %>%
  mutate(height_3mo = mean(c(height_1_3mo,height_2_3mo,height_3_3mo), na.rm = T)) %>%
  mutate(wcmid_3mo = mean(c(wc_mid_1_3mo,wc_mid_2_3mo,wc_mid_3_3mo), na.rm = T)) %>%
  mutate(wcic_3mo = mean(c(wc_ic_1_3mo,wc_ic_2_3mo,wc_ic_2_3mo), na.rm = T)) %>%
  mutate(hip_3mo = mean(c(hip_1_3mo,hip_2_3mo,hip_3_3mo), na.rm = T)) %>%
  mutate(sysbp_3mo = mean(c(sys_bp_1_3mo,sys_bp_2_3mo,sys_bp_3_3mo), na.rm = T)) %>%
  mutate(diabp_3mo = mean(c(dia_bp_1_3mo,dia_bp_2_3mo,dia_bp_3_3mo), na.rm = T)) %>%
  mutate(pulse_3mo = mean(c(pulse_1_3mo,pulse_2_3mo,pulse_3_3mo), na.rm = T)) %>%
  select(-c(height_1_3mo:pulse_3_3mo))

# create new vars
anthro <- anthro %>%
  mutate(bmi = weight/(height/100)^2) %>% # calculate BMI
  mutate(bmi_3mo = weight_3mo/(height_3mo/100)^2) %>%
  mutate(whr = wcmid/hip) %>%
  mutate(whr_3mo = wcmid_3mo/hip_3mo) %>%
  mutate(map = diabp+((1/3)*(sysbp-diabp))) %>% # calc MAP
  mutate(map_3mo = diabp_3mo+((1/3)*(sysbp_3mo-diabp_3mo))) %>% # calc MAP
  mutate(age = time_length(interval(DOB, exam_date), "years")) %>% # calc current age (interval between date and dob)
  mutate(age = round(age)) %>%
  mutate(duration = (age-diab_debut_age)) %>%
  select(-DOB) %>%
  mutate(dyslip = case_when(str_detect(comorbidities,'hypercholesterolemia') | (HDL < 1.2 & Sex == 'Female' | HDL < 1.0 & Sex == 'Male') | Triglyceride > 1.7 | LDL > 2.5 ~ 'Yes', 
                            TRUE ~ 'No')) %>%
  mutate(hypt = case_when(str_detect(comorbidities, 'hypertension') | sysbp > 140 | diabp > 90 ~ 'Yes',
                          TRUE ~ 'No'),
         first_deg_rel = case_when(fam_hist_bin == 1 ~ 'Yes',
                                   fam_hist_bin == 0 ~ 'No'),
         treat_current = case_when(treat_current == 'metformin' ~ 'Metformin', 
                                   treat_current == 'diet' ~ 'Diet',
                                   str_detect(treat_current,',') ~ 'Combination',
                                   TRUE ~ ''),
         bmi_cat = case_when(bmi >= 25 & bmi < 30 ~ 'Overweight',
                             bmi >= 30 ~ 'Obese',
                             bmi >= 18.5 & bmi < 25 ~ 'Normal weight',
                             bmi < 18.5 ~ 'Underweight'),
         ACMG3 = case_when(ACMG == 'P' | ACMG == 'LP' ~ 'P/LP',
                           ACMG == 'VUS' | ACMG == 'LB' | ACMG == 'B' ~ 'VUS/B'))

# ----------- #
# Sensor data #
# ----------- #
allsens <- list.files()

allsensordata <- data.frame()
for(allsen in allsens){
  tempdata <- read.table(paste0(unique(allsen)), skip=3) #load
  tempdata$id <- allsen
  allsensordata <- rbind(allsensordata, tempdata) # merge
  rm(tempdata) # cleanup
}

allsensordata <- allsensordata %>%
  group_by(id) %>%
  rename('date' = V2, 'time' = V3, 'value' = V5) %>%
  select(-c(V1,V4)) %>%
  mutate(date = ymd(date)) %>%
  mutate(datetime = ymd_hm(paste(date,time))) %>%
  mutate(time = as.POSIXct(time,format="%H:%M")) %>%
  mutate(consdate = date - first(date)) %>%
  separate(id, into=c('project_id', 'seq', NA), sep='[_.]',remove=F) %>%
  mutate(project_id = tolower(project_id))

MAGE <- allsensordata %>%
  mutate(nextvalue = shift(value, n=1, fill=NA, type='lead'),
         prevvalue = shift(value, n=1, fill=NA, type='lag')) %>%
  mutate(diff_nextvalue = value - nextvalue,
         diff_prevvalue = value - prevvalue) %>%
  group_by(consdate,id) %>%
  mutate(sdbg = sd(value)) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(nodiff = case_when(diff_nextvalue==0 ~ 'nodiff',
                            TRUE ~ 'diff')) %>%
  filter(nodiff == 'diff') %>%
  select(-c(nodiff, nextvalue, prevvalue, diff_nextvalue, diff_prevvalue)) %>%
  group_by(id) %>%
  mutate(nextvalue = shift(value, n=1, fill=NA, type='lead'),
         prevvalue = shift(value, n=1, fill=NA, type='lag')) %>%
  mutate(diff_nextvalue = value - nextvalue,
         diff_prevvalue = value - prevvalue) %>%
  mutate(pv = case_when(diff_nextvalue > 0 & diff_prevvalue > 0 ~ 'peak',
                        diff_nextvalue < 0 & diff_prevvalue < 0 ~ 'valley')) %>%
  select(-c(nextvalue, prevvalue, diff_nextvalue, diff_prevvalue)) %>%
  drop_na(pv) %>%
  mutate(adjacent_pvval = shift(value, n=1, fill=NA, type='lead')) %>%
  mutate(diff_pv = abs(value-adjacent_pvval)) %>%
  mutate(true_pv = case_when(diff_pv > sdbg ~ 'true',
                             TRUE ~ 'false')) %>%
  filter(true_pv == 'true') %>%
  group_by(project_id,seq) %>%
  summarise(MAGE = mean(diff_pv))

cgm_agg <- allsensordata %>%
  group_by(project_id, seq) %>%
  filter(consdate != 0) %>%
  summarise(glucosemean = mean(value),
            glucosesd = sd(value),
            glucosecv = sd(value)/mean(value)*100)

cgm_agg <- full_join(cgm_agg, MAGE, by=c('seq','project_id')) %>%
  pivot_wider(names_from = seq, values_from = glucosemean:MAGE)

#cleanup
rm(MAGE, allsen, allsens)
