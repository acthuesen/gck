library(tidyverse)
library(lubridate)

# set wd to place where all sensor files are
allsens <- list.files()

# load all files
allsensordata <- data.frame()
for(allsen in allsens){
  tempdata <- read.table(paste0(unique(allsen)), skip=3) #load
  tempdata$id <- allsen
  allsensordata <- rbind(allsensordata, tempdata) # merge
  rm(tempdata) # cleanup
}

# cleanup variables and dates, calculate daily sd, create next/prev obs
cgm <- allsensordata %>%
  group_by(id) %>%
  rename('date' = V2, 'time' = V3, 'value' = V5) %>%
  select(-c(V1,V4)) %>%
  mutate(date = ymd(date)) %>%
  mutate(datetime = ymd_hm(paste(date,time))) %>%
  mutate(time = as.POSIXct(time,format="%H:%M")) %>%
  mutate(consdate = date - first(date)) %>%
  separate(id, into=c('patid', 'seq', NA), sep='[_.]',remove=F) %>%
  mutate(nextvalue = shift(value, n=1, fill=NA, type='lead'),
         prevvalue = shift(value, n=1, fill=NA, type='lag')) %>%
  mutate(diff_nextvalue = value - nextvalue,
         diff_prevvalue = value - prevvalue) %>%
  group_by(consdate,id) %>%
  mutate(sdbg = sd(value)) %>%
  ungroup()

# remove obs's that have no difference between next/prev
cgm_filterednodiff <- cgm %>%
  group_by(id) %>%
  mutate(nodiff = case_when(diff_nextvalue==0 ~ 'nodiff',
                            TRUE ~ 'diff')) %>%
  filter(nodiff == 'diff') %>%
  select(-c(nodiff, nextvalue, prevvalue, diff_nextvalue, diff_prevvalue))

# calculate glycemic excursion, filter only those > daily sd
cgm_pv <- cgm_filterednodiff %>%
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
  filter(true_pv == 'true')

# write output
MAGEs <- cgm_pv %>%
  group_by(patid,seq) %>%
  summarise(MAGE = mean(diff_pv))
write.table(MAGEs, wd, row.names = F)
