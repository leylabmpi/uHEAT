####Import####
##Read Files - QC
seqstats1 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats2 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC_pool3/reports/final/seqkit_stats.tsv", header=TRUE) 

##Read Files - rel abund
#output of llmgp #update: instead, for re-analysis, read in the ids only or the filtered table only
#profile_dir <- "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken" #updated: no subsampling, GTDB map, all samples
# listing files
#brk_cls_files = list_files(profile_dir, 'all-combined-bracken.tsv') 
#brk_cls_files %>% length
# reading tables
#brk_cls = brk_cls_files %>%
#  plyr::llply(read_bracken) %>%
#  data.table::rbindlist(use.names=TRUE, idcol='dataset')
#brk_cls
#names(brk_cls)

#fix the double participant ID problem 
#this is very slow: what might be faster would be to pivot, and then fix the names, and then pivot back.
#H179 is the 2nd set for H051, 
#brk_cls$Sample <- gsub("H179_1", "H051_8", brk_cls$Sample)
#brk_cls$Sample <- gsub("H179_2", "H051_9", brk_cls$Sample)
#brk_cls$Sample <- gsub("H179_3", "H051_10", brk_cls$Sample)
#brk_cls$Sample <- gsub("H179_4", "H051_11", brk_cls$Sample)
#brk_cls$Sample <- gsub("H179_5", "H051_12", brk_cls$Sample)
#brk_cls$Sample <- gsub("H179_6", "H051_13", brk_cls$Sample)

#& H120 is considered the 2nd set for H182 (see metadata) #no vaccine, covid! exclude entirely
#brk_cls$Sample <- gsub("H120_1", "H182_8", brk_cls$Sample)
#brk_cls$Sample <- gsub("H120_2", "H182_9", brk_cls$Sample)
#brk_cls$Sample <- gsub("H120_3", "H182_10", brk_cls$Sample)
#brk_cls$Sample <- gsub("H120_4", "H182_11", brk_cls$Sample)
#brk_cls$Sample <- gsub("H120_5", "H182_12", brk_cls$Sample)
#brk_cls$Sample <- gsub("H120_6", "H182_13", brk_cls$Sample)

#re-parse participant IDs for later merging
#brk_cls$Participant_ID <- sapply(strsplit(brk_cls$Sample,"_"), `[`, 1)
#brk_cls$Participant_ID <- gsub("H", "HEAT_", brk_cls$Participant_ID)

#check the Abundance metric here, is it definitely relative?
#brk_cls %>% group_by(Sample) %>% summarize(sum(Abundance)) #yes! :)

#check the sample IDs. are they all here?
#ids <- brk_cls %>% select(Sample) %>% filter(!duplicated(Sample))
#ids$Participant_ID <- sapply(strsplit(ids$Sample,"_"), `[`, 1)
#ids$Participant_ID <- gsub("H", "HEAT_", ids$Participant_ID)
#ids %>% group_by(Participant_ID) %>% count() %>% View() #yes #looks very nice and complete. #in this version I kept 120
#write.csv(ids, "/ebio/abt3_projects2/uHEAT/data/metadata/ids.csv")
ids <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/ids.csv", row.names=1)

#keep a tax table
#tax <- brk_cls %>% select(Domain, Phylum, Class, Order, Family, Genus, Species, taxonomy_id) %>%
#  filter(!duplicated(taxonomy_id))
#write.csv(tax, "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken/tax_table.csv")

###Read the raw values (counts, not abund)###
#raw = brk_cls_files %>%
#  plyr::llply(read_bracken, keep_frac=FALSE) %>% #keeps the raw counts instead
#  data.table::rbindlist(use.names=TRUE, idcol='dataset')

#fix the double participant ID problem 
#H179 is the 2nd set for H051, 
#raw$Sample <- gsub("H179_1", "H051_8", raw$Sample)
#raw$Sample <- gsub("H179_2", "H051_9", raw$Sample)
#raw$Sample <- gsub("H179_3", "H051_10", raw$Sample)
#raw$Sample <- gsub("H179_4", "H051_11", raw$Sample)
#raw$Sample <- gsub("H179_5", "H051_12", raw$Sample)
#raw$Sample <- gsub("H179_6", "H051_13", raw$Sample)

#& H120 is considered the 2nd set for H182 (see metadata)
#raw$Sample <- gsub("H120_1", "H182_8", raw$Sample)
#raw$Sample <- gsub("H120_2", "H182_9", raw$Sample)
#raw$Sample <- gsub("H120_3", "H182_10", raw$Sample)
#raw$Sample <- gsub("H120_4", "H182_11", raw$Sample)
#raw$Sample <- gsub("H120_5", "H182_12", raw$Sample)
#raw$Sample <- gsub("H120_6", "H182_13", raw$Sample)

#add metadata
#meta <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_fever_combinedR_2023.01.22_doubles2_fevercircadian_72h.csv") #updated with z score and circadian rhythm new Dec 30
meta <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_fever_combinedR_2023.01.22_doubles2_fevercircadian_72h.csv")
meta <- meta %>% filter(Status.=="completed") %>%
  filter(!(Participant_ID%in%c("HEAT_179", "HEAT_120"))) #they are now doubled up as 051, 182
#instead, 120 should be kicked out entirely, they never got vaccine. just use 182.
extra_fever <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/fevers_comparewindows_2023.03.17.csv") #H013, H170 already kicked out as untrustworthy
meta_extra <- meta[,c(1:99)] %>% full_join(extra_fever) #so, keep the new temperature data only
#if I want more: go back and re-define it, and double check it again. So there.
meta_extra
meta <- meta_extra

#fix certain variables in the metadata
#alcohol - jeez excel, why
meta$Alcohol_category_per_week <- gsub("03-Feb", "2-3", meta$Alcohol_category_per_week)
meta$Alcohol_category_per_week <- gsub("09-Apr", "4-9", meta$Alcohol_category_per_week)
meta$Alcohol_category_per_week <- ifelse(meta$Alcohol=="Never", yes="Never", no=meta$Alcohol_category_per_week)
#quickly add variable: "booster fever"
#meta$fever_booster <- ifelse(!is.na(meta$fever2)&meta$Vaccine_StudyDose==1,
#                             yes=meta$fever2,
#                             no=meta$fever)
#redefine fever_aboveavg due to rounding issues #already done
#meta$fever_aboveavg <- ifelse(round(meta$fever, 3) > 0.800, "yes", "no") #strictly above 0.8
#summary(meta$fever1d)
#summary(meta$fever2d)
#summary(meta$fever3d)
#summary(meta$fever7d)

#fever diet interaction by z score
#meta$Diet_fever_interactionZ <- paste(meta$Diet2, meta$FZ_cat, sep="_")
#meta$Diet_fever_interactionZ <- ifelse(meta$Diet_fever_interaction %in% c("NA_NA", "Omni_NA", "Veg_NA"),
#                                      yes=NA, no=meta$Diet_fever_interaction)

#add IgG and CRP to the meta file
crp <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/CRP.csv", row.names=1)
igg <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/IgG_titres.csv", row.names=1)

meta <- Reduce(full_join, list(meta, crp, igg))

#add bss + technical metadata (plate, pool, time at room temp) #actually this also contains the bss so it's fine
batch <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/batch-and-technical-metadata.csv")
batch <- batch %>% rename(SampleID=sample_id, Date_Sampled=Date, Time_Sampled=Time) %>% select(-c(sample_no, sample_type, X, X.1, Participant_ID))
#and: sample o.d. as a further technical parameter
od <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/SampleID-OD_final_forR.csv")

#longitudinal info to calculate the 'preceding temperature' / current temp for each sample
tempstool_long2 <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/temperature_and_stool_data_cleaned_long3_Nov2022.csv") #merge error fixed
tempstool_long2$SampleID <- paste(gsub("HEAT_", "H", tempstool_long2$Participant_ID), "_", tempstool_long2$Stool_ID, sep="")

#I need to separate out the time stamp for the Stool Entry, and then calculate the time window minus 24h, and then take the mean temp
#also pull in the Relative Day so I can use it as a tag for stool samples
stool_timestamps <- tempstool_long2 %>% filter(Entry=="Stool") %>% 
  group_by(SampleID) %>%
  select(Participant_ID, SampleID, Time0, Time01, Date) #Relative Day is NA for the Stools but I want to have Date info to include it
stool_timestamps$Time_previous <- stool_timestamps$Time0 - 1 #because Time0 is in days #so the previous time is 1 day before
stool_timestamps <- stool_timestamps %>% filter(!duplicated(SampleID)) %>% rename(Time0_stool=Time0,
                                                                                  Time01_stool=Time01,
                                                                                  Time_previous_stool=Time_previous)
stool_timestamps
stool_timestamps$SampleID

#pull the relative day for later, otherwise I will lose it
stool_relative_day <- tempstool_long2 %>% filter(Entry=="Temperature") %>% 
  select(Participant_ID, Temperature, Time0, Time01, Date.YY.MM.DD, Relative_Day) %>% 
  rename(Date=Date.YY.MM.DD) %>% 
  full_join(stool_timestamps) %>% 
  filter(!duplicated(SampleID)) %>% 
  select(SampleID, Relative_Day, Time0_stool, Time01_stool) %>% #also keep time zero for plotting
  rename(Time0=Time0_stool, Time01=Time01_stool) 
stool_relative_day 

#and now that it's clear per sample, clean up the doubles.
stool_relative_day$SampleID <- gsub("_v2", "", stool_relative_day$SampleID, fixed=TRUE)
stool_relative_day$SampleID <- gsub(".2", "", stool_relative_day$SampleID, fixed=TRUE)
#H179 is the 2nd set for H051, 
stool_relative_day$SampleID <- gsub("H179_1", "H051_8", stool_relative_day$SampleID)
stool_relative_day$SampleID <- gsub("H179_2", "H051_9", stool_relative_day$SampleID)
stool_relative_day$SampleID <- gsub("H179_3", "H051_10", stool_relative_day$SampleID)
stool_relative_day$SampleID <- gsub("H179_4", "H051_11", stool_relative_day$SampleID)
stool_relative_day$SampleID <- gsub("H179_5", "H051_12", stool_relative_day$SampleID)
stool_relative_day$SampleID <- gsub("H179_6", "H051_13", stool_relative_day$SampleID)

#define previous temp
previous_temp <- tempstool_long2 %>% filter(Entry=="Temperature") %>% 
  select(Participant_ID, Temperature, Time0, Time01) %>% #Date.YY.MM.DD, Relative_Day) %>% #include both date and Relative Day #this should mean the Relative Day is OK for doubles
  #rename(Date=Date.YY.MM.DD) %>% #however by keeping Date and Relative Day... I lose the 'previous' part to merge...
  full_join(stool_timestamps) #joins by ID only
previous_temp


#some samples are missing mainly due to odd days, actual missing dates, or transcription errors #now fixed?
missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
missing_days$Time01 <- missing_days$Time0 - 7.5 #assume they got their impfung in the exact middle. 
stool_relative_day <- full_join(stool_relative_day, missing_days) %>% 
  select(SampleID, Relative_Day, Time0, Time01) #note however, missing Time0 for the new relative days
stool_relative_day <- stool_relative_day %>% filter(!is.na(SampleID)) %>% filter(!is.na(Relative_Day))
#are any now duplicated
dim(stool_relative_day) #it's rather long #also, 3 samples are apparently duplicated
head(stool_relative_day)
tail(stool_relative_day)
#stool_relative_day %>% filter(duplicated(SampleID))
#stool_relative_day %>% filter(SampleID=="H118_4") #H033_9, H104_3, H118_4
stool_relative_day <- stool_relative_day %>% filter(!duplicated(SampleID))

#define previous temp #so the window is: temperatures taken before or at the stool time, but no longer than a day ago
previous_temp$window <- ifelse(previous_temp$Time0<=previous_temp$Time0_stool & previous_temp$Time0>previous_temp$Time_previous_stool,
                               "yes", "no")
previous_temp$window_recent <- ifelse(previous_temp$Time0<=previous_temp$Time0_stool & previous_temp$Time0>(previous_temp$Time0_stool - 0.5), #half a day instead of a day
                                      "yes", "no")
previous_temp %>% count(SampleID, window, window_recent) #now we found all of the temp samples per person #except that it is double counted for the doubles like H015 #so maybe keep the separate IDs for now..
previous_temp <- previous_temp %>% filter(window=="yes"&!is.na(Temperature))
previous_temp <- previous_temp %>% 
  group_by(SampleID) %>% 
  summarize(Previous_Stool_Temp=mean(Temperature),
            Previous_Stool_Temp_max=max(Temperature),
            Previous_Stool_Temp_recent=mean(Temperature[window_recent=="yes"]))
previous_temp
hist(previous_temp$Previous_Stool_Temp) 
hist(previous_temp$Previous_Stool_Temp_max) 
hist(previous_temp$Previous_Stool_Temp_recent) 

#and now that it's clear per sample, clean up the doubles.
previous_temp$SampleID <- gsub("_v2", "", previous_temp$SampleID, fixed=TRUE)
previous_temp$SampleID <- gsub(".2", "", previous_temp$SampleID, fixed=TRUE)
#H179 is the 2nd set for H051, 
previous_temp$SampleID <- gsub("H179_1", "H051_8", previous_temp$SampleID)
previous_temp$SampleID <- gsub("H179_2", "H051_9", previous_temp$SampleID)
previous_temp$SampleID <- gsub("H179_3", "H051_10", previous_temp$SampleID)
previous_temp$SampleID <- gsub("H179_4", "H051_11", previous_temp$SampleID)
previous_temp$SampleID <- gsub("H179_5", "H051_12", previous_temp$SampleID)
previous_temp$SampleID <- gsub("H179_6", "H051_13", previous_temp$SampleID)

#missing sample temps e.g. 009_1, 013_6. why? #generally because there is actually no temp information in the right window
previous_temp %>% filter(duplicated(SampleID)) #good.

#extra info: drugs, symptoms.
drugs_daily <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/drug_info_all_long.csv")
symptoms <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/symptoms_info.csv", row.names=1)

#similar to the temperature: I would need to associate drugs with a particular stool sample or not
drugs_timing <- drugs_daily %>% select(-c(Time_of_Day_v2, Time_hours_drugs, Time_minutes_drugs)) %>%
  full_join(stool_timestamps) %>% 
  filter(!is.na(Medication))

drugs_timing$window <- ifelse(drugs_timing$Time0<=drugs_timing$Time0_stool & drugs_timing$Time0>drugs_timing$Time_previous_stool,
                              "yes", "no")
drugs_timing <- drugs_timing %>% filter(window=="yes")
drugs_timing <- drugs_timing %>% select(SampleID, Medication, Name_of_Medication, Type_of_Medication, 
                                        Antipyretic, Dose_mg)
drugs_timing <- drugs_timing %>% rename(Medication_24h_prior=Medication)
drugs_timing #this is now a cleaned table of samples with a drug taken in the 24h prior to them
#H030_4 comes up twice for two different NSAIDs, which isn't really useful for the question
drugs_timing <- drugs_timing %>% filter(!duplicated(SampleID))

#and conversely to make a cleaner summary, e.g. with antipyretics or not...
drugs_summary <- drugs_daily %>% group_by(Participant_ID) %>%
  summarise(Total_Meds_Taken=sum(Medication=="yes"),
            Total_VaccineAntipyretics_Taken=sum(Antipyretic_vaccine=="yes"),
            Total_Antipyretics_Taken=sum(Antipyretic=="yes"),
            Total_NSAID_Taken=sum(Type_of_Medication=="NSAID")
  )
drugs_summary$vaccine_antipyretics <- ifelse(drugs_summary$Total_VaccineAntipyretics_Taken>0, "yes", "no")

#add sequencing depth as a variable for the meta_seq table
#seq_depth <- raw %>% group_by(Sample) %>% 
#  summarize(seq_depth=sum(Abundance)) %>%
#  as.data.frame() %>% 
#  rename(SampleID=Sample)

#or: take it directly from the seq_stats table
seqstats <- full_join(seqstats1, seqstats2) %>% 
  filter(Read==1) %>%
  rename(SampleID=Sample, seq_depth=num_seqs) %>% 
  select(SampleID, seq_depth) %>%
  mutate(seq_depth=seq_depth*2)
seqstats$SampleID <- gsub("H179_1", "H051_8", seqstats$SampleID)
seqstats$SampleID <- gsub("H179_2", "H051_9", seqstats$SampleID)
seqstats$SampleID <- gsub("H179_3", "H051_10", seqstats$SampleID)
seqstats$SampleID <- gsub("H179_4", "H051_11", seqstats$SampleID)
seqstats$SampleID <- gsub("H179_5", "H051_12", seqstats$SampleID)
seqstats$SampleID <- gsub("H179_6", "H051_13", seqstats$SampleID)

#clean and merge a combined metadata file
meta_seq_raw <- ids #raw %>% select(Sample) %>% group_by(Sample) %>% filter(!duplicated(Sample))
meta_seq_raw <- meta_seq_raw %>% rename(SampleID=Sample) %>% as.data.frame()
#plus, add actual metadata to combine 
meta_seq <- Reduce(full_join, list(meta_seq_raw, meta, batch, od))
#and seq depth
meta_seq <- full_join(meta_seq, seqstats)
#and then add; previous temperature (previous 24h) & Relative Day
meta_seq <- Reduce(full_join, list(meta_seq, previous_temp, stool_relative_day))
#and drugs and symptoms
meta_seq <- Reduce(full_join, list(meta_seq, drugs_timing, drugs_summary, symptoms))
#cleanup
meta_seq <- meta_seq %>%
  filter(!is.na(meta_seq$SampleID)&!is.na(meta_seq$Participant_ID)&!(meta_seq$Participant_ID %in% c("HEAT_179", "HEAT_120")))
row.names(meta_seq) <- meta_seq$SampleID

#expected number of samples per participant?
dplyr::count(meta_seq, Participant_ID) #nice nice nice
#anything weird with : missing Relative Day, missing ID, etc.
meta_seq %>% filter(is.na(SampleID))
meta_seq %>% filter(is.na(Participant_ID))

#define a few new metadata variables to use later
#better BSS definition
meta_seq$Bristol_Stool_Score2 <- case_when(meta_seq$Bristol_Stool_Score%in%c('1', '1.5', '2', '2.5') ~ "1-2.5",
                                           meta_seq$Bristol_Stool_Score%in%c('3', '3.5', '4', '4.5') ~ "3-4.5",
                                           meta_seq$Bristol_Stool_Score%in%c('5', '5.5', '6', '6.5', '7') ~ "5-7") 

#fix missing Stool ID for the double participants
meta_seq$Stool_ID <- ifelse(!is.na(meta_seq$Participant_ID2)&is.na(meta_seq$Stool_ID), #double participants who are missing stool ID
                            yes=sapply(strsplit(meta_seq$SampleID,"_"), `[`, 2),
                            no=meta_seq$Stool_ID)

#Timing variables.
meta_seq %>% filter(is.na(Relative_Day)) %>% select(SampleID) #only weird ones are missing
meta_seq %>% filter(is.na(Time01)) %>% select(SampleID) #slightly more are missing. 

#relative day numeric
meta_seq$Relative_Day_numeric <- gsub("t+", "", meta_seq$Relative_Day)
meta_seq$Relative_Day_numeric <- gsub("t", "", meta_seq$Relative_Day_numeric)
meta_seq$Relative_Day_numeric <- as.numeric(meta_seq$Relative_Day_numeric)
#If the sample has a relative day but no timestamp, impute Time01 as relative day
summary(meta_seq$Relative_Day_numeric)
summary(meta_seq$Time01)
meta_seq$Time01 <- ifelse(is.na(meta_seq$Time01), yes=meta_seq$Relative_Day_numeric, no=meta_seq$Time01)
meta_seq %>% filter(is.na(Time01)) %>% select(SampleID, Relative_Day, Time0, Time01) #only weird ones missing as intended

#timing variable #define it the same way as fever. #try to capture all possible acute samples.
meta_seq$Stool_timing1 <- ifelse(meta_seq$Relative_Day %in% c("t0", "t+1", "t+2") & meta_seq$Time01>0, yes="acute", no="baseline")
meta_seq$Stool_timing2 <- case_when(meta_seq$Relative_Day%in%c("t-7", "t-8", "t-6", "t-5", "t-4", "t-3", "t-2", "t-1", "t0") & meta_seq$Time01<=0 ~ "before",
                                    meta_seq$Relative_Day %in%c("t0", "t+1", "t+2") & meta_seq$Time01>0  ~ "acute",
                                    meta_seq$Relative_Day%in%c("t+3", "t+4", "t+5", "t+6", "t+7", "t+8") ~ "after") 
meta_seq$Stool_timing3 <- ifelse(meta_seq$Time01<=0, "before", "after") 

#assign the correct average temperature / fever data per stool ID (series 1, series 2)
meta_seq$sample_no <- as.numeric(sapply(strsplit(meta_seq$SampleID,"_"), `[`, 2))
meta_seq$fever7d1 <- meta_seq$fever7d
meta_seq$fever7d <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever7d2),
                           yes=meta_seq$fever7d2, no=meta_seq$fever7d)
meta_seq$fever2d1 <- meta_seq$fever2d
meta_seq$fever2d <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever2d2),
                           yes=meta_seq$fever2d2, no=meta_seq$fever2d)

#meta_seq$fever481 <- meta_seq$fever48
#meta_seq$fever48 <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever482),
#                         yes=meta_seq$fever482, no=meta_seq$fever48)

#do the same for baseline & other temperatures
meta_seq$median_base1 <- meta_seq$median_base
meta_seq$median_base <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$median_base2),
                               yes=meta_seq$median_base2, no=meta_seq$median_base)

meta_seq$max_48h_after1 <- meta_seq$max_48h_after
meta_seq$max_48h_after <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$max_48h_after2),
                                 yes=meta_seq$max_48h_after2, no=meta_seq$max_48h_after)

meta_seq$dummy_fever1 <- meta_seq$fever7db #the entire week before
meta_seq$dummy_fever <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever7db2),
                               yes=meta_seq$fever7db2, no=meta_seq$fever7db)


#meta_seq$fever_aboveavg1 <- meta_seq$fever_aboveavg
#meta_seq$fever_aboveavg <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever_aboveavg2),
#                               yes=meta_seq$fever_aboveavg2, no=meta_seq$fever_aboveavg)

#meta_seq$fever_zscore1 <- meta_seq$fever_zscore
#meta_seq$fever_zscore <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$fever_zscore2),
#                                  yes=meta_seq$fever_zscore2, no=meta_seq$fever_zscore)

#meta_seq$FZ_cat1 <- meta_seq$FZ_cat
#meta_seq$FZ_cat <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$FZ_cat2),
#                                yes=meta_seq$FZ_cat2, no=meta_seq$FZ_cat)

#meta_seq$FZ_cat_lmh1 <- meta_seq$FZ_cat_lmh
#meta_seq$FZ_cat_lmh <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$FZ_cat_lmh2),
#                          yes=meta_seq$FZ_cat_lmh2, no=meta_seq$FZ_cat_lmh)

#comfortingly, this makes NO DIFFERENCE to summary fever category because everyone was consistent across both doses
#whichever definition I pick #update: new z score difference for 053
meta_seq %>% filter(!is.na(fever7d2)) %>% select(SampleID, 
                                                 fever7d, fever7d1, fever7d2, 
                                                 fever7d_aboveavg)

#now, define the fever categories here based on correct fever per sample
#the third quartile doesn't change with time. Maybe it's the better category. 
meta_seq$fever7d_aboveavg <- ifelse(meta_seq$fever7d>median(meta_seq$fever7d, na.rm=TRUE), "hi", "lo") #change to median to have more even distribution?
meta_seq$fever7d_aboveavg <- factor(meta_seq$fever7d_aboveavg, levels=c("lo", "hi"))
dplyr::count(meta_seq, fever7d_aboveavg)
meta_seq$fever2d_aboveavg <- ifelse(meta_seq$fever2d>median(meta_seq$fever2d, na.rm=TRUE), "hi", "lo")
meta_seq$fever2d_aboveavg <- factor(meta_seq$fever2d_aboveavg, levels=c("lo", "hi"))
dplyr::count(meta_seq, fever2d_aboveavg)

meta_seq$fever7d_quart <- ifelse(meta_seq$fever7d>quantile(meta_seq$fever7d, na.rm=TRUE)[4], "hi", "lo")
meta_seq$fever7d_quart <- factor(meta_seq$fever7d_quart, levels=c("lo", "hi"))
dplyr::count(meta_seq, fever7d_quart)

#and: fever diet interaction #need to decide if we use 7d or 2d or what?
meta_seq$Diet_fever_interaction <- paste(meta_seq$Diet2, meta_seq$fever7d_aboveavg, sep="_")
meta_seq$Diet_fever_interaction <- ifelse(meta_seq$Diet_fever_interaction %in% c("NA_NA", "Omni_NA", "Veg_NA"),
                                          yes=NA, no=meta_seq$Diet_fever_interaction)

#Similarly: make corrected per-sample for IgG doubles (CRP I use the avg so it doesn't really matter)
meta_seq %>% filter(Second_diary_participation!="") %>% select(SampleID, IgG_B1, IgG_B2, IgG_B3)
meta_seq$IgG_B1 <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$IgG_B3),
                          yes=meta_seq$IgG_B2, no=meta_seq$IgG_B1) #if it's a second series with new blood, replace 1 with 2
meta_seq$IgG_B2 <- ifelse(meta_seq$sample_no>7 & !is.na(meta_seq$IgG_B3),
                          yes=meta_seq$IgG_B3, no=meta_seq$IgG_B2) #if it's a second series with new blood, replace 2 with 3
meta_seq %>% filter(Second_diary_participation!="") %>% select(SampleID, IgG_B1, IgG_B2, IgG_B3) #right
#and then: for 051 / 179
meta_seq %>% filter(Participant_ID=="HEAT_051") %>% select(SampleID, IgG_B1, IgG_B2, IgG_B3)
meta_seq[which(meta_seq$Participant_ID=="HEAT_051"&meta_seq$sample_no>7),]$IgG_B1 <- c(igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B1)
meta_seq[which(meta_seq$Participant_ID=="HEAT_051"&meta_seq$sample_no>7),]$IgG_B2 <- c(igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B2)
meta_seq %>% filter(Participant_ID=="HEAT_051") %>% select(SampleID, IgG_B1, IgG_B2, IgG_B3) #right

#now: re-define correct antibody categories?
meta_seq$IgG_B1_aboveavg <- ifelse(meta_seq$IgG_B1>median(meta_seq$IgG_B1, na.rm=TRUE), "hi", "lo") #mean or median??
dplyr::count(meta_seq, IgG_B1_aboveavg)
meta_seq$IgG_B2_aboveavg <- ifelse(meta_seq$IgG_B2>median(meta_seq$IgG_B2, na.rm=TRUE), "hi", "lo") #mean or median??
dplyr::count(meta_seq, IgG_B2_aboveavg)

#&: re-define correct immune ratios?
meta_seq$Immune_response <- meta_seq$IgG_B2 / meta_seq$fever7d
meta_seq$Immune_response_aboveavg <- ifelse(meta_seq$Immune_response>median(meta_seq$Immune_response, na.rm=TRUE), "hi", "lo") #mean or median??
dplyr::count(meta_seq, Immune_response_aboveavg)

#a 'relative temperature' (minus baseline) for the stool temp
meta_seq$Previous_Stool_Temp_Relative <- meta_seq$Previous_Stool_Temp - meta_seq$median_base
meta_seq$Previous_Stool_Temp_max_Relative <- meta_seq$Previous_Stool_Temp_max - meta_seq$median_base
meta_seq$Previous_Stool_Temp_Relative_recent <- meta_seq$Previous_Stool_Temp_recent - meta_seq$median_base

#define a "fever plus symptoms"
#meta_seq$fever_symptoms <- ifelse(meta_seq$fever>mean(meta_seq$fever[!is.na(meta_seq$fever)]) & meta_seq$num_symptoms>0,
#                                  "yes", "no")
#meta_seq$feverz_symptoms <- ifelse(meta_seq$fever_zscore48>mean(meta_seq$fever_zscore48[!is.na(meta_seq$fever_zscore48)]) & meta_seq$num_symptoms>0,
#                                  "yes", "no")
#conversely define "any symptoms"
meta_seq$Symptoms_Any <- ifelse(meta_seq$num_symptoms==0 & meta_seq$fever7d_quart=="lo",
                                "no", "yes")

#define a time of day per sample
meta_seq$Time_Sampled_numeric <- as.numeric(gsub(":", ".", meta_seq$Time_Sampled, fixed=TRUE))
meta_seq$Time_of_Day <- case_when(meta_seq$Time_Sampled_numeric<12~1,
                                  meta_seq$Time_Sampled_numeric>=12&meta_seq$Time_Sampled_numeric<18~2,
                                  meta_seq$Time_Sampled_numeric>=18~3)


#write out the final meta_seq so that it is easier to read in afterwards
#write.csv(meta_seq, "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_23.01.20_72h.csv")

#and if so: read it back in
#meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_23.01.20.csv", row.names=1)

####QC####
##first QC - read depth
#pools 1-3 combined
seqstats <- full_join(seqstats1, seqstats2)
seqstats$sample_type <- seqstats$Sample
seqstats[grep("blank", seqstats$sample_type),]$sample_type <- c("blank")
seqstats[grep("mock", seqstats$sample_type),]$sample_type <- c("mock")
seqstats$sample_type <- ifelse(seqstats$sample_type%in%(c("blank", "mock")),
                               yes=seqstats$sample_type, no="feces")
seqstats$sample_type
seqstats$sample_no <- sapply(strsplit(seqstats$Sample,"_"), `[`, 2)
seqstats$sample_type2 <- ifelse(seqstats$sample_type=="feces",
                                yes=paste(seqstats$sample_type, seqstats$sample_no, sep=""),
                                no=seqstats$sample_type)
seqstats_filt <- seqstats %>% filter(Read==1) #already done
seqstats_filt$num_seqs2 <- seqstats_filt$num_seqs*2 #already done

pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots/seq_depth_histogram_pools1-3_postqc.pdf")
p <- ggplot(seqstats, aes(x=num_seqs, fill=sample_type, colour=sample_type)) + geom_histogram(alpha=0.1)
p <- p + xlab("Number of paired reads")
p <- p + theme_bw()
p
dev.off()

pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots2/seq_depth_histogram_pools1-3_postqc_cutoffs2.pdf")
p <- ggplot(seqstats, aes(x=num_seqs, fill=sample_type, colour=sample_type)) + geom_histogram(alpha=0.1)
p <- p + xlab("Number of paired reads") + ylab("count")
p <- p + annotate("segment", x=500000, xend=500000, y=0, yend=185, linetype="dashed")
p <- p + annotate("text", label="0.5M \n exclude", x=500000, y=200)
p <- p + annotate("segment", x=5000000, xend=5000000, y=0, yend=185, linetype="dashed")
p <- p + annotate("text", label="5M \n rarefy", x=5000000, y=200)
p <- p + annotate("segment", x=10000000, xend=10000000, y=0, yend=185, linetype="dashed")
p <- p + annotate("text", label="10M \n rarefy2", x=9000000, y=200)
p <- p + annotate("segment", x=20000000, xend=20000000, y=0, yend=185, linetype="dashed")
p <- p + annotate("text", label="20M \n subsample", x=20000000, y=200)
p <- p + theme_bw()
p
dev.off()

#consider: are there any samples that should be filtered out entirely at this stage?
exclude <- seqstats_filt %>% filter(num_seqs<500000) %>% select(Sample) ##all of the blanks, plus 4 samples (some might actually be blank)
exclude

#for MAGS building: stop here to re-make the sample table excluding blanks, mocks and almost-blank samples; and incomplete series
samples_table <- read.delim("/ebio/abt3_projects2/uHEAT/data/output_kraken_all/samples.txt")
samples_table
dim(samples_table)
samples_table_filt <- samples_table %>% filter(!(Sample %in% exclude$Sample))
dim(samples_table_filt) #and now look for some misplaced samples e.g. mocks and incomplete series
exclude_mags <- samples_table_filt[grep(c("blank|mock|pilot|H012|H024|H025|H026|H032|H047|H063|H068|H124|H148"), samples_table_filt$Sample),]$Sample
exclude_mags

samples_table_filt <- samples_table_filt %>% filter(!(Sample %in% exclude_mags))
dim(samples_table)
dim(samples_table_filt)
head(samples_table_filt)
tail(samples_table_filt)

write.table(samples_table_filt, "/ebio/abt3_projects2/uHEAT/data/MAGS/samples_table.txt", sep="\t", quote=FALSE, row.names=FALSE)

#also exclude all the 120 samples for vaccine-related analyses
exclude <- c(exclude$Sample, c("H120_1", "H120_2", "H120_3", "H120_4", "H120_5", "H120_6"))
exclude
#and additionally: exclude samples from infected individuals; and other samples excessively old; and anyone without vaccine
infected <- meta_seq %>% filter(COVID_infection_during_study%in%c("y", "maybe"))
old <- meta_seq %>% filter(Time_at_RT>35)
no_vaccine <- meta_seq %>% filter(Vaccine_StudyDose=="no vaccine")

exclude_complete <- c(exclude, infected$SampleID, old$SampleID, no_vaccine$SampleID)
exclude_complete

#write out a table of the excluded samples, so I remember why I excluded them.
exclude1 <- seqstats_filt %>% filter(num_seqs<500000) %>% select(Sample) %>% rename(SampleID=Sample)##all of the blanks, plus 4 samples (some might actually be blank)
exclude1$Reason <- c("Low_depth_<0.5M")
exclude2 <- data.frame(SampleID=c("H120_1", "H120_2", "H120_3", "H120_4", "H120_5", "H120_6"),
                       Reason=c("No_vaccine_duplicate"))
exclude2
meta_seq_old <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_FILT_22.11.21.csv")
infected <- meta_seq_old %>% filter(COVID_infection_during_study%in%c("y", "maybe")) %>% 
  select(SampleID) %>% mutate(Reason=c("infected")) #right 
old <- meta_seq_old %>% filter(Time_at_RT>35) %>% 
  select(SampleID) %>% mutate(Reason=c("old_RT>35"))
no_vaccine <- meta_seq_old %>% filter(Vaccine_StudyDose=="no vaccine") %>% 
  select(SampleID) %>% mutate(Reason=c("no_vaccine"))

excluded_all <- Reduce(full_join, list(exclude1, exclude2, infected, old, no_vaccine))
excluded_all$Excluded <- c("Excluded")
excluded_all

meta_seq_exclude <- meta_seq_old %>% full_join(excluded_all) %>% 
  filter(SampleID %in% excluded_all$SampleID) %>%
  filter(!duplicated(SampleID))

write.csv(meta_seq_exclude, "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_EXCLUDED.csv")

###filter all tables for read depth, before proceeding.
exclude_complete <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_EXCLUDED.csv")
meta_seq <- meta_seq %>% filter(!(SampleID %in% exclude_complete$SampleID))

#also exclude mocks, pilot samples, etc. for the rest of the analysis
meta_seq <- meta_seq %>% filter(Status.=="completed")
#brk_cls <- brk_cls %>% filter(Sample %in% meta_seq$SampleID)

#write out the CLEAN FILTERED COPIES
#write.csv(meta_seq, "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_FILT_23.01.20.csv")
write.csv(meta_seq, "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_FILT_23.05.23.csv")
#write.csv(brk_cls, "/ebio/abt3_projects2/uHEAT/data/brk_cls_FILT_23.01.20.csv")

#*##second QC - unclassified reads
brk_uncls_files = list_files(profile_dir, 'all-combined_kraken-unclassified.tsv') 
unclassified <- read.table(brk_uncls_files, header=TRUE)
head(unclassified)
hist(unclassified$percent_reads)

#re-parse participant IDs for later merging
unclassified$Participant_ID <- sapply(strsplit(unclassified$sample,"_"), `[`, 1)
unclassified$Participant_ID <- gsub("H", "HEAT_", unclassified$Participant_ID)
unclassified$sample_no <- sapply(strsplit(unclassified$sample,"_"), `[`, 2)
unclassified$sample_type <- ifelse(unclassified$Participant_ID %in% c("mock", "blank"),
                                   yes=unclassified$Participant_ID,
                                   no="feces")


#make a nicer histogram
pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots/unclassified_histogram_pools1-3_postqc_withSTRUOGTDB.pdf")
p <- ggplot(unclassified, aes(x=percent_reads, fill=sample_type, colour=sample_type)) + geom_histogram(alpha=0.1)
p <- p + xlab("% of unclassified reads")
p <- p + theme_bw()
p
dev.off()

##QC - final read distribution of mapped reads
#check the Abundance metric here, it is raw, but not rarefied it seems
counts <- raw %>% group_by(Sample) %>% summarize(Depth=sum(Abundance)) %>% as.data.frame()#raw but not rarefied?
summary(counts$Depth)
hist(counts$Depth) #now with GTDB mapping we are back to similar read counts to the original

exclude2 <- counts %>% filter(Depth<100000) %>% select(Sample) #it would be the same samples as above


#QC: does the read depth change systemically with e.g. fever or relative day?
wilcox.test(meta_seq$seq_depth~meta_seq$fever_aboveavg) #NS
cor.test(meta_seq$seq_depth, meta_seq$fever, method='spearman') #NS
cor.test(meta_seq$seq_depth, meta_seq$Relative_Day_numeric, method='spearman') #NS

ggplot(meta_seq, aes(x=fever_aboveavg, y=seq_depth)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.1, height=0)
ggplot(meta_seq, aes(x=fever, y=seq_depth)) + geom_point()
ggplot(meta_seq, aes(x=Relative_Day_numeric, y=seq_depth)) + geom_point() #okay.

pdata <- meta_seq %>% filter(!is.na(fever_aboveavg)&Relative_Day_numeric<8&Relative_Day_numeric>-8)
pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots2/seq_depth_by_fever_overtime.pdf")
p <- ggplot(pdata, aes(x=Relative_Day_numeric, y=seq_depth, colour=fever_aboveavg)) + geom_point() + geom_smooth()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.title.y = element_text(margin = margin(r = 10)))
p <- p + xlab("Day Relative to Vaccine")
p <- p + scale_x_continuous(breaks=c(seq(-8, 8, by=1)))
p <- p + scale_colour_manual(values=c("grey50", "firebrick3"))
p
dev.off()
#OH SHIT

pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots2/seq_depth_by_fever_Time01.pdf")
p <- ggplot(pdata, aes(x=Time01, y=seq_depth, colour=fever_aboveavg)) + geom_point() + geom_smooth()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.title.y = element_text(margin = margin(r = 10)))
p <- p + xlab("Day Relative to Vaccine")
p <- p + scale_x_continuous(breaks=c(seq(-8, 8, by=1)))
p <- p + scale_colour_manual(values=c("grey50", "firebrick3"))
p
dev.off()

pdata1 <- pdata %>% filter(Relative_Day%in%c("t+1"))
cor.test(pdata1$seq_depth, pdata1$fever, method='spearman')
#well, shit #is it strong?

p <- ggplot(pdata1, aes(x=seq_depth, y=fever)) + geom_point() + geom_smooth(method='lm')
p
#I mean barely. I don't really believe that.

#for comparison, check: sample O.D. instead of seq depth (they are correlated with each other)
wilcox.test(meta_seq$OD~meta_seq$fever_aboveavg) #NS
cor.test(meta_seq$OD, meta_seq$fever, method='spearman') #NS
cor.test(meta_seq$OD, meta_seq$Relative_Day_numeric, method='spearman') #NS

ggplot(meta_seq, aes(x=fever_aboveavg, y=OD)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.1, height=0)
ggplot(meta_seq, aes(x=fever, y=OD)) + geom_point()
ggplot(meta_seq, aes(x=Relative_Day_numeric, y=OD)) + geom_point() #okay.

#plot 1 #hey it's not confounded yayyyy! :) 
p <- ggplot(pdata, aes(x=Relative_Day_numeric, y=OD, colour=fever_aboveavg)) + geom_point() + geom_smooth()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.title.y = element_text(margin = margin(r = 10)))
p <- p + xlab("Day Relative to Vaccine")
p <- p + scale_x_continuous(breaks=c(seq(-8, 8, by=1)))
p <- p + scale_colour_manual(values=c("grey50", "firebrick3"))
p

pdf(width=6, height=5, "/ebio/abt3_projects2/uHEAT/data/plots2/od_by_fever_Time01.pdf")
p <- ggplot(pdata, aes(x=Time01, y=log(OD), colour=fever_aboveavg)) + geom_point() + geom_smooth()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.title.y = element_text(margin = margin(r = 10)))
p <- p + xlab("Day Relative to Vaccine")
p <- p + scale_x_continuous(breaks=c(seq(-8, 8, by=1)))
p <- p + scale_colour_manual(values=c("grey50", "firebrick3"))
p
dev.off()

#read in the samples.txt table and make one that specifically includes only 'complete' samples for next rarefying q's
samples1 <- read.delim("/ebio/abt3_projects2/uHEAT/data/output_kraken_all/samples.txt")
samples1
#ooh! fix the sample IDs at this point.
samples1$Sample <- gsub("H179_1", "H051_8", samples1$Sample)
samples1$Sample <- gsub("H179_2", "H051_9", samples1$Sample)
samples1$Sample <- gsub("H179_3", "H051_10", samples1$Sample)
samples1$Sample <- gsub("H179_4", "H051_11", samples1$Sample)
samples1$Sample <- gsub("H179_5", "H051_12", samples1$Sample)
samples1$Sample <- gsub("H179_6", "H051_13", samples1$Sample)

#& H120 is considered the 2nd set for H182 (see metadata)
samples1$Sample <- gsub("H120_1", "H182_8", samples1$Sample)
samples1$Sample <- gsub("H120_2", "H182_9", samples1$Sample)
samples1$Sample <- gsub("H120_3", "H182_10", samples1$Sample)
samples1$Sample <- gsub("H120_4", "H182_11", samples1$Sample)
samples1$Sample <- gsub("H120_5", "H182_12", samples1$Sample)
samples1$Sample <- gsub("H120_6", "H182_13", samples1$Sample)

samples1_meta <- samples1 %>% rename("SampleID"="Sample") %>%
  full_join(meta_seq)
samples2 <- samples1_meta %>% filter(Status.=="completed") %>%
  select(SampleID, Read1, Read2, Notes) %>%
  rename("Sample"="SampleID")
samples2 #quick greps suggest that mocks, blanks, pilot, etc. are all gone.

write.table(samples2, "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/samples2.txt",
            sep="\t", quote=FALSE)

####Re-Import####
#meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_FILT_22.11.21.csv") #old version not circadian rhythm
#meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_FILT_23.01.20.csv") #new version circadian rhythm 48 & exclude infected & old samples
brk_cls <- fread("/ebio/abt3_projects2/uHEAT/data/brk_cls_FILT_23.01.20.csv") #exclude low depth, blank, mock, infected & old samples
tax <- read.csv("/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken/tax_table.csv",row.names=1)
#meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.03.14.csv") #starting back from here.
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.04.01.csv", row.names=1) #fever7d

#incorporate the new ELISA and TLR5 data, Sept 25 2023
elisa <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/ELISA_TLR5_results.csv", row.names=1)
meta_seq <- full_join(meta_seq, elisa)
meta_seq %>% filter(is.na(Status.)) #so here is one issue: all those extra relative-day data are missing Participant ID information

#define immune response from beginning. ...
meta_seq$Immune_response <- meta_seq$IgG_B2 / meta_seq$fever7d
meta_seq$Immune_response_aboveavg <- ifelse(meta_seq$Immune_response>median(meta_seq$Immune_response, na.rm=TRUE), "hi", "lo") #mean or median??
dplyr::count(meta_seq, Immune_response_aboveavg)

#and make a clean meta version, all the info but per person.
meta <- meta_seq %>% filter(!duplicated(Participant_ID))

meta_seq %>% filter(!duplicated(Participant_ID)) %>% dim()
meta_seq %>% filter(!duplicated(Participant_ID)) %>% count(fever_aboveavg) %>% as.data.frame() %>% mutate(per=n/sum(n))

#above: re-read "exclude", and exclude them.


#additionally add: CRP and IgG data to the post-cleaned meta_seq version. (to avoid changing this file anymore and causing confusion)
#crp <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/CRP.csv", row.names=1)
#igg <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/IgG_titres.csv", row.names=1)
#and: sample o.d. as a further technical parameter
#od <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/SampleID-OD_final_forR.csv")

#meta_seq <- Reduce(full_join, list(meta_seq, crp, igg, od))
#meta_seq <- meta_seq %>% filter(!is.na(Status.))

#double check the missing relative days here and re-define them. super annoying.
#I think actually those samples are missing, like, everything? not just relative days?
meta_seq %>% filter(is.na(Relative_Day)) %>% select(SampleID, Relative_Day, Age) #now it's only some that more or less make sense
#missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
#missing_id <- meta_seq %>% filter(is.na(Relative_Day)) %>% select(SampleID, Relative_Day, Time0)
#missing_id$SampleID %in% missing_days$SampleID #yes, these are the exact same samples. So there was a joining error somewhere.
#question is, should I fix it before, or now. Ach I can't concentrate today.
#with the new start again it's only some of the samples missing which is a little weird.
#temp <- missing_days %>% select(SampleID, Relative_Day, Time0)
#temp$Participant_ID <- sapply(strsplit(temp$SampleID,"_"), `[`, 1)
#temp$Participant_ID <- gsub("H", "HEAT_", temp$Participant_ID)
#temp <- temp %>% full_join(meta_seq)
#temp %>% filter(is.na(Participant_ID)|is.na(Relative_Day))
#so. maybe what I want to do is: first, get rid the column in meta seq. ? then, join?

#extra fever data: to compare with the drugs.
extra_fever <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/fevers_comparewindows_2023.03.17.csv")
meta_seq_extra <- extra_fever[,c(1, 8:17)] %>% full_join(meta_seq)

#an additional variable for the antibody analysis is: time between vaccine and B2
meta_seq$Date_StudyVacc_Days #already in days
#use the same formula I used before to convert visit 2 date to days
date_to_days <- function(x) {years <- strsplit(x, "-")[[1]][1]
years <- ifelse(years=="2021", yes=0, no=1) #assuming no=2022 
months <- strsplit(x, "-")[[1]][2]
days <- strsplit(x, "-")[[1]][3]
days2 <- as.numeric(years)*365 + as.numeric(months)*30 + as.numeric(days)
return(days2)

}
meta_seq$Visit_2_Date_Days <- sapply(meta_seq$Visit_2_Date, date_to_days)
meta_seq$Visit_1_Date_Days <- sapply(meta_seq$Visit_1_Date, date_to_days)
meta_seq$Visit_3_Date_Days <- sapply(meta_seq$Visit_3_Date, date_to_days)

meta_seq$B2_time_since_last_vaccine <- meta_seq$Visit_2_Date_Days - meta_seq$Date_StudyVacc_Days
summary(meta_seq$B2_time_since_last_vaccine) #so some of this is definitely wrong
meta_seq %>% filter(B2_time_since_last_vaccine<0) %>% select(Participant_ID, Visit_1_Date, Visit_2_Date, Visit_3_Date, Date_StudyVacc)
#okay it's a clear data copy error, visit 2 should be January 2022 not 2021
meta_seq[which(meta_seq$Participant_ID=="HEAT_052"),]$Visit_2_Date <- c("2022-01-03")
meta_seq$Visit_2_Date_Days <- sapply(meta_seq$Visit_2_Date, date_to_days)
meta_seq$B2_time_since_last_vaccine <- meta_seq$Visit_2_Date_Days - meta_seq$Date_StudyVacc_Days
summary(meta_seq$B2_time_since_last_vaccine) #the zero still seems weird. did someone accidentally copy the same date twice?
meta_seq %>% filter(B2_time_since_last_vaccine==0) %>% select(Participant_ID, Visit_1_Date, Visit_2_Date, Visit_3_Date, Date_StudyVacc)
#but for most people it's almost exactly 2 weeks (14-18 days), that's nice.
hist(meta_seq$B2_time_since_last_vaccine)

#export the updated meta_seq file.
#paste to avoid date irritations when I'mcleaning in excel.
#meta_seq$Bristol_Stool_Score2 <- paste("b", meta_seq$Bristol_Stool_Score2)
#meta_seq$Alcohol_category_per_week <- paste(meta_seq$Alcohol_category_per_week, "drinks", sep="_")
#meta_seq$BowelMovements_per_day  #garbage but unimportant


write.csv(meta_seq, "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.04.01.csv")
#reimport



#make a cleaner version: no antipyretics, no infection
meta_seq_clean1 <- meta_seq %>% 
  filter(!is.na(fever48) & vaccine_antipyretics %in% c(NA, "no") & COVID_infection_during_study=="") 
meta_seq_clean1 %>% filter(!duplicated(Participant_ID)) %>% dim() #n=140 people remaining
meta_seq_clean1 %>% filter(!duplicated(Participant_ID)) %>% count(fever_aboveavg) %>% as.data.frame() %>% mutate(per=n/sum(n))

#make a cleaner+ version: and they all got a BioNTech booster
meta_seq_clean2 <- meta_seq_clean1 %>% 
  filter(StudyVaccine_Type=="BioNTech" & Vaccine_StudyDose!="1")
meta_seq_clean2 %>% filter(!duplicated(Participant_ID)) %>% dim() #n=125 people remaining
meta_seq_clean2 %>% filter(!duplicated(Participant_ID)) %>% count(fever_aboveavg) %>% as.data.frame() %>% mutate(per=n/sum(n))

#make a cleanest version: only BioNTech 3rd dose for all-BioNTech people
meta_seq_clean3 <- meta_seq_clean1 %>% 
  filter(Vaccines_all=="all_BioNTech" & Vaccine_StudyDose=="3")
meta_seq_clean3 %>% filter(!duplicated(Participant_ID)) %>% dim() #n=72 people remaining
meta_seq_clean3 %>% filter(!duplicated(Participant_ID)) %>% count(fever_aboveavg) %>% as.data.frame() %>% mutate(per=n/sum(n))

####metadata clean 3 Aug 2024 - no doubles, keep series 2####
#re-processing a metaseq file including a fever7d max and excluding double participant series.
fever7d <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/24.08.13.Fever7d.Recleaned.csv") #the actual new one
#per-sample metadata #note: mocks, blanks & samples with low read depth already excluded #so are old & infected samples
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.04.01.csv", 
                     row.names=1)

#per-person metadata, matching
meta <- meta_seq %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta)

#do the new and old fever7ds match??
fever7d <- fever7d %>% filter(Participant_ID %in% meta$Participant_ID)
#fever7d <- fever7d[order(fever7d$Participant_ID),]
meta <- meta[order(meta$Participant_ID),]
identical(fever7d$Participant_ID, meta$Participant_ID)
identical(fever7d$fever7d, meta$fever7d) #omg thank god

#and now let's for the doubles (the ones with 11-12 samples)
count <- as.data.frame(dplyr::count(meta_seq, Participant_ID))
doubles <- count %>% filter(n>=11)
doubles <- meta_seq %>% filter(Participant_ID %in% doubles$Participant_ID & !is.na(Participant_ID))
dplyr::count(doubles, Participant_ID) %>% as.data.frame()

meta_seq %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179", "HEAT_051", "HEAT_120", "HEAT_182")) %>%
  select(SampleID, Vaccine_StudyDose, Vaccine1_Type, Vaccine2_Type, Vaccine3_Type, StudyVaccine_Type) %>%
  as.data.frame() %>% View()

#174 is also the (only) strange one that got Novovax for their study dose?
dplyr::count(meta, StudyVaccine_Type) %>% as.data.frame() #eben
meta %>% filter(Participant_ID=="HEAT_174") %>% as.data.frame() #and low fever
#I think we should exclude HEAT174 entirely from clean analysis.

#so: HEAT 114 was already filtered, 182 also (instead of 120)

#for the rest: take the series #2 samples, and rename study dose as dose 2 instead of dose 1
dplyr::count(doubles, Participant_ID) %>% as.data.frame()
doubles_to_exclude <- doubles %>% filter(Participant_ID %in% c(doubles$Participant_ID)&
                                           sample_no<7) #exclude the first series

doubles_to_exclude$SampleID
doubles_to_exclude2 <- doubles %>% filter(Participant_ID=="HEAT_174") #exclude the whole series
doubles_to_exclude2$SampleID

#write.csv(doubles_to_exclude, "/ebio/abt3_projects2/uHEAT/data/metadata/doubles_to_exclude2024.08.csv")

#add the fever information back from the new df
overlap_names <- names(fever7d)[which(names(fever7d)%in%names(meta_seq))]
overlap_names <- overlap_names[2:length(overlap_names)] #but keep Participant_ID, #1

meta_seq_nodoubles <- meta_seq %>% 
  dplyr::select(-c(overlap_names)) %>% #exclude all the temperature data that already overlap
  full_join(fever7d, by=c("Participant_ID"))
#and now exclude the double samples, i.e. 'series 1'
meta_seq_nodoubles <- meta_seq_nodoubles %>% 
  filter(!(SampleID %in% c(doubles_to_exclude$SampleID, doubles_to_exclude2$SampleID)))
#and for those doubles, update the vaccine dose and the recorded fever, as well as baseline / max.
meta_seq_nodoubles$Vaccine_StudyDose <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                               yes=meta_seq_nodoubles$Vaccine_StudyDose+1,
                                               no=meta_seq_nodoubles$Vaccine_StudyDose)
names(fever7d)
meta_seq_nodoubles$fever7d <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                     yes=meta_seq_nodoubles$fever7d2,
                                     no=meta_seq_nodoubles$fever7d)
meta_seq_nodoubles$max7d <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                   yes=meta_seq_nodoubles$max7d2,
                                   no=meta_seq_nodoubles$max7d)
meta_seq_nodoubles$max_base <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                      yes=meta_seq_nodoubles$max_base2,
                                      no=meta_seq_nodoubles$max_base)
meta_seq_nodoubles$median_base <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                         yes=meta_seq_nodoubles$median_base2,
                                         no=meta_seq_nodoubles$median_base)
meta_seq_nodoubles$median_after <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                          yes=meta_seq_nodoubles$median_after2,
                                          no=meta_seq_nodoubles$median_after)
meta_seq_nodoubles$Time_of_Day_Fever <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                               yes=meta_seq_nodoubles$Time_of_Day_Fever2,
                                               no=meta_seq_nodoubles$Time_of_Day_Fever)

meta_seq_nodoubles %>% filter(Participant_ID %in% doubles$Participant_ID) %>% 
  dplyr::select(SampleID, Vaccine_StudyDose, fever7d, fever7d2) %>% as.data.frame()

#for H179 vs H051: rename H051 as H179, because right now only H051 is there.
meta_seq_nodoubles[which(meta_seq_nodoubles$Participant_ID=="HEAT_179"),]$Participant_ID
meta_seq_nodoubles[which(meta_seq_nodoubles$Participant_ID=="HEAT_051"),]$Participant_ID <- c("HEAT_179")

#correct the sample ids and sample numbers for doubles - but only for H179 - then I never have to do it again.
meta_seq_nodoubles$SampleID <- gsub("H051_8", "H179_1", meta_seq_nodoubles$SampleID)
meta_seq_nodoubles$SampleID <- gsub("H051_9", "H179_2", meta_seq_nodoubles$SampleID)
meta_seq_nodoubles$SampleID <- gsub("H051_10", "H179_3", meta_seq_nodoubles$SampleID)
meta_seq_nodoubles$SampleID <- gsub("H051_11", "H179_4",  meta_seq_nodoubles$SampleID)
meta_seq_nodoubles$SampleID <- gsub("H051_12", "H179_5", meta_seq_nodoubles$SampleID)
meta_seq_nodoubles$SampleID <- gsub("H051_13", "H179_6", meta_seq_nodoubles$SampleID)
#make sure I go through and fix all the markdowns, after...
#do NOT rename the others, sequencing data will not match then lol
#and for H182 vs  H120: H120 is already excluded everywhere

#the antibody titres must also be corrected for the doubles!! note also for somalogic, meta.
#hm.. some of the antibody titres are missing, however, e.g. from H179
#bah. re-import blood files too. 
crp <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/CRP.csv", row.names=1)
igg <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/IgG_titres.csv", row.names=1)
igg %>% filter(Participant_ID %in% c(doubles$Participant_ID,"HEAT_179")) %>% as.data.frame()
#OK fine, it's really just H179 that I'm missing then.

meta_seq_nodoubles$IgG_B1 <- ifelse(meta_seq_nodoubles$Participant_ID %in% doubles$Participant_ID & !is.na(meta_seq_nodoubles$IgG_B3),
                                    yes=meta_seq_nodoubles$IgG_B2,
                                    meta_seq_nodoubles$IgG_B1)
meta_seq_nodoubles$IgG_B2 <- ifelse(meta_seq_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_nodoubles$IgG_B3),
                                    yes=meta_seq_nodoubles$IgG_B3,
                                    meta_seq_nodoubles$IgG_B2)
meta_seq_nodoubles$CRP_B1 <- ifelse(meta_seq_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_nodoubles$CRP_B3),
                                    yes=meta_seq_nodoubles$CRP_B2,
                                    meta_seq_nodoubles$CRP_B1)
meta_seq_nodoubles$CRP_B2 <- ifelse(meta_seq_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_nodoubles$CRP_B3),
                                    yes=meta_seq_nodoubles$CRP_B3,
                                    meta_seq_nodoubles$CRP_B2)

#aaaand special case 179
meta_seq_nodoubles$IgG_B1 <- ifelse(meta_seq_nodoubles$Participant_ID=="HEAT_179",
                                    yes=igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B1,
                                    meta_seq_nodoubles$IgG_B1)
meta_seq_nodoubles$IgG_B2 <- ifelse(meta_seq_nodoubles$Participant_ID=="HEAT_179",
                                    yes=igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B2,
                                    meta_seq_nodoubles$IgG_B2)
meta_seq_nodoubles$CRP_B1 <- ifelse(meta_seq_nodoubles$Participant_ID=="HEAT_179",
                                    yes=crp[which(crp$Participant_ID=="HEAT_179"),]$CRP_B1,
                                    meta_seq_nodoubles$CRP_B1)
meta_seq_nodoubles$CRP_B2 <- ifelse(meta_seq_nodoubles$Participant_ID=="HEAT_179",
                                    yes=crp[which(crp$Participant_ID=="HEAT_179"),]$CRP_B2,
                                    meta_seq_nodoubles$CRP_B2)

meta_seq_nodoubles %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2, CRP_B1, CRP_B2) %>% as.data.frame()

#so this meta_seq file should now be correct by doubles, and also has a max raw temp per person.

#does categorization change at all if I use the updated no-doubles, undup data to define fever7d above avg?
meta_nodoubles <- meta_seq_nodoubles %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta_nodoubles)
meta_nodoubles$fever7d_aboveavg_NEW <- ifelse(meta_nodoubles$fever7d>mean(meta_nodoubles$fever7d, na.rm=TRUE),
                                              "hi", "lo")
identical(meta_nodoubles$fever7d_aboveavg, meta_nodoubles$fever7d_aboveavg_NEW)
meta_nodoubles[which(meta_nodoubles$fever7d_aboveavg!=meta_nodoubles$fever7d_aboveavg_NEW),] %>%
  dplyr::select(Participant_ID, fever7d, fever7d_aboveavg, fever7d_aboveavg_NEW) %>%
  as.data.frame()
#the only difference is HEAT_033, which is reasonable, because now I'm taking fever2 (hi) instead of 1 (lo)

#so: best would be indeed to redefine the fever7d category cleanly here.
#add the new category back to meta.seq
new_fever <- meta_nodoubles %>% dplyr::select(Participant_ID, fever7d_aboveavg_NEW) %>%
  rename(fever7d_aboveavg=fever7d_aboveavg_NEW)
as.data.frame(new_fever)
meta_seq_nodoubles_newfever <- meta_seq_nodoubles %>% 
  dplyr::select(-c(fever7d_aboveavg)) %>% 
  full_join(new_fever, by=c("Participant_ID"))

dplyr::count(meta_seq_nodoubles_newfever, Participant_ID) %>% as.data.frame()
meta_seq_nodoubles_newfever %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Time01, Stool_timing3, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2) %>% as.data.frame()


#but there have been some problem samples retained with no fever.
meta_seq_nodoubles_newfever %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame() %>% View()

#I think it makes sense.

#and now, export it.
write.csv(meta_seq_nodoubles_newfever, 
          "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_24.08.13.csv")

####metadata clean 15 Aug 2024 - fill in missing NA####
#as above except adding back in missing samples. boo. 

#re-processing a metaseq file including a fever7d max and excluding double participant series.
fever7d <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/24.08.13.Fever7d.Recleaned.csv") #the actual new one
#per-sample metadata #note: mocks, blanks & samples with low read depth already excluded #so are old & infected samples
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.04.01.csv", 
                     row.names=1)

#some specific samples are mysteriously missing metadata, still. 
missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
still_missing_days <- meta_seq %>% filter(is.na(fever7d_aboveavg))
still_missing_days$SampleID %in% missing_days$SampleID #OMG yes, they are mostly there! fuck's sake man I thought I fixed it.
#I guess the question is whether I really want those data but I think generally I do, mostly they are just 'out of range' samples
#the real problem children (i.e. infections etc) were already deliberately fully excluded
missing_days$Time01 <- missing_days$Time0 - 7.5 #assume they got their impfung in the exact middle. 
#so, add additional sample data for merging.
missing_days$sample_no <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 2)
missing_days$sample_no <- as.numeric(missing_days$sample_no)
missing_days$Participant_ID <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 1)
missing_days$Participant_ID <-  gsub("H", "HEAT_", missing_days$Participant_ID)
#make sure these samples aren't the ones that were deliberately excluded e.g. for low seq depth or infection
excluded <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_EXCLUDED.csv")
excluded %>% dplyr::select(SampleID, Reason)
missing_days$SampleID %in% excluded$SampleID
excluded %>% filter(SampleID %in% missing_days$SampleID) %>% as.data.frame() %>% dplyr::select(SampleID, Reason)
#only one of them should be deliberately excluded, the rest were just accidentally missing for typos, timing etc.
missing_days <- missing_days %>% filter(!(SampleID %in% excluded$SampleID))
#add sequencing depth information for the missing samples.
seqstats1 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats2 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC_pool3/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats <- full_join(seqstats1, seqstats2) %>% 
  filter(Read==1) %>%
  rename(SampleID=Sample, seq_depth=num_seqs) %>% 
  dplyr::select(SampleID, seq_depth) %>%
  mutate(seq_depth=seq_depth*2)

missing_days <- missing_days %>% full_join(seqstats) %>% filter(!is.na(Participant_ID))
missing_days

#bristol stool score, time at RT
#add bss + technical metadata (plate, pool, time at room temp) #actually this also contains the bss so it's fine
batch <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/batch-and-technical-metadata.csv")
batch <- batch %>% rename(SampleID=sample_id, Date_Sampled=Date, Time_Sampled=Time) %>% 
  dplyr::select(-c(sample_no, sample_type, X, X.1, Participant_ID)) %>%
  rename(Comment1=Comment)
batch %>% filter(SampleID %in% missing_days$SampleID)
missing_days <- missing_days %>% full_join(batch) %>%
  filter(!is.na(Participant_ID))
missing_days

#now add the rest of the metadata.
names_exclude <- names(missing_days)
names_exclude <- names_exclude[!names_exclude%in%c("Participant_ID", "Date", "Time", "Comment1")]
missing_days_meta <- meta_seq %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  as.data.frame() %>%
  dplyr::select(-c(names_exclude)) %>%
  full_join(missing_days) %>%
  filter(!is.na(Comment1))
View(missing_days_meta)
missing_days_meta %>% dplyr::select(SampleID, Relative_Day, Comment, Participant_ID,
                                    seq_depth, Bristol_Stool_Score, 
                                    Age, Sex, fever7d_aboveavg) %>% 
  as.data.frame() #okay, nice, this is relatively complete.
#and now... to somehow add it back to the normal meta_seq file...
meta_seq_full <- meta_seq %>%
  filter(!is.na(fever7d_aboveavg)) %>%
  full_join(missing_days_meta)
dim(meta_seq_full)
meta_seq_full %>% filter(SampleID %in% missing_days$SampleID) %>%
  as.data.frame()
head(meta_seq_full)
tail(meta_seq_full)
duplicated(meta_seq_full$SampleID)
meta_seq_full %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame()
#I think it's ok now. continue with 'meta_seq_full' as my input table to create the new fever one.

#per-person metadata, matching
meta <- meta_seq_full %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta)

#do the new and old fever7ds match??
fever7d <- fever7d %>% filter(Participant_ID %in% meta$Participant_ID)
#fever7d <- fever7d[order(fever7d$Participant_ID),]
meta <- meta[order(meta$Participant_ID),]
identical(fever7d$Participant_ID, meta$Participant_ID)
identical(fever7d$fever7d, meta$fever7d) #omg thank god

#and now let's for the doubles (the ones with 11-12 samples)
count <- as.data.frame(dplyr::count(meta_seq_full, Participant_ID))
doubles <- count %>% filter(n>=11)
doubles <- meta_seq_full %>% filter(Participant_ID %in% doubles$Participant_ID & !is.na(Participant_ID))
dplyr::count(doubles, Participant_ID) %>% as.data.frame()

meta_seq_full %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179", "HEAT_051", "HEAT_120", "HEAT_182")) %>%
  select(SampleID, Vaccine_StudyDose, Vaccine1_Type, Vaccine2_Type, Vaccine3_Type, StudyVaccine_Type) %>%
  as.data.frame() %>% View()

#174 is also the (only) strange one that got Novovax for their study dose?
dplyr::count(meta, StudyVaccine_Type) %>% as.data.frame() #eben
meta %>% filter(Participant_ID=="HEAT_174") %>% as.data.frame() #and low fever
#I think we should exclude HEAT174 entirely from clean analysis.

#so: HEAT 114 was already filtered, 182 also (instead of 120)

#for the rest: take the series #2 samples, and rename study dose as dose 2 instead of dose 1
dplyr::count(doubles, Participant_ID) %>% as.data.frame()
doubles_to_exclude <- doubles %>% filter(Participant_ID %in% c(doubles$Participant_ID)&
                                           sample_no<7) #exclude the first series

doubles_to_exclude$SampleID
doubles_to_exclude2 <- doubles %>% filter(Participant_ID=="HEAT_174") #exclude the whole series
doubles_to_exclude2$SampleID

#write.csv(doubles_to_exclude, "/ebio/abt3_projects2/uHEAT/data/metadata/doubles_to_exclude2024.08.csv")

#add the fever information back from the new df
overlap_names <- names(fever7d)[which(names(fever7d)%in%names(meta_seq_full))]
overlap_names <- overlap_names[2:length(overlap_names)] #but keep Participant_ID, #1

meta_seq_full_nodoubles <- meta_seq_full %>% 
  dplyr::select(-c(overlap_names)) %>% #exclude all the temperature data that already overlap
  full_join(fever7d, by=c("Participant_ID"))
#and now exclude the double samples, i.e. 'series 1'
meta_seq_full_nodoubles <- meta_seq_full_nodoubles %>% 
  filter(!(SampleID %in% c(doubles_to_exclude$SampleID, doubles_to_exclude2$SampleID)))
#and for those doubles, update the vaccine dose and the recorded fever, as well as baseline / max.
meta_seq_full_nodoubles$Vaccine_StudyDose <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                                    yes=meta_seq_full_nodoubles$Vaccine_StudyDose+1,
                                                    no=meta_seq_full_nodoubles$Vaccine_StudyDose)
names(fever7d)
meta_seq_full_nodoubles$fever7d <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                          yes=meta_seq_full_nodoubles$fever7d2,
                                          no=meta_seq_full_nodoubles$fever7d)
meta_seq_full_nodoubles$max7d <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                        yes=meta_seq_full_nodoubles$max7d2,
                                        no=meta_seq_full_nodoubles$max7d)
meta_seq_full_nodoubles$max_base <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                           yes=meta_seq_full_nodoubles$max_base2,
                                           no=meta_seq_full_nodoubles$max_base)
meta_seq_full_nodoubles$median_base <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                              yes=meta_seq_full_nodoubles$median_base2,
                                              no=meta_seq_full_nodoubles$median_base)
meta_seq_full_nodoubles$median_after <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                               yes=meta_seq_full_nodoubles$median_after2,
                                               no=meta_seq_full_nodoubles$median_after)
meta_seq_full_nodoubles$Time_of_Day_Fever <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID,
                                                    yes=meta_seq_full_nodoubles$Time_of_Day_Fever2,
                                                    no=meta_seq_full_nodoubles$Time_of_Day_Fever)

meta_seq_full_nodoubles %>% filter(Participant_ID %in% doubles$Participant_ID) %>% 
  dplyr::select(SampleID, Vaccine_StudyDose, fever7d, fever7d2) %>% as.data.frame()

#for H179 vs H051: rename H051 as H179, because right now only H051 is there.
meta_seq_full_nodoubles[which(meta_seq_full_nodoubles$Participant_ID=="HEAT_179"),]$Participant_ID
meta_seq_full_nodoubles[which(meta_seq_full_nodoubles$Participant_ID=="HEAT_051"),]$Participant_ID <- c("HEAT_179")

#correct the sample ids and sample numbers for doubles - but only for H179 - then I never have to do it again.
meta_seq_full_nodoubles$SampleID <- gsub("H051_8", "H179_1", meta_seq_full_nodoubles$SampleID)
meta_seq_full_nodoubles$SampleID <- gsub("H051_9", "H179_2", meta_seq_full_nodoubles$SampleID)
meta_seq_full_nodoubles$SampleID <- gsub("H051_10", "H179_3", meta_seq_full_nodoubles$SampleID)
meta_seq_full_nodoubles$SampleID <- gsub("H051_11", "H179_4",  meta_seq_full_nodoubles$SampleID)
meta_seq_full_nodoubles$SampleID <- gsub("H051_12", "H179_5", meta_seq_full_nodoubles$SampleID)
meta_seq_full_nodoubles$SampleID <- gsub("H051_13", "H179_6", meta_seq_full_nodoubles$SampleID)
#make sure I go through and fix all the markdowns, after...
#do NOT rename the others, sequencing data will not match then lol
#and for H182 vs  H120: H120 is already excluded everywhere
##BUT ADD BACK THE CORRECT SEQUENCING DEPTH FOR H179###
seqstats %>% filter(SampleID %in% c("H179_1", "H179_2", "H179_3", "H179_4", "H179_5", "H179_6"))
meta_seq_full_nodoubles %>% filter(Participant_ID=="HEAT_179") %>% 
  dplyr::select(SampleID, seq_depth) %>% 
  as.data.frame()
#Oh no it's okay, it's already correct

#the antibody titres must also be corrected for the doubles!! note also for somalogic, meta.
#hm.. some of the antibody titres are missing, however, e.g. from H179
#bah. re-import blood files too. 
crp <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/CRP.csv", row.names=1)
igg <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/IgG_titres.csv", row.names=1)
igg %>% filter(Participant_ID %in% c(doubles$Participant_ID,"HEAT_179")) %>% as.data.frame()
#OK fine, it's really just H179 that I'm missing then.

meta_seq_full_nodoubles$IgG_B1 <- ifelse(meta_seq_full_nodoubles$Participant_ID %in% doubles$Participant_ID & !is.na(meta_seq_full_nodoubles$IgG_B3),
                                         yes=meta_seq_full_nodoubles$IgG_B2,
                                         meta_seq_full_nodoubles$IgG_B1)
meta_seq_full_nodoubles$IgG_B2 <- ifelse(meta_seq_full_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_full_nodoubles$IgG_B3),
                                         yes=meta_seq_full_nodoubles$IgG_B3,
                                         meta_seq_full_nodoubles$IgG_B2)
meta_seq_full_nodoubles$CRP_B1 <- ifelse(meta_seq_full_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_full_nodoubles$CRP_B3),
                                         yes=meta_seq_full_nodoubles$CRP_B2,
                                         meta_seq_full_nodoubles$CRP_B1)
meta_seq_full_nodoubles$CRP_B2 <- ifelse(meta_seq_full_nodoubles$Participant_ID%in% doubles$Participant_ID & !is.na(meta_seq_full_nodoubles$CRP_B3),
                                         yes=meta_seq_full_nodoubles$CRP_B3,
                                         meta_seq_full_nodoubles$CRP_B2)

#aaaand special case 179
meta_seq_full_nodoubles$IgG_B1 <- ifelse(meta_seq_full_nodoubles$Participant_ID=="HEAT_179",
                                         yes=igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B1,
                                         meta_seq_full_nodoubles$IgG_B1)
meta_seq_full_nodoubles$IgG_B2 <- ifelse(meta_seq_full_nodoubles$Participant_ID=="HEAT_179",
                                         yes=igg[which(igg$Participant_ID=="HEAT_179"),]$IgG_B2,
                                         meta_seq_full_nodoubles$IgG_B2)
meta_seq_full_nodoubles$CRP_B1 <- ifelse(meta_seq_full_nodoubles$Participant_ID=="HEAT_179",
                                         yes=crp[which(crp$Participant_ID=="HEAT_179"),]$CRP_B1,
                                         meta_seq_full_nodoubles$CRP_B1)
meta_seq_full_nodoubles$CRP_B2 <- ifelse(meta_seq_full_nodoubles$Participant_ID=="HEAT_179",
                                         yes=crp[which(crp$Participant_ID=="HEAT_179"),]$CRP_B2,
                                         meta_seq_full_nodoubles$CRP_B2)

meta_seq_full_nodoubles %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2, CRP_B1, CRP_B2) %>% as.data.frame()

#so this meta_seq_full file should now be correct by doubles, and also has a max raw temp per person.

#does categorization change at all if I use the updated no-doubles, undup data to define fever7d above avg?
meta_nodoubles <- meta_seq_full_nodoubles %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta_nodoubles)
meta_nodoubles$fever7d_aboveavg_NEW <- ifelse(meta_nodoubles$fever7d>mean(meta_nodoubles$fever7d, na.rm=TRUE),
                                              "hi", "lo")
identical(meta_nodoubles$fever7d_aboveavg, meta_nodoubles$fever7d_aboveavg_NEW)
meta_nodoubles[which(meta_nodoubles$fever7d_aboveavg!=meta_nodoubles$fever7d_aboveavg_NEW),] %>%
  dplyr::select(Participant_ID, fever7d, fever7d_aboveavg, fever7d_aboveavg_NEW) %>%
  as.data.frame()
#the only difference is HEAT_033, which is reasonable, because now I'm taking fever2 (hi) instead of 1 (lo)

#so: best would be indeed to redefine the fever7d category cleanly here.
#add the new category back to meta.seq
new_fever <- meta_nodoubles %>% dplyr::select(Participant_ID, fever7d_aboveavg_NEW) %>%
  rename(fever7d_aboveavg=fever7d_aboveavg_NEW)
as.data.frame(new_fever)
meta_seq_full_nodoubles_newfever <- meta_seq_full_nodoubles %>% 
  dplyr::select(-c(fever7d_aboveavg)) %>% 
  full_join(new_fever, by=c("Participant_ID"))

dplyr::count(meta_seq_full_nodoubles_newfever, Participant_ID) %>% as.data.frame()
meta_seq_full_nodoubles_newfever %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Time01, Stool_timing3, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2) %>% as.data.frame()


#but there have been some problem samples retained with no fever.
meta_seq_full_nodoubles_newfever %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame() %>% View()
#nothing is NA anymore
dim(meta_seq_full_nodoubles_newfever)
#I think it makes sense.

#wait, was this also a problem for the meta_seq_dna?
dplyr::count(meta_seq_full_nodoubles_newfever, Participant_ID) #one person has 7 samples
duplicated(meta_seq_full_nodoubles_newfever$SampleID) #H133_3 got counted twice

meta_seq_full_nodoubles_newfever <- meta_seq_full_nodoubles_newfever %>% filter(!duplicated(SampleID))

#another remaining problem: missing accurate BSS data for the second series people
missing_bss <- meta_seq_full_nodoubles_newfever %>% filter(is.na(Bristol_Stool_Score))
missing_bss
#do I easily have that data?
batch %>% filter(SampleID %in% missing_bss$SampleID) %>% View()
#no...not really... only for H179
#what about in the other bss file?
bss <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/bss_temp_long_for_individual_plots.csv")
bss <- bss %>% filter(Entry=="Stool")
bss <- bss %>% filter(Participant_ID %in% missing_bss$Participant_ID)
View(bss) #but it's wrong because of the doubles ID... and still confusing... just do this manually
#write out the missing table so I know what to look for.
fill_me <- batch %>% filter(SampleID %in% missing_bss$SampleID) %>% as.data.frame()
write.csv(fill_me, "/ebio/abt3_projects2/uHEAT/data/metadata/missing_bss_doubles.csv")
#reimport the filled data
bss_fill <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_bss_doubles_manual2024.08.15.csv")
bss_fill
#add the time at rt via date to days
date_to_days <- function(x) {years <- strsplit(x, "-")[[1]][1]
years <- ifelse(years=="21", yes=0, no=1) #assuming no=2022 
months <- strsplit(x, "-")[[1]][2]
days <- strsplit(x, "-")[[1]][3]
days2 <- as.numeric(years)*365 + as.numeric(months)*30 + as.numeric(days)
return(days2)
}
bss_fill$Date_Sampled <- gsub(".", "-", bss_fill$Date_Sampled, fixed=TRUE)
bss_fill$Date_Sampled_Days <- sapply(bss_fill$Date_Sampled, date_to_days)
bss_fill
#ok. nowwww merge it back. would have made sense to add it earlier but whatever.
bss_fill$CommentBSS <- c("Manually_added_BSS")
names_exclude <- names(bss_fill)
names_exclude <- names_exclude[!names_exclude%in%c("SampleID", "CommentBSS")]
#is there other PER-SAMPLE metadata that is NOT missing for those samples?
names_exclude2 <- names(missing_days)
per_sample_bss <- missing_bss %>% dplyr::select(c(names_exclude2)) %>% as.data.frame() #yes, there is plenty that exists.
#okay, first merge by adding the per-sample metadata.

bss_fill2 <- per_sample_bss %>% dplyr::select(-c(Date, Time, names_exclude)) %>%
  full_join(bss_fill)
dim(bss_fill)
dim(bss_fill2)
View(bss_fill2)
bss_fill2$Time_at_RT <- bss_fill2$Date_Frozen_Days - bss_fill2$Date_Sampled_Days

#OK, great.
names_exclude3 <- names(bss_fill2)
names_exclude3 <- names_exclude3[!names_exclude3%in%c("Participant_ID", "CommentBSS")]
names_exclude3
#also exclude "Stool Timing", because this is incorrect after merge and I think the safest is to redefine later based on Time01
names_exclude3 <- c(names_exclude3, "Stool_timing1", "Stool_timing2", "Stool_timing3",
                    "X", "X.1", "X.2")

#and now merge by ParticipantID to obtain the rest of the metadata
bss_fill2$Participant_ID
#this should yield complete metadata, ONLY for those samples..
missing_bss_meta <- meta_seq_full_nodoubles_newfever %>% 
  filter(!duplicated(Participant_ID)) %>%
  as.data.frame() %>%
  dplyr::select(-c(names_exclude3)) %>%
  full_join(bss_fill2) %>% 
  filter(SampleID %in% bss_fill2$SampleID)
View(missing_bss_meta) #now it's complete
dim(missing_bss_meta) #just the missing samples
names(missing_bss_meta)

#and now... to somehow add it back to the normal meta_seq file...
meta_seq_full_nodoubles_newfever$Bristol_Stool_Score <- as.numeric(meta_seq_full_nodoubles_newfever$Bristol_Stool_Score)
missing_bss_meta$CommentBSS <- NULL

meta_seq_full_nodoubles_newfever_bss <- meta_seq_full_nodoubles_newfever %>%
  filter(!is.na(fever7d_aboveavg)) %>%
  filter(!(SampleID %in% missing_bss_meta$SampleID)) %>%
  full_join(missing_bss_meta)
dim(meta_seq_full_nodoubles_newfever_bss) #1021 samples
meta_seq_full_nodoubles_newfever_bss %>% filter(SampleID %in% missing_bss$SampleID) %>%
  as.data.frame() %>% View()
head(meta_seq_full_nodoubles_newfever_bss)
tail(meta_seq_full_nodoubles_newfever_bss)
duplicated(meta_seq_full_nodoubles_newfever_bss$SampleID)
dplyr::count(meta_seq_full_nodoubles_newfever_bss, Participant_ID) #now it is not doubled.
meta_seq_full_nodoubles_newfever_bss %>% filter(Participant_ID=="HEAT_179") %>% View()
#right so the problem for example, is differences in recorded Times / technical metadata.
meta_seq_full_nodoubles_newfever_bss %>%
  dplyr::select(names_exclude3) %>%
  View()
meta_seq_full_nodoubles_newfever_bss$Time_at_RT
meta_seq_full_nodoubles_newfever_bss$Bristol_Stool_Score

#God. OK and now, export it.
write.csv(meta_seq_full_nodoubles_newfever, 
          "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_24.08.15.csv")
####metadata clean 18 Aug 2024 - no doubles, keep series *1* & fill missing NA ####
#as above except filter to keep series 1 instead of 2

#re-processing a metaseq file including a fever7d max and excluding double participant series.
fever7d <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/24.08.13.Fever7d.Recleaned.csv") #the actual new one
#per-sample metadata #note: mocks, blanks & samples with low read depth already excluded #so are old & infected samples
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_23.04.01.csv", 
                     row.names=1)

#some specific samples are mysteriously missing metadata, still. 
missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
still_missing_days <- meta_seq %>% filter(is.na(fever7d_aboveavg))
still_missing_days$SampleID %in% missing_days$SampleID #OMG yes, they are mostly there! fuck's sake man I thought I fixed it.
#I guess the question is whether I really want those data but I think generally I do, mostly they are just 'out of range' samples
#the real problem children (i.e. infections etc) were already deliberately fully excluded
missing_days$Time01 <- missing_days$Time0 - 7.5 #assume they got their impfung in the exact middle. 
#so, add additional sample data for merging.
missing_days$sample_no <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 2)
missing_days$sample_no <- as.numeric(missing_days$sample_no)
missing_days$Participant_ID <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 1)
missing_days$Participant_ID <-  gsub("H", "HEAT_", missing_days$Participant_ID)
#make sure these samples aren't the ones that were deliberately excluded e.g. for low seq depth or infection
excluded <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_EXCLUDED.csv")
excluded %>% dplyr::select(SampleID, Reason)
missing_days$SampleID %in% excluded$SampleID
excluded %>% filter(SampleID %in% missing_days$SampleID) %>% as.data.frame() %>% dplyr::select(SampleID, Reason)
#only one of them should be deliberately excluded, the rest were just accidentally missing for typos, timing etc.
missing_days <- missing_days %>% filter(!(SampleID %in% excluded$SampleID))
#add a categorization for before vs after
missing_days$Relative_Day_numeric <- gsub("t+", "", missing_days$Relative_Day, fixed=TRUE)
missing_days$Relative_Day_numeric <- gsub("t", "", missing_days$Relative_Day_numeric, fixed=TRUE)
missing_days$Relative_Day_numeric <- as.numeric(missing_days$Relative_Day_numeric)
missing_days$Stool_timing3 <- ifelse(missing_days$Relative_Day_numeric<=0, yes="before", no="after")
missing_days

#add sequencing depth information for the missing samples.
seqstats1 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats2 <- read.table("/ebio/abt3_projects2/uHEAT/data/output_QC_pool3/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats <- full_join(seqstats1, seqstats2) %>% 
  filter(Read==1) %>%
  rename(SampleID=Sample, seq_depth=num_seqs) %>% 
  dplyr::select(SampleID, seq_depth) %>%
  mutate(seq_depth=seq_depth*2)

missing_days <- missing_days %>% full_join(seqstats) %>% filter(!is.na(Participant_ID))
missing_days

#bristol stool score, time at RT
#add bss + technical metadata (plate, pool, time at room temp) #actually this also contains the bss so it's fine
batch <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/batch-and-technical-metadata.csv")
batch <- batch %>% rename(SampleID=sample_id, Date_Sampled=Date, Time_Sampled=Time) %>% 
  dplyr::select(-c(sample_no, sample_type, X, X.1, Participant_ID)) %>%
  rename(Comment1=Comment)
batch %>% filter(SampleID %in% missing_days$SampleID)
missing_days <- missing_days %>% full_join(batch) %>%
  filter(!is.na(Participant_ID))
missing_days

#now add the rest of the metadata.
names_exclude <- names(missing_days)
names_exclude <- names_exclude[!names_exclude%in%c("Participant_ID", "Date", "Time", "Comment1")]
missing_days_meta <- meta_seq %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  as.data.frame() %>%
  dplyr::select(-c(names_exclude)) %>%
  full_join(missing_days) %>%
  filter(!is.na(Comment1))
View(missing_days_meta)
missing_days_meta %>% dplyr::select(SampleID, Relative_Day, Comment, Participant_ID,
                                    seq_depth, Bristol_Stool_Score, 
                                    Age, Sex, fever7d_aboveavg) %>% 
  as.data.frame() #okay, nice, this is relatively complete.
#and now... to somehow add it back to the normal meta_seq file...
meta_seq_full <- meta_seq %>%
  filter(!is.na(fever7d_aboveavg)) %>%
  full_join(missing_days_meta)
dim(meta_seq_full)
meta_seq_full %>% filter(SampleID %in% missing_days$SampleID) %>%
  as.data.frame()
head(meta_seq_full)
tail(meta_seq_full)
duplicated(meta_seq_full$SampleID)
meta_seq_full %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame()
#I think it's ok now. continue with 'meta_seq_full' as my input table to create the new fever one.

#per-person metadata, matching
meta <- meta_seq_full %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta)

#do the new and old fever7ds match??
fever7d <- fever7d %>% filter(Participant_ID %in% meta$Participant_ID)
#fever7d <- fever7d[order(fever7d$Participant_ID),]
meta <- meta[order(meta$Participant_ID),]
identical(fever7d$Participant_ID, meta$Participant_ID)
identical(fever7d$fever7d, meta$fever7d) #omg thank god

#and now let's for the doubles (the ones with 11-12 samples)
count <- as.data.frame(dplyr::count(meta_seq_full, Participant_ID))
doubles <- count %>% filter(n>=11)
doubles <- meta_seq_full %>% filter(Participant_ID %in% doubles$Participant_ID & !is.na(Participant_ID))
dplyr::count(doubles, Participant_ID) %>% as.data.frame() #plus 179/051

#meta_seq_full %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179", "HEAT_051", "HEAT_120", "HEAT_182")) %>%
#  select(SampleID, Vaccine_StudyDose, Vaccine1_Type, Vaccine2_Type, Vaccine3_Type, StudyVaccine_Type) %>%
#  as.data.frame() %>% View()

#174 is also the (only) strange one that got Novovax for their study dose?
dplyr::count(meta, StudyVaccine_Type) %>% as.data.frame() #eben
meta %>% filter(Participant_ID=="HEAT_174") %>% as.data.frame() #and low fever
#I think we should exclude HEAT174 entirely from clean analysis.

#so: HEAT 114 was already filtered, 182 also (instead of 120)

#for the rest: take the series #2 samples, and rename study dose as dose 2 instead of dose 1
dplyr::count(doubles, Participant_ID) %>% as.data.frame()
doubles_to_exclude <- doubles %>% filter(Participant_ID %in% c(doubles$Participant_ID)&
                                           sample_no>7) #exclude the SECOND series 
#I'm very undecided for H051 vs H179 but I think either is fine, there's no particular reason either way
#I guess for simplicity also just keep the first available series, i.e. H051, as before

doubles_to_exclude$SampleID
doubles_to_exclude2 <- doubles %>% filter(Participant_ID %in% c("HEAT_174")) #exclude the whole series
doubles_to_exclude2$SampleID
doubles_to_exclude_exp <- meta_seq %>% filter(SampleID %in% c(doubles_to_exclude$SampleID, doubles_to_exclude2$SampleID))

#write.csv(doubles_to_exclude_exp, "/ebio/abt3_projects2/uHEAT/data/metadata/doubles_to_exclude2024.08.18.csv")

#add the fever information back from the new df
overlap_names <- names(fever7d)[which(names(fever7d)%in%names(meta_seq_full))]
overlap_names <- overlap_names[2:length(overlap_names)] #but keep Participant_ID, #1

meta_seq_full_nodoubles <- meta_seq_full %>% 
  dplyr::select(-c(overlap_names)) %>% #exclude all the temperature data that already overlap
  full_join(fever7d, by=c("Participant_ID"))
#and now exclude the double samples, i.e. 'series 1'
meta_seq_full_nodoubles <- meta_seq_full_nodoubles %>% 
  filter(!(SampleID %in% c(doubles_to_exclude$SampleID, doubles_to_exclude2$SampleID)))
#unlike in previous version: now all the fever info should be correct since it was already there for series 1
meta_seq_full_nodoubles %>% filter(Participant_ID %in% doubles$Participant_ID) %>% 
  dplyr::select(SampleID, Vaccine_StudyDose, fever7d, fever7d2) %>% as.data.frame()

#for H179 vs H051: this time we are keeping H051 so it's okay.
#all the H051_8 or H179_1 samples are now just excluded from meta_seq so it's okay, don't need to rename later either
#and for H182 vs  H120: H120 is already excluded everywhere

meta_seq_full_nodoubles %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2, CRP_B1, CRP_B2) %>% as.data.frame()

#so this meta_seq_full file should now be correct by doubles, and also has a max raw temp per person.
#okay looks reasonable

#does categorization change at all if I use the updated no-doubles, undup data to define fever7d above avg?
meta_nodoubles <- meta_seq_full_nodoubles %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta_nodoubles)
meta_nodoubles$fever7d_aboveavg_NEW <- ifelse(meta_nodoubles$fever7d>mean(meta_nodoubles$fever7d, na.rm=TRUE),
                                              "hi", "lo")
identical(meta_nodoubles$fever7d_aboveavg, meta_nodoubles$fever7d_aboveavg_NEW) #now it's the same as before
#because e.g. for H033 I went back to taking series 1, low fever.

#so: no need to redefine the fever7d category #wait except why do I have 179 instead of 051
meta_nodoubles$fever7d_aboveavg_NEW <- NULL

dplyr::count(meta_seq_full_nodoubles, Participant_ID) %>% as.data.frame()
meta_seq_full_nodoubles %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Time01, Stool_timing3, Vaccine_StudyDose, fever7d, fever7d2,
                IgG_B1, IgG_B2) %>% as.data.frame()


#but there have been some problem samples retained with no fever.
meta_seq_full_nodoubles %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame() %>% View()
#nothing is NA anymore
dim(meta_seq_full_nodoubles)
#I think it makes sense.

#wait, was this also a problem for the meta_seq_dna?
dplyr::count(meta_seq_full_nodoubles, Participant_ID) #one person has 7 samples
duplicated(meta_seq_full_nodoubles$SampleID) 
meta_seq_full_nodoubles %>% filter(duplicated(SampleID)) %>% as.data.frame() #H133_3 and H040_3 got counted twice
sample_dups <- meta_seq_full_nodoubles %>% filter(duplicated(SampleID))
meta_seq_full_nodoubles %>% filter(SampleID %in% sample_dups$SampleID) %>% as.data.frame() %>% View()
#for H040_3, keep the one where I imputed stool timing as 'before', even though Time01 is missing
meta_seq_full_nodoubles <- meta_seq_full_nodoubles %>% filter(!(SampleID=="H040_3"&is.na(Stool_timing3)))
#for H133_3, simply exclude one, they are meaningfully identical
meta_seq_full_nodoubles <- meta_seq_full_nodoubles %>% filter(!duplicated(SampleID))
meta_seq_full_nodoubles %>% filter(duplicated(SampleID)) %>% as.data.frame() #now it is clean.

#double check, is BSS then also still there?
meta_seq_full_nodoubles %>% filter(Participant_ID %in% c(doubles$Participant_ID, "HEAT_179")) %>% 
  dplyr::select(Participant_ID, SampleID, Time01, Stool_timing3, seq_depth, Bristol_Stool_Score) %>% as.data.frame()
meta_seq_full_nodoubles[order(meta_seq_full_nodoubles$Participant_ID),] %>% filter(Participant_ID %in% missing_days$Participant_ID) %>% 
  dplyr::select(Participant_ID, SampleID, Time01, Stool_timing3, seq_depth, Bristol_Stool_Score) %>% as.data.frame()
#yes.

#I should have just kept series 1 from the beginning lol, so much easier. 
dplyr::count(meta_seq_full_nodoubles, Participant_ID)
#OK and now, export it.
write.csv(meta_seq_full_nodoubles, 
          "/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_series1.24.08.18.csv")
