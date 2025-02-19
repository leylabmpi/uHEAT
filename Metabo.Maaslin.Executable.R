#Metabolomics code executable
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(Maaslin2)
####chunk 1 import####
#import metadata
#updated: this version also excludes SERIES 2 of all double series.
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_series1.24.08.18.csv", 
                     row.names=1)

#define DIET 3
meta_seq$Diet3 <- ifelse(meta_seq$Diet %in% c("Vegan", "Vegetarian"), yes="Veg", no="Omni")

#import metabolomics
pos <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metabolomics/MS019_pos_final_Nov2024.csv") #updated with higher confidence metabolites
neg <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metabolomics/MS019_neg_final_Nov2024.csv") #updated with higher confidence metabolites

#per-person metadata, matching
meta <- meta_seq %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta)

#retrieve annotations directly from the data file
names(pos)
dim(pos)
annotations_pos <- pos[,c(1:32)]
names(neg)
dim(neg)
annotations_neg <- neg[,c(1:32)]

#extract the data itself, separately., as I had it before
pos <- pos[,c(1:3, 33:length(names(pos)))]
neg <- neg[,c(1:3, 33:length(names(neg)))]

#make a unified annotation column - pick the first name if it is exists, otherwise the second
annotations_pos[annotations_pos==""] <- NA
annotations_pos$annotation <- ifelse(!is.na(annotations_pos$spectral_db_matches.compound_name),
                                     yes=annotations_pos$spectral_db_matches.compound_name,
                                     no=annotations_pos$name)
#and then if both of those names are missing, use the class.
annotations_pos$annotation <- ifelse(is.na(annotations_pos$annotation),
                                     yes=annotations_pos$ClassyFire.most.specific.class,
                                     no=annotations_pos$annotation)
annotations_pos$annotation
#make a small version for merging
ap <- annotations_pos %>% dplyr::select(id, rt, mz, annotation, molecularFormula, 
                                        NPC.pathway, NPC.pathway.Probability, 
                                        NPC.superclass, NPC.superclass.Probability)
ap #and now it is clear that some of them are dups, so only do the annotation merging aftper Maaslin2.
#make a non-numeric id for rownames and column names.
ap$id <- paste("Mp", ap$id, sep="_")

#make a unified annotation column - pick the first name if it is exists, otherwise the second
annotations_neg[annotations_neg==""] <- NA
annotations_neg$annotation <- ifelse(!is.na(annotations_neg$spectral_db_matches.compound_name),
                                     yes=annotations_neg$spectral_db_matches.compound_name,
                                     no=annotations_neg$name)
#and then if both of those names are missing, use the class.
annotations_neg$annotation <- ifelse(is.na(annotations_neg$annotation),
                                     yes=annotations_neg$ClassyFire.most.specific.class,
                                     no=annotations_neg$annotation)
annotations_neg$annotation
#make a small version for merging
an <- annotations_neg %>% dplyr::select(id, rt, mz, annotation, molecularFormula, 
                                        NPC.pathway, NPC.pathway.Probability, 
                                        NPC.superclass, NPC.superclass.Probability)
an #and now it is clear that some of them are dups, so only do the annotation merging after Maaslin2.
#make a non-numeric id for rownames and column names.
an$id <- paste("Mn", an$id, sep="_")

#cleaning it up cont. #weirdly in the new positive version there are 'duplicate ids'.
row.names(pos) <- paste("Mp", pos$id, sep="_")
ftp <- pos %>% dplyr::select(-c("id", "rt", "mz")) #start with positive ion mode - name based on Petras lab script
ftp <- ftp[,order(colnames(ftp)), drop=F]

row.names(neg) <- paste("Mn", neg$id, sep="_")
ftn <- neg %>% dplyr::select(-c("id", "rt", "mz")) #start with positive ion mode - name based on Petras lab script
ftn <- ftn[,order(colnames(ftn)), drop=F]

#make a 'metabolomics' metadata file which contains all the sample info also for blanks and QCs
#can I do this identically for the neg and pos ones? #technically not but the only diff is pos/neg.
identical(names(ftp), names(ftn)) #OK still not identical even with pos/neg designation removed. just make two separate metadata files for now.
mdp <- data.frame(Metabo_SampleID=c(names(ftp)))
mdp$Participant_ID=sapply(strsplit(as.character(mdp$Metabo_SampleID),"_"), `[`, 2)
mdp$Participant_ID <- gsub("H", "HEAT_", mdp$Participant_ID)
mdp$Blood <- sapply(strsplit(as.character(mdp$Metabo_SampleID),"_"), `[`, 3)
mdp$Blood <- ifelse(mdp$Participant_ID %in% c("Pos", "Neg"), yes=NA, no=mdp$Blood)
mdp$Participant_ID <- ifelse(mdp$Participant_ID%in%c("Pos", "Neg"), yes=NA, no=mdp$Participant_ID)
mdp$SampleType <- sapply(strsplit(as.character(mdp$Metabo_SampleID),"_"), `[`, 1)
mdp$SampleType <- ifelse(mdp$SampleType %in% c("Blank", "QC"), yes=mdp$SampleType, no="Serum")
#meta_small <- meta %>% dplyr::select(Participant_ID, fever7d_aboveavg) #start with just minimal metadata to keep things simple
mdp <- mdp %>% full_join(meta) %>% filter(!is.na(SampleType)) #keep full metadata for all the samples
head(mdp)
row.names(mdp) <- mdp$Metabo_SampleID
mdp <- mdp[order(row.names(mdp)),,drop=F]
#did any samples still get run twice? I think it is fixed now.
mdp$Metabo_SampleNo=sapply(strsplit(as.character(mdp$Metabo_SampleID),"_"), `[`, 1)
mdp[which(duplicated(mdp$Metabo_SampleNo)),]$Metabo_SampleNo #yes, #71 and #97 are duplicated?
mdp[which(duplicated(mdp$Metabo_SampleNo)),]$Metabo_SampleID
dups <- mdp %>% filter(Metabo_SampleNo %in% c("71", "97"))
ftp[,which(names(ftp) %in% dups$Metabo_SampleID)]
#remove the duplicates from the ftp table according to Dennis' instructions on which one is a better intensity.
mdp <- mdp %>% filter(!(Metabo_SampleID %in% c("71_H065_B2_Pos", "97_H107_B1_Pos")))
mdp[which(duplicated(mdp$Metabo_SampleNo)),]$Metabo_SampleNo #Now it's clean.
ftp <- ftp[,which(names(ftp)%in%mdp$Metabo_SampleID)]

# how many files in the metadata are also present in the feature table
table(rownames(mdp) %in% colnames(ftp))
identical(rownames(mdp), colnames(ftp))
#transpose the feature table and ensure numeric
ftp_t <- as.data.frame(t(ftp)) #transposing the ftp
ftp_t <- ftp_t %>% mutate_all(as.numeric)  #converting all values to numeric
print("identical(rownames(mdp),rownames(ftp_t))") #for log file
identical(rownames(mdp),rownames(ftp_t)) #should return TRUE now

##repeat for negative mode.
mdn <- data.frame(Metabo_SampleID=c(names(ftn)))
mdn$Participant_ID=sapply(strsplit(as.character(mdn$Metabo_SampleID),"_"), `[`, 2)
mdn$Participant_ID <- gsub("H", "HEAT_", mdn$Participant_ID)
mdn$Blood <- sapply(strsplit(as.character(mdn$Metabo_SampleID),"_"), `[`, 3)
mdn$Blood <- ifelse(mdn$Participant_ID %in% c("Pos", "Neg"), yes=NA, no=mdn$Blood)
mdn$Participant_ID <- ifelse(mdn$Participant_ID%in%c("Pos", "Neg"), yes=NA, no=mdn$Participant_ID)
mdn$SampleType <- sapply(strsplit(as.character(mdn$Metabo_SampleID),"_"), `[`, 1)
mdn$SampleType <- ifelse(mdn$SampleType %in% c("Blank", "QC"), yes=mdn$SampleType, no="Serum")
#meta_small <- meta %>% dplyr::select(Participant_ID, fever7d_aboveavg) #start with just minimal metadata to keep things simple
mdn <- mdn %>% full_join(meta) %>% filter(!is.na(SampleType)) #keep full metadata for all the samples
head(mdn)
row.names(mdn) <- mdn$Metabo_SampleID
mdn <- mdn[order(row.names(mdn)),,drop=F]
#did any samples still get run twice? I think it is fixed now.
mdn$Metabo_SampleNo=sapply(strsplit(as.character(mdn$Metabo_SampleID),"_"), `[`, 1)
mdn[which(duplicated(mdn$Metabo_SampleNo)),]$Metabo_SampleNo #yes, #149 and 93 and 95 dups
mdn[which(duplicated(mdn$Metabo_SampleNo)),]$Metabo_SampleID #yes, #149 and 93 and 95 dups
dups <- mdn %>% filter(Metabo_SampleNo %in% c("93", "95", "149"))
ftn[,which(names(ftn) %in% dups$Metabo_SampleID)]
#remove the duplicates from the ftp table according to Dennis' instructions on which one is a better intensity.
#later noticed 2 additional technical duplicates, filtered based on which one had fewer NA
mdn <- mdn %>% filter(!(Metabo_SampleID %in% c("149_H038_B2_Neg_1", "93_H154_B1_Neg_1", "93_H154_B1_Neg_2", "95_H155_B1_Neg_1",
                                               "X114_H051_B1_Neg", "X154_H104_B1_Neg_1")))
mdn[which(duplicated(mdn$Metabo_SampleNo)),]$Metabo_SampleNo #now it's clean.
ftn <- ftn[,which(names(ftn)%in%mdn$Metabo_SampleID)]

# how many files in the metadata are also present in the feature table
table(rownames(mdn) %in% colnames(ftn))
identical(rownames(mdn), colnames(ftn))
#transpose the feature table and ensure numeric
ftn_t <- as.data.frame(t(ftn)) #transposing the ftn
ftn_t <- ftn_t %>% mutate_all(as.numeric)  #converting all values to numeric
print("identical(rownames(mdn),rownames(ftn_t))") #for log file
identical(rownames(mdn),rownames(ftn_t)) #should return TRUE now

####chunk2 filter####
#first, I think that all NAs should be considered zeroes for the purposes of means, prevalence, etc.
ftp_t[is.na(ftp_t)] <- 0
ftn_t[is.na(ftn_t)] <- 0

#POS.
#separate samples and blanks.
blank_ids_p <- mdp[which(mdp$SampleType=="Blank"),]$Metabo_SampleID
blank_ids_p
qc_ids_p <- mdp[which(mdp$SampleType=="QC"),]$Metabo_SampleID
qc_ids_p
sample_ids_p <- mdp[which(mdp$SampleType=="Serum"),]$Metabo_SampleID
sample_ids_p

# Getting the corresponding rows from ftp_t
Blank_p <- ftp_t[which(rownames(ftp_t) %in% (blank_ids_p)), , drop=F]
Samples_p <- ftp_t[which(rownames(ftp_t) %in% (sample_ids_p)), , drop=F]

#Define a cutoff for blank features removal.
#When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
Cutoff <- 0.3 # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3

#Perform Blank Removal
#Getting mean for every feature in blank and Samples in a data frame named 'Avg_ftp'
#note - this doesn't help for prevalence. Later, we shoudl also filter for prevalence.
Avg_ftp <- data.frame(Avg_blank=colMeans(Blank_p, na.rm = T)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ftp$`Avg_samples` <- colMeans(Samples_p, na.rm= T) # adding another column 'Avg_samples' for feature means of samples
#this still contains many NaN from when the ENTIRE COLUMN is blank.
#all NaNs now replaced by zeroes to begin with

#Getting the ratio of blank vs Sample
Avg_ftp$`Ratio_blank_Sample` <- (Avg_ftp$`Avg_blank`+1)/(Avg_ftp$`Avg_samples`+1)

# Creating a bin with 1s when the ratio>Cutoff (EXCLUDE), else put 0s (INCLUDE)
Avg_ftp$`Bg_bin` <- ifelse(Avg_ftp$`Ratio_blank_Sample` > Cutoff, 1, 0 )
Avg_ftp #there are a lot of NA!!! from things that appear in samples but not at all in blanks.
Blank_p$Mp_3 #for example, this whole thing is entirely NA. Or zero.
#Avg_ftp$Bg_bin <- ifelse(is.na(Avg_ftp$Avg_blank), yes=0, no=Avg_ftp$Bg_bin) #correct the NA blanks to 'do not remove' zeros


#Calculating the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ftp)))
print(paste("No.of Background or noise features:",sum(Avg_ftp$`Bg_bin` ==1,na.rm = T))) 
print(paste("No.of features after excluding noise:",(ncol(Samples_p) - sum(Avg_ftp$`Bg_bin` ==1,na.rm = T))))

#do this in two steps to avoid lo
identical(names(Samples_p), row.names(Avg_ftp)) #the metabolite ids are the same
blk_rem_p <- merge(as.data.frame(t(Samples_p)), Avg_ftp, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  dplyr::select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) #removing the last 4 columns
row.names(blk_rem_p) <- blk_rem_p$Row.names #extra step bc columns-to-rownames didn't work
blk_rem_p <- as.data.frame(t(blk_rem_p %>% dplyr::select(-c(Row.names))))
dim(blk_rem_p) #now the samples are retained
names(blk_rem_p) #and so are the metabolite names :)

#NEG
#repeat for neg
#separate samples and blanks.
blank_ids_n <- mdn[which(mdn$SampleType=="Blank"),]$Metabo_SampleID
blank_ids_n
qc_ids_n <- mdn[which(mdn$SampleType=="QC"),]$Metabo_SampleID
qc_ids_n
sample_ids_n <- mdn[which(mdn$SampleType=="Serum"),]$Metabo_SampleID
sample_ids_n

# Getting the corresponding rows from ftp_t
Blank_n <- ftn_t[which(rownames(ftn_t) %in% (blank_ids_n)), , drop=F]
Samples_n <- ftn_t[which(rownames(ftn_t) %in% (sample_ids_n)), , drop=F]

#Define a cutoff for blank features removal.
#When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
Cutoff <- 0.3 # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3

#Perform Blank Removal
#Getting mean for every feature in blank and Samples in a data frame named 'Avg_ftn'
Avg_ftn <- data.frame(Avg_blank=colMeans(Blank_n, na.rm= T)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ftn$`Avg_samples` <- colMeans(Samples_n, na.rm= T) # adding another column 'Avg_samples' for feature means of samples


#Getting the ratio of blank vs Sample
Avg_ftn$`Ratio_blank_Sample` <- (Avg_ftn$`Avg_blank`+1)/(Avg_ftn$`Avg_samples`+1)

# Creating a bin with 1s when the ratio>Cutoff, else put 0s
Avg_ftn$`Bg_bin` <- ifelse(Avg_ftn$`Ratio_blank_Sample` > Cutoff, 1, 0 )
Avg_ftn #there are a lot of NA!!! from things that appear in samples but not at all in blanks.
#Blank$V1 #for example, this whole thing is entirely NA.
#Avg_ftn$Bg_bin <- ifelse(is.na(Avg_ftn$Avg_blank), yes=0, no=Avg_ftn$Bg_bin) #correct the NA blanks to 'do not remove' zeros
#Avg_ftn

#Calculating the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ftn)))
print(paste("No.of Background or noise features:",sum(Avg_ftn$`Bg_bin` ==1,na.rm = T))) 
print(paste("No.of features after excluding noise:",(ncol(Samples_n) - sum(Avg_ftn$`Bg_bin` ==1,na.rm = T))))

blk_rem_n <- merge(as.data.frame(t(Samples_n)), Avg_ftn, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  dplyr::select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) #removing the last 4 columns
row.names(blk_rem_n) <- blk_rem_n$Row.names #extra step bc columns-to-rownames didn't work
blk_rem_n <- as.data.frame(t(blk_rem_n %>% dplyr::select(-c(Row.names))))
dim(blk_rem_n) #now the samples are retained
names(blk_rem_n) #and so are the metabolite names :)

####chunk 3 maaslin2 no covariates####
#repeat on the BLNK REMOVED DATA - skip imputation, this was a mess - baseline only
mdp_clean <- mdp %>% filter(mdp$Metabo_SampleID %in% row.names(blk_rem_p)
                            & Blood=="B1") #baseline samples only!
blk_rem_p <- blk_rem_p[which(row.names(blk_rem_p)%in%row.names(mdp_clean)),]
identical(row.names(mdp_clean), row.names(blk_rem_p))
blk_rem_pt <- as.data.frame(t(blk_rem_p))
print("identical(row.names(mdp_clean), names(blk_rem_pt))") #for log file
identical(row.names(mdp_clean), names(blk_rem_pt))
row.names(blk_rem_pt) #and the metabolite ids are still retained.

mdp_clean$fever7d_aboveavg <- factor(mdp_clean$fever7d_aboveavg, levels=c("lo", "hi"))


#fever7d (pos data) - no covariates
fit_data = Maaslin2(
  input_data = blk_rem_pt, 
  input_metadata = mdp_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("fever7d_aboveavg"),
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
pos_fever <- fit_data$results
pos_fever %>% filter(qval<0.1) #centered scaled data also no association, this seems fair.
pos_fever %>% filter(pval<0.05) 
write.csv(pos_fever, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboPos_Maaslin2_Fever7d.csv")

#ALSO WRITE OUT the final clean files used for Maaslin2, so that I can plot them later.
write.csv(mdp_clean, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/MetaboPos_MdpClean.csv")
write.csv(blk_rem_pt, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/MetaboPos_blk_rem_pt.csv")

##NEG
mdn_clean <- mdn %>% filter(mdn$Metabo_SampleID %in% row.names(blk_rem_n)
                            &Blood=="B1") #baseline samples only!
blk_rem_n <- blk_rem_n[which(row.names(blk_rem_n)%in%row.names(mdn_clean)),]
identical(row.names(mdn_clean), row.names(blk_rem_n))
blk_rem_nt <- as.data.frame(t(blk_rem_n))
print("identical(row.names(mdn_clean), names(blk_rem_nt))")
identical(row.names(mdn_clean), names(blk_rem_nt))
row.names(blk_rem_nt) #and the metabolite ids are still retained.

mdn_clean$fever7d_aboveavg <- factor(mdn_clean$fever7d_aboveavg, levels=c("lo", "hi"))

#fever7d (neg data) - no covariates
fit_data = Maaslin2(
  input_data = blk_rem_nt, 
  input_metadata = mdn_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("fever7d_aboveavg"),
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
neg_fever <- fit_data$results #no hits
neg_fever %>% filter(qval<0.1) #no hits
neg_fever %>% filter(pval<0.05) #nothing even remotely impressive, actually
write.csv(neg_fever, "/ebio/abt3_projects2/uHEAT/code/R_code/Rscripts_server/maaslin2_output/MetaboNeg_Maaslin2_Fever7d.csv")

#ALSO WRITE OUT the final clean files used for Maaslin2, so that I can plot them later.
write.csv(mdn_clean, "/ebio/abt3_projects2/uHEAT/code/R_code/Rscripts_server/maaslin2_input/MetaboNeg_MdnClean.csv")
write.csv(blk_rem_nt, "/ebio/abt3_projects2/uHEAT/code/R_code/Rscripts_server/maaslin2_input/MetaboNeg_blk_rem_nt.csv")


####chunk 4 maaslin2 with covariates Sex Age####
identical(row.names(mdp_clean), names(blk_rem_pt))
#fever7d (pos data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_pt, 
  input_metadata = mdp_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "fever7d_aboveavg"),
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
pos_fever_corr1 <- fit_data$results
pos_fever_corr1 %>% filter(qval<0.05)
pos_fever_corr1 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(pos_fever_corr1, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboPos_Maaslin2_SexAgeFever7d.csv")

identical(row.names(mdn_clean), names(blk_rem_nt))
#fever7d (neg data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_nt, 
  input_metadata = mdn_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "fever7d_aboveavg"),
  reference=c("Sex", "m"), 
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
neg_fever_corr1 <- fit_data$results #no hits
neg_fever_corr1 %>% filter(qval<0.05) #no hits
neg_fever_corr1 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(neg_fever_corr1, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboNeg_Maaslin2_SexAgeFever7d.csv")


####chunk 4 maaslin2 with covariates Sex Age Diet####
identical(row.names(mdp_clean), names(blk_rem_pt))
#fever7d (pos data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_pt, 
  input_metadata = mdp_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "Diet3", "fever7d_aboveavg"),
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
pos_fever_corr2 <- fit_data$results
pos_fever_corr2 %>% filter(qval<0.05)
pos_fever_corr2 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(pos_fever_corr2, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboPos_Maaslin2_SexAgeDiet3Fever7d.csv")

identical(row.names(mdn_clean), names(blk_rem_nt))
#fever7d (neg data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_nt, 
  input_metadata = mdn_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "Diet3", "fever7d_aboveavg"),
  reference=c("Sex", "m"), 
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
neg_fever_corr2 <- fit_data$results #no hits
neg_fever_corr2 %>% filter(qval<0.05) #no hits
neg_fever_corr2 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(neg_fever_corr2, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboNeg_Maaslin2_SexAgeDiet3Fever7d.csv")


####chunk 4 maaslin2 with covariates Sex Age Diet BMI####
identical(row.names(mdp_clean), names(blk_rem_pt))
#fever7d (pos data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_pt, 
  input_metadata = mdp_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "Diet3", "BMI_kg_m2", "fever7d_aboveavg"),
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
pos_fever_corr3 <- fit_data$results
pos_fever_corr3 %>% filter(qval<0.05)
pos_fever_corr3 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(pos_fever_corr3, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboPos_Maaslin2_SexAgeDiet3BMIFever7d.csv")

identical(row.names(mdn_clean), names(blk_rem_nt))
#fever7d (neg data) - with covariates
fit_data = Maaslin2(
  input_data = blk_rem_nt, 
  input_metadata = mdn_clean, 
  output = "maaslin2_output", 
  fixed_effects = c("Sex", "Age", "Diet3", "BMI_kg_m2", "fever7d_aboveavg"),
  reference=c("Sex", "m"), 
  normalization="tss", #returned to default recommendation
  transform="log", #& log transform
  #min_abundance=0.0001, #no reasonable min abund exists
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, 
  plot_scatter=FALSE
)
neg_fever_corr3 <- fit_data$results #no hits
neg_fever_corr3 %>% filter(qval<0.05) #no hits
neg_fever_corr3 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 

write.csv(neg_fever_corr3, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/MetaboNeg_Maaslin2_SexAgeDiet3BMIFever7d.csv")

