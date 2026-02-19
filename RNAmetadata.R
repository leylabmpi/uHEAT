####Process a new clean *RNA* metadata - keeping series 2 ####
#the original RNA metadata
meta_seq_rna <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d.csv") 
dim(meta_seq_rna)

#the new DNA metadata
meta_seq_dna <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_24.08.15.csv", 
                     row.names=1)

#the RNA sequencing depth 
seqstats1 <- read.table("/ebio/abt3_projects2/uHEAT2/data/output_QC_RNAseq/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats2 <- read.table("/ebio/abt3_projects2/uHEAT2/data/output_QC_RNAseq_reps23/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats <- full_join(seqstats1, seqstats2) %>% 
  filter(Read==1) %>%
  rename(SampleID=Sample, seq_depth=num_seqs) %>% 
  mutate(RNA_seq_depth=seq_depth*2) %>%
  dplyr::select(SampleID, RNA_seq_depth)

#are any of the missing_days people participants in RNAseq?
missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
missing_days$Participant_ID <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 1)
missing_days$Participant_ID <-  gsub("H", "HEAT_", missing_days$Participant_ID)

missing_days$Participant_ID %in% meta_seq_rna$Participant_ID #yes
missing_days[which(missing_days$Participant_ID %in% meta_seq_rna$Participant_ID),]
#lol okay crazy.
meta_seq_rna[which(meta_seq_rna$SampleID %in% missing_days$SampleID),] %>% View()
meta_seq_rna[which(meta_seq_rna$SampleID %in% missing_days$SampleID),] %>% 
  dplyr::select(SampleID, sample_no, Relative_Day, Time_at_RT) %>%
  as.data.frame()
#but for those samples the data has not been missing, anyways, for whatever reason. Sigh. But good.

#before I take that approach, I need to re-correct the H051 to H179
meta_seq_rna$SampleID <- gsub("H051_8", "H179_1", meta_seq_rna$SampleID)
meta_seq_rna$SampleID <- gsub("H051_9", "H179_2", meta_seq_rna$SampleID)
meta_seq_rna$SampleID <- gsub("H051_10", "H179_3", meta_seq_rna$SampleID)
meta_seq_rna$SampleID <- gsub("H051_11", "H179_4",  meta_seq_rna$SampleID)
meta_seq_rna$SampleID <- gsub("H051_12", "H179_5", meta_seq_rna$SampleID)
meta_seq_rna$SampleID <- gsub("H051_13", "H179_6", meta_seq_rna$SampleID)

#so (now?) I can just filter the new DNA table and add correct RNA_seq_depth
meta_seq_rna_new <- meta_seq_dna %>% filter(SampleID %in% meta_seq_rna$SampleID)
dim(meta_seq_rna) #n=242
dim(meta_seq_rna_new) #smaller (n=218) because some double series are now excluded.
#add the rna seq depth - but actually add it from the old table - because I forgot about naming issues
rna_seq_depth <- meta_seq_rna %>% dplyr::select(SampleID, RNA_seq_depth) %>%
  filter(SampleID %in% meta_seq_dna$SampleID)


meta_seq_rna_new <- meta_seq_rna_new %>%
  rename(DNA_seq_depth=seq_depth) %>%
  full_join(rna_seq_depth)
dim(meta_seq_rna_new) #n=223

#and H179 is there, yes?
meta_seq_rna_new %>% filter(Participant_ID=="HEAT_179") #yes, now it is

dplyr::count(meta_seq_rna_new, Participant_ID) #looks good
meta_seq_rna_new %>% filter(Participant_ID=="HEAT_133") %>% View()

#so now it's fine?
dplyr::count(meta_seq_rna_new, Participant_ID) 
meta_seq_rna_new %>% 
  dplyr::select(SampleID, Time01, Bristol_Stool_Score, RNA_seq_depth, DNA_seq_depth, fever7d)

#now it looks quite nice
meta_seq_rna_new %>% filter(is.na(fever7d_aboveavg)) #none, good
meta_rna_new <- meta_seq_rna_new %>% filter(!duplicated(Participant_ID))
dplyr::count(meta_rna_new, fever7d_aboveavg)

#okay woo!! write it out.
write.csv(meta_seq_rna_new, "/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d_24.08.15.csv")

####posthoc check####
#checking why I have one fewer participant than before, who is it?
meta_seq_rna_new <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d_24.08.15.csv") 
meta_seq_rna <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d.csv") 
dplyr::count(meta_seq_rna, Participant_ID)
dplyr::count(meta_seq_rna_new, Participant_ID)

setdiff(meta_seq_rna$Participant_ID, meta_seq_rna_new$Participant_ID)
#H051 replaced with H179, that's fine.
#but H074 got kicked out entirely??
meta_seq_rna %>% filter(Participant_ID=="HEAT_074") %>% dplyr::select(SampleID) %>% as.data.frame()
#ah... because I only actually DID RNAseq on the first series of H074...

####Process a new clean *RNA* metadata - keeping series 1 ####
#the original RNA metadata
meta_seq_rna <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d.csv") 
dim(meta_seq_rna)

#the new DNA metadata
#updated: this version also excludes SERIES 2 of all double series.
meta_seq_dna <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_series1.24.08.18.csv", 
                     row.names=1)

#the RNA sequencing depth 
seqstats1 <- read.table("/ebio/abt3_projects2/uHEAT2/data/output_QC_RNAseq/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats2 <- read.table("/ebio/abt3_projects2/uHEAT2/data/output_QC_RNAseq_reps23/reports/final/seqkit_stats.tsv", header=TRUE) 
seqstats <- full_join(seqstats1, seqstats2) %>% 
  filter(Read==1) %>%
  rename(SampleID=Sample, seq_depth=num_seqs) %>% 
  mutate(RNA_seq_depth=seq_depth*2) %>%
  dplyr::select(SampleID, RNA_seq_depth)

#are any of the missing_days people participants in RNAseq?
missing_days <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/missing_relative_days_manual_22.11.05.csv")
missing_days$Participant_ID <-  sapply(strsplit(missing_days$SampleID,"_"), `[`, 1)
missing_days$Participant_ID <-  gsub("H", "HEAT_", missing_days$Participant_ID)

missing_days$Participant_ID %in% meta_seq_rna$Participant_ID #yes
missing_days[which(missing_days$Participant_ID %in% meta_seq_rna$Participant_ID),]
#lol okay crazy.
meta_seq_rna[which(meta_seq_rna$SampleID %in% missing_days$SampleID),] %>% View()
meta_seq_rna[which(meta_seq_rna$SampleID %in% missing_days$SampleID),] %>% 
  dplyr::select(SampleID, sample_no, Relative_Day, Time_at_RT) %>%
  as.data.frame()
#but for those samples the data has not been missing, anyways, for whatever reason. Sigh. But good.

#so (now?) I can just filter the new DNA table and add correct RNA_seq_depth
meta_seq_rna_new <- meta_seq_dna %>% filter(SampleID %in% meta_seq_rna$SampleID)
dim(meta_seq_rna) #n=242
dim(meta_seq_rna_new) #smaller (n=230) because some double series are now excluded, but now fewer
#add the rna seq depth - but actually add it from the old table - because I forgot about naming issues
rna_seq_depth <- meta_seq_rna %>% dplyr::select(SampleID, RNA_seq_depth) %>%
  filter(SampleID %in% meta_seq_dna$SampleID)

meta_seq_rna_new <- meta_seq_rna_new %>%
  rename(DNA_seq_depth=seq_depth) %>%
  full_join(rna_seq_depth)
dim(meta_seq_rna_new) #n=230

#and H179 is NOT there
meta_seq_rna_new %>% filter(Participant_ID=="HEAT_179") 
#instead we have H051
dplyr::count(meta_seq_rna_new, Participant_ID) %>% as.data.frame() #looks good

#so now it's fine?
meta_seq_rna_new %>% 
  dplyr::select(SampleID, Time01, Bristol_Stool_Score, RNA_seq_depth, DNA_seq_depth, fever7d) %>%
  as.data.frame()

#now it looks quite nice
meta_seq_rna_new %>% filter(is.na(fever7d_aboveavg)) %>% as.data.frame() #none, good
meta_rna_new <- meta_seq_rna_new %>% filter(!duplicated(Participant_ID))
dplyr::count(meta_rna_new, fever7d_aboveavg) %>% as.data.frame() #this version is also better balanced

#okay woo!! write it out.
write.csv(meta_seq_rna_new, "/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d_series1_24.08.18.csv")
