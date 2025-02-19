#Transcriptomics Code Executable
#libraries
library(tidyr)
library(data.table)
library(dplyr)
library(Maaslin2)
library(ggplot2)
source("/ebio/abt3_projects2/uHEAT/code/R_code/LeyLabRMisc/R/init.R")

####function read bracken####
#' Function for reading in a bracken taxonomy table
#'
#' The table will be converted to long form (sample ~ abundance).
#' Only "_frac" or "_num" columns will be kept (see "keep_frac").
#' Taxonomy will be split into separate levels (see "tax_levs").
#' tidytable (w/ data.table) used to speed the process up.
#'
#' @param infile Path to bracken table file
#' @param nrows Number of table rows to read. If Inf, all lines will be read.
#' @param keep_frac If TRUE, keep all columns ending in "_frac"; otherwise, keep "_num" columns.
#' @param tax_levs Taxonomic levels to separate the taxonomy column into.
#' @param ... Params passed to fread()
#' @return data.table
#' @export
#' @import tidytable
#' @importFrom data.table fread
#' 
read_bracken = function(infile, nrows=Inf, keep_frac=TRUE,
                        tax_levs = c('Domain', 'Phylum', 'Class', 'Order',
                                     'Family', 'Genus', 'Species'),
                        nThread = 4, ...){
  if(keep_frac){
    to_rm = '_num'
    to_keep = '_frac'
  } else {
    to_rm = '_frac'
    to_keep = '_num'
  }
  dt = data.table::fread(infile, sep='\t', nrows=nrows, check.names=TRUE,
                         nThread=nThread, ...) %>%
    tidytable::select(-taxIDs, -ends_with(!!to_rm)) %>%
    tidytable::mutate(taxonomy = gsub(';[pcofgs]__', ';', taxonomy),
                      taxonomy = gsub('^d__', '', taxonomy)) %>%
    tidytable::separate(taxonomy, tax_levs, sep=';') %>%
    tidytable::pivot_longer(cols=ends_with(!!to_keep),
                            names_to='Sample',
                            values_to='Abundance') %>%
    tidytable::mutate(Sample = gsub('(_frac|_num)$', '', Sample))
  
  return(dt)
}

####chunk 3 import####
#read in the brk_cls file
##Read Files - rel abund
#output of llmgp #update: instead, for re-analysis, read in the ids only or the filtered table only
profile_dir <- "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken" #updated: no subsampling, GTDB map, all samples
# listing files
brk_cls_files = c("/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken/all-combined-bracken.tsv")
brk_cls_files %>% length()
# reading tables
brk_cls = brk_cls_files %>%
  plyr::llply(read_bracken) %>%
  data.table::rbindlist(use.names=TRUE, idcol='dataset')
brk_cls
names(brk_cls)

#re-parse participant IDs for later merging
brk_cls$Participant_ID <- sapply(strsplit(brk_cls$Sample,"_"), `[`, 1)
brk_cls$Participant_ID <- gsub("H", "HEAT_", brk_cls$Participant_ID)

#check the Abundance metric here, is it definitely relative?
brk_cls %>% group_by(Sample) %>% summarize(sum(Abundance)) #yes! :)

#check the sample IDs. are they all here?
ids <- brk_cls %>% select(Sample) %>% filter(!duplicated(Sample))
ids$Participant_ID <- sapply(strsplit(ids$Sample,"_"), `[`, 1)
ids$Participant_ID <- gsub("H", "HEAT_", ids$Participant_ID)
ids %>% group_by(Participant_ID) %>% count() #yes #looks very nice and complete
#write.csv(ids, "/ebio/abt3_projects2/uHEAT/data/metadata/ids.csv")
#ids <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/ids.csv", row.names=1)

#keep a tax table
tax <- brk_cls %>% select(Domain, Phylum, Class, Order, Family, Genus, Species, taxonomy_id) %>%
  filter(!duplicated(taxonomy_id))
#write.csv(tax, "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/kraken/tax_table.csv")

#per-sample metadata #note: mocks, blanks & samples with low read depth already excluded #so are old & infected samples
#updated: this version also excludes SERIES 2 of all double series.
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT/data/metadata/meta_seq_clean_series1.24.08.18.csv", 
                     row.names=1)

#number of clean samples
length(meta_seq$SampleID)

#per-person metadata, matching
meta <- meta_seq %>% filter(!duplicated(Participant_ID)&!is.na(fever7d_aboveavg)) %>%
  dplyr::select(-c(SampleID, Relative_Day, Time0, Time01))
dim(meta)

#define DIET 3
meta$Diet3 <- ifelse(meta$Diet %in% c("Vegan", "Vegetarian"), yes="Veg", no="Omni")
meta_seq$Diet3 <- ifelse(meta_seq$Diet %in% c("Vegan", "Vegetarian"), yes="Veg", no="Omni")

####chunk 4 maaslin2####
#make dataframes compatible for Maaslin2, at different taxonomic levels
feces_df <- brk_cls %>% select(Phylum, Family, Genus, Species, Sample, Abundance) %>% 
  pivot_wider(names_from=Sample, values_from=Abundance)
head(feces_df) 

#species #essentially strain level
feces_df_species <- feces_df %>% select(-c(Phylum, Family, Genus)) %>% as.data.frame()
row.names(feces_df_species) <- feces_df_species$Species
feces_df_species$Species <- NULL
meta_seq <- meta_seq[order(meta_seq$SampleID),] #order
feces_df_species <- feces_df_species[,which(names(feces_df_species)%in% meta_seq$SampleID)] #filter - no mocks, blanks
feces_df_species <- feces_df_species[,order(names(feces_df_species))] #order
identical(names(feces_df_species), meta_seq$SampleID)

#genus
feces_df_genus <- brk_cls %>% select(Genus, Sample, Abundance) %>%
  group_by(Sample, Genus) %>% summarize(Abundance=sum(Abundance)) %>%
  pivot_wider(names_from=Sample, values_from=Abundance) %>%
  as.data.frame()
head(feces_df_genus)
#cleanup
feces_df_genus$Genus[1] <- "Unknown"
row.names(feces_df_genus) <- feces_df_genus$Genus
feces_df_genus$Genus <- NULL
feces_df_genus <- feces_df_genus[,which(names(feces_df_genus)%in% meta_seq$SampleID)] #filter - no mocks, blanks
feces_df_genus <- feces_df_genus[,order(names(feces_df_genus))] #order
identical(names(feces_df_genus), meta_seq$SampleID)
colSums(feces_df_genus) #rounding errors
row.names(feces_df_genus)

#maybe I want this, e.g. for RNA-seq analysis
#write.csv(feces_df_genus, "/ebio/abt3_projects2/uHEAT/data/output_kraken_all/feces_df_genus.csv")
#feces_df_genus <- read.csv("/ebio/abt3_projects2/uHEAT/data/output_kraken_all/feces_df_genus.csv",
#                           row.names=1)

#later, for plotting
plot_gen <- data.frame(t(feces_df_genus), meta_seq)
plot_spec <- data.frame(t(feces_df_species), meta_seq)

row.names(meta_seq) <- meta_seq$SampleID
meta_seq$fever7d_aboveavg <- factor(meta_seq$fever7d_aboveavg, levels=c("lo", "hi"))
#meta_seq$Bristol_Stool_Score2 <- ifelse(meta_seq$Bristol_Stool_Score2=="bNA",
#                                        NA, meta_seq$Bristol_Stool_Score2) #BSS grouped
meta_seq$Bristol_Stool_Score <- as.numeric(meta_seq$Bristol_Stool_Score)

##Fever - overall (genus) #7d fever with technical covariates only
fit_data = Maaslin2(
  input_data = feces_df_genus, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  #reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #min abundance 10^-3; anything lower gets unbelievable
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergen <- fit_data$results
res_long_fevergen %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg") 

#export results
write.csv(res_long_fevergen, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAg_seqBSSTimeRT_fever7d.csv")

#ALSO EXPORT the final clean input table so I can plot them later
write.csv(feces_df_genus, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/DNA_feces_df_genus.csv")

##Fever - overall (species) #7d fever with technical covariates only
fit_data = Maaslin2(
  input_data = feces_df_species, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  #reference=c("Sex", "m"), #tell it to use plate 1 as the ref for plate compare / male as ref for Sex
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #min abundance 10^-3; anything lower gets unbelievable
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_feverspec <- fit_data$results
res_long_feverspec %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg") 
tax %>% dplyr::filter(Genus=="Acetatifactor") %>% as.data.frame() %>% head()

#export results
write.csv(res_long_feverspec, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAs_seqBSSTimeRT_fever7d.csv")

#ALSO EXPORT the final clean input table so I can plot them later
write.csv(feces_df_species, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/DNA_feces_df_species.csv")

#baseline only
meta_seq_baseline <- meta_seq %>% filter(Stool_timing3=="before")
feces_df_genus_baseline <- feces_df_genus[,which(names(feces_df_genus)%in%meta_seq_baseline$SampleID)]
identical(names(feces_df_genus_baseline), meta_seq_baseline$SampleID)
feces_df_species_baseline <- feces_df_species[,which(names(feces_df_species)%in%meta_seq_baseline$SampleID)]
identical(names(feces_df_species_baseline), meta_seq_baseline$SampleID)

row.names(meta_seq_baseline) <- meta_seq_baseline$SampleID

fit_data = Maaslin2(
  input_data = feces_df_genus_baseline, 
  input_metadata = meta_seq_baseline, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"),
  #reference=c("Stool_timing_fever", "before"), 
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #expand the threshold from 10^-3 to 10%-4 ?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05,
  plot_scatter=FALSE
)
res_long_fevergenbefore <- fit_data$results
res_long_fevergenbefore %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg")

#export results
write.csv(res_long_fevergenbefore, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAg_before_seqBSSTimeRT_fever7d.csv")
#ALSO EXPORT the final clean input table so I can plot them later
write.csv(feces_df_species, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/DNA_feces_df_species.csv")



fit_data = Maaslin2(
  input_data = feces_df_species_baseline, 
  input_metadata = meta_seq_baseline, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"),
  #reference=c("Stool_timing_fever", "before"), 
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #expand the threshold from 10^-3 to 10%-4 ?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05,
  plot_scatter=FALSE
)
res_long_feverspecbefore <- fit_data$results
res_long_feverspecbefore %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg")

#export results
write.csv(res_long_feverspecbefore, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAs_before_seqBSSTimeRT_fever7d.csv")


####chunk 5 maaslin2 corrected####

##Fever - overall (genus) #7d fever with additional covariates only
fit_data = Maaslin2(
  input_data = feces_df_genus, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age", "Diet3", 
                    "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #min abundance 10^-3; anything lower gets unbelievable
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergen_corr <- fit_data$results
res_long_fevergen_corr %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg") 

#export results
write.csv(res_long_fevergen_corr, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAg_seqBSSTimeRT_SexAgeDiet3_fever7d.csv")

##Fever - overall (species) #7d fever with technical covariates only
fit_data = Maaslin2(
  input_data = feces_df_species, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age", "Diet3", 
                    "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #min abundance 10^-3; anything lower gets unbelievable
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_feverspec_corr <- fit_data$results
res_long_feverspec_corr %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg") 

#export results
write.csv(res_long_feverspec_corr, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAs_seqBSSTimeRT_SexAgeDiet3_fever7d.csv")

#baseline only
meta_seq_baseline <- meta_seq %>% filter(Stool_timing3=="before")
feces_df_genus_baseline <- feces_df_genus[,which(names(feces_df_genus)%in%meta_seq_baseline$SampleID)]
identical(names(feces_df_genus_baseline), meta_seq_baseline$SampleID)
feces_df_species_baseline <- feces_df_species[,which(names(feces_df_species)%in%meta_seq_baseline$SampleID)]
identical(names(feces_df_species_baseline), meta_seq_baseline$SampleID)

row.names(meta_seq_baseline) <- meta_seq_baseline$SampleID

fit_data = Maaslin2(
  input_data = feces_df_genus_baseline, 
  input_metadata = meta_seq_baseline, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age", "Diet3", 
                    "fever7d_aboveavg"),
  random_effects = c("Participant_ID"),
  reference=c("Sex", "m"), 
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #expand the threshold from 10^-3 to 10%-4 ?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05,
  plot_scatter=FALSE
)
res_long_fevergenbefore_corr <- fit_data$results
res_long_fevergenbefore_corr %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg")

#export results
write.csv(res_long_fevergenbefore_corr, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAg_before_seqBSSTimeRT_SexAgeDiet3_fever7d.csv")

fit_data = Maaslin2(
  input_data = feces_df_species_baseline, 
  input_metadata = meta_seq_baseline, 
  output = "maaslin2_output", 
  fixed_effects = c("seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age", "Diet3", 
                    "fever7d_aboveavg"),
  random_effects = c("Participant_ID"),
  reference=c("Sex", "m"), 
  normalization="tss", #total sum scaling
  min_abundance=0.0001, #expand the threshold from 10^-3 to 10%-4 ?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05,
  plot_scatter=FALSE
)
res_long_feverspecbefore_corr <- fit_data$results
res_long_feverspecbefore_corr %>% dplyr::filter(qval<0.05&metadata=="fever7d_aboveavg")

#export results
write.csv(res_long_feverspecbefore_corr, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/DNAs_before_seqBSSTimeRT_SexAgeDiet3_fever7d.csv")
