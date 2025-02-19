#Transcriptomics Code Executable
#libraries
library(tidyr)
library(data.table)
library(dplyr)
library(Maaslin2)
library(ggplot2)
source("/ebio/abt3_projects2/uHEAT/code/R_code/LeyLabRMisc/R/init.R")
####chunk 1 import####
#7d fever with RNA seq depth, etc - updated
#this version now includes ONLY series1
meta_seq <- read.csv("/ebio/abt3_projects2/uHEAT2/data/metadata/meta_seq_rna_fever7d_series1_24.08.18.csv") 

#define DIET 3
meta_seq$Diet3 <- ifelse(meta_seq$Diet %in% c("Vegan", "Vegetarian"), yes="Veg", no="Omni")

#genes - rna - unstratified
hm3_genes <- fread("/ebio/abt3_projects2/uHEAT2/data/output_LLMGP_RNA_all/humann3/normalized/unstratified/regroup/annot/genefamilies_uniref50_default.tsv")
names(hm3_genes) <- gsub("_Abundance-RPKs", "", names(hm3_genes), fixed=TRUE)

hm3_genes <- hm3_genes %>% dplyr::rename(Gene=`# Gene Family`)
row.names(hm3_genes) <- hm3_genes$Gene
gene_key <- hm3_genes %>% select(c(Gene, annotation)) #but then correct it to go with Maaslin2
gene_key$Gene_unfixed <- gene_key$Gene
gene_key$Gene <- gsub("-", ".", gene_key$Gene, fixed=TRUE)
gene_key$Gene <- gsub(";", ".", gene_key$Gene, fixed=TRUE)
gene_key$Gene <- gsub(":", ".", gene_key$Gene, fixed=TRUE)
gene_key$Gene <- gsub(" ", ".", gene_key$Gene, fixed=TRUE)

hm3_genes_M <- hm3_genes %>% as.data.frame() %>% dplyr::select(-c("annotation", "Gene")) #to refilter, start again here.
head(hm3_genes_M)
row.names(hm3_genes_M) <- hm3_genes$Gene

meta_seq <- meta_seq[order(meta_seq$SampleID),] #order
meta_seq <- meta_seq[which(meta_seq$SampleID %in% names(hm3_genes_M)),] #if I didn't filter above - make sure I have only the RNA samples, not all of them
hm3_genes_M <- hm3_genes_M[,which(names(hm3_genes_M) %in% meta_seq$SampleID)] #and if it's the stringent metadata, kick out the infected samples too
hm3_genes_M <- hm3_genes_M[,order(names(hm3_genes_M))] #order
print("identical(names(hm3_genes_M), meta_seq$SampleID)") #for the log file
identical(names(hm3_genes_M), meta_seq$SampleID) #cleaned up again
row.names(meta_seq) <- meta_seq$SampleID

meta_seq$fever7d_aboveavg <- factor(meta_seq$fever7d_aboveavg, levels=c("lo", "hi"))

####chunk 2 rna test####
##Fever - overall #7d fever with technical covariates only
fit_data = Maaslin2(
  input_data = hm3_genes_M, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("RNA_seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  #reference=c("Sex", "m"), #tell it to use plate 1 as the ref for plate compare / male as ref for Sex
  normalization="tss", #total sum scaling
  min_abundance=0.00001, #expand the threshold from 10^-3 to 10%-4 because genes?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergenes3 <- fit_data$results
res_long_fevergenes3 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") 
hits3b <- fit_data$results %>% filter(pval<0.05&metadata=="fever7d_aboveavg")
hm3_genes[which(hm3_genes$Gene %in% hits3b$feature),] %>% select(c(Gene, annotation)) %>% as.data.frame()
#export 
f7genes <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(hits3b) %>% filter(!is.na(coef))
f7genes
#write.csv(f7genes, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNAUniref50_Fever7d_aboveavg_hits.csv")

#vs: write out entire table, e.g. for a volcano plot.
res_long_fevergenes3 <- fit_data$results
rna_all <- fit_data$results %>% filter(metadata=="fever7d_aboveavg")
#export 
f7genes_all <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(rna_all) %>% filter(!is.na(coef))
f7genes_all
#write.csv(f7genes_all, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNAUniref50_Fever7d_aboveavg_all.csv")

#additionally write out the input table for later
#write.csv(hm3_genes_M, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/RNAUniref50_hm3_genes_M.csv")

#baseline only
meta_seq_baseline <- meta_seq %>% filter(Stool_timing3=="before")
hm3_genes_M_baseline <- hm3_genes_M[,which(names(hm3_genes_M)%in%meta_seq_baseline$SampleID)]
identical(names(hm3_genes_M_baseline), meta_seq_baseline$SampleID)
row.names(meta_seq_baseline) <- meta_seq_baseline$SampleID


####chunk 3 dna vs rna create####
#matching DNA models run on the same samples.
#genes - unstratified
hm3_genes_dna <- fread("/ebio/abt3_projects2/uHEAT/data/output_LLMGP_humann_all_rep/output_LLMGP_humann_pool3/humann3/normalized/unstratified/regroup/annot/genefamilies_uniref50_default.tsv")
names(hm3_genes_dna) <- gsub("_Abundance-RPKs", "", names(hm3_genes_dna), fixed=TRUE)
head(hm3_genes_dna)
hm3_genes_dna <- hm3_genes_dna %>% rename(Gene=`# Gene Family`)

hm3_genes_dna_M <- hm3_genes_dna %>% select(-c("annotation", "Gene"))
head(hm3_genes_dna_M)
hm3_genes_dna_M <- as.data.frame(hm3_genes_dna_M)
row.names(hm3_genes_dna_M) <- hm3_genes_dna$Gene

meta_seq <- meta_seq[order(meta_seq$SampleID),] #order
#meta_seq <- meta_seq[which(meta_seq$SampleID %in% names(hm3_genes_dna_M)),] #if I didn't filter above - make sure I have only the RNA samples, not all of them
hm3_genes_dna_M <- hm3_genes_dna_M[,which(names(hm3_genes_dna_M) %in% meta_seq$SampleID)] #kick out all samples for which no RNA-seq.
hm3_genes_dna_M <- hm3_genes_dna_M[,order(names(hm3_genes_dna_M))] #order
paste("identical(names(hm3_genes_dna_M), meta_seq$SampleID)")
identical(names(hm3_genes_dna_M), meta_seq$SampleID) #cleaned up again
row.names(meta_seq) <- meta_seq$SampleID
row.names(hm3_genes_dna_M)

hm3_genes_dna <- NULL #because it's so big, get rid of it!
###OK. Relative RNA to DNA expression, across all genes. Los gehts!
print("identical(names(hm3_genes_M), names(hm3_genes_dna_M))") #for the log file
identical(names(hm3_genes_M), names(hm3_genes_dna_M))
hm3_genes_dna_M_filt <- hm3_genes_dna_M[which(row.names(hm3_genes_dna_M)%in%row.names(hm3_genes_M)),] #select only genes which also appear in RNA
hm3_genes_M_filt <- hm3_genes_M[which(row.names(hm3_genes_M)%in%row.names(hm3_genes_dna_M_filt)),] #and vice versa
print("identical(row.names(hm3_genes_dna_M_filt), row.names(hm3_genes_M_filt))")
identical(row.names(hm3_genes_dna_M_filt), row.names(hm3_genes_M_filt))
hm3_genes_rel <- hm3_genes_M_filt/(hm3_genes_dna_M_filt+0.0000001) #entire dataframe at once with pseudocount
print("identical(names(hm3_genes_rel), meta_seq$SampleID)")
identical(names(hm3_genes_rel), meta_seq$SampleID)

#for sequencing depth: what would I control for, here? both? relative sequencing depth?
meta_seq$RNA_DNA_seq_depth <- meta_seq$RNA_seq_depth / meta_seq$DNA_seq_depth 
meta_seq <- meta_seq[which(meta_seq$SampleID %in% names(hm3_genes_rel)),] #if I didn't filter above - make sure I have only the RNA samples, not all of them
meta_seq$RNA_DNA_seq_depth

print("identical(names(hm3_genes_rel), meta_seq$SampleID)")
identical(names(hm3_genes_rel), meta_seq$SampleID)
row.names(meta_seq) <- meta_seq$SampleID
print("row.names(hm3_genes_rel)")
row.names(hm3_genes_rel)

####chunk 4 dna vs rna test####
#now, test via Maaslin2: without TSS, because it makes no sense for ratios
#WARNING: this takes ≥10 hours
fit_data = Maaslin2(
  input_data = hm3_genes_rel, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("RNA_DNA_seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"), #use relative seq depth, RNA vs DNA
  random_effects = c("Participant_ID"), 
  normalization="NONE", #no additional normalizaation for ratios
  #min_abundance=0.00001, #no reasonable min abund for ratios
  min_prevalence=0.25, #must be present in 25% of samples #may want to increase this
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergenes_rnadna <- fit_data$results 
res_long_fevergenes_rnadna %>% filter(qval<0.05&metadata=="fever7d_aboveavg") #lots.. .and mostly the same as before

#export immediately
hits_rnadna <- res_long_fevergenes_rnadna %>% filter(pval<0.05&metadata=="fever7d_aboveavg")
hits_rnadna
genes_rnadna <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(hits_rnadna) %>% filter(!is.na(coef))
genes_rnadna
write.csv(genes_rnadna, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNADNAUniref50_Fever7d_aboveavg_hits.csv") 

#versus: write out entire table, e.g. for a volcano plot.
all_rnadna <- res_long_fevergenes_rnadna %>% filter(metadata=="fever7d_aboveavg")
all_rnadna
genes_rnadna_all <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(all_rnadna) %>% filter(!is.na(coef))
genes_rnadna_all
write.csv(genes_rnadna_all, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNADNAUniref50_Fever7d_aboveavg_all.csv") 

#additionally write out the input table for later
#write.csv(hm3_genes_rel, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/RNADNAUniref50_hm3_genes_rel.csv")


####chunk 4 dna vs rna Age, Sex ####
#now, test via Maaslin2: without TSS, because it makes no sense for ratios
#WARNING: this takes ≥10 hours
fit_data = Maaslin2(
  input_data = hm3_genes_rel, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("RNA_DNA_seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age",
                    "fever7d_aboveavg"), #use relative seq depth, RNA vs DNA
  random_effects = c("Participant_ID"), 
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="NONE", #no additional normalizaation for ratios
  #min_abundance=0.00001, #no reasonable min abund for ratios
  min_prevalence=0.25, #must be present in 25% of samples #may want to increase this
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergenes2 <- fit_data$results 
res_long_fevergenes2 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") #lots.. .and mostly the same as before

#versus: write out entire table, e.g. for a volcano plot.
all_rnadna2 <- res_long_fevergenes2 %>% filter(metadata=="fever7d_aboveavg")

genes_rnadna2 <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(all_rnadna2) %>% filter(!is.na(coef))
genes_rnadna2
write.csv(genes_rnadna2, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNADNAUniref50_SexAgeFever7d_all.csv") 


####chunk 4 dna vs rna Age, Sex, Diet ####
#now, test via Maaslin2: without TSS, because it makes no sense for ratios
#WARNING: this takes ≥10 hours
fit_data = Maaslin2(
  input_data = hm3_genes_rel, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("RNA_DNA_seq_depth", "Bristol_Stool_Score", "Time_at_RT", 
                    "Sex", "Age", "Diet3",
                    "fever7d_aboveavg"), #use relative seq depth, RNA vs DNA
  random_effects = c("Participant_ID"), 
  reference=c("Sex", "m"), #tell it to use male as ref for Sex
  normalization="NONE", #no additional normalizaation for ratios
  #min_abundance=0.00001, #no reasonable min abund for ratios
  min_prevalence=0.25, #must be present in 25% of samples #may want to increase this
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
res_long_fevergenes3 <- fit_data$results 
res_long_fevergenes3 %>% filter(qval<0.05&metadata=="fever7d_aboveavg") #lots.. .and mostly the same as before

#versus: write out entire table, e.g. for a volcano plot.
all_rnadna3 <- res_long_fevergenes3 %>% filter(metadata=="fever7d_aboveavg")

genes_rnadna3 <- gene_key %>% dplyr::rename(feature=Gene) %>% 
  full_join(all_rnadna3) %>% filter(!is.na(coef))
genes_rnadna3
write.csv(genes_rnadna3, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/RNADNAUniref50_SexAgeDiet3Fever7d_all.csv") 

####chunk 5 stratified flagellins####

#genes - stratified
hm3_genes_strat <- fread("/ebio/abt3_projects2/uHEAT2/data/output_LLMGP_RNA_all/humann3/normalized/stratified/regroup/annot/genefamilies_uniref50_default.tsv")
names(hm3_genes_strat) <- gsub("_Abundance-RPKs", "", names(hm3_genes_strat), fixed=TRUE)

hm3_genes_strat <- hm3_genes_strat %>% rename(Gene=`# Gene Family`)

#pick (real) flagellins only
flag <- hm3_genes_strat %>% filter(annotation=="Flagellin")
flag %>% select(Gene, annotation) %>% as.data.frame()
flag[grep("FliC", flag$annotation),] 

#save annotations
flag_key <- flag %>% select(c(Gene, annotation)) #but then correct it to go with Maaslin2
flag_key$rownames <- paste("V", row.names(flag_key), sep="")
flag_key$Gene_unfixed <- flag_key$Gene
flag_key$Gene <- gsub("-", ".", flag_key$Gene, fixed=TRUE)
flag_key$Gene <- gsub(";", ".", flag_key$Gene, fixed=TRUE)
flag_key$Gene <- gsub(":", ".", flag_key$Gene, fixed=TRUE)
flag_key$Gene <- gsub(" ", ".", flag_key$Gene, fixed=TRUE)
flag_key$Gene <- gsub("|", ".", flag_key$Gene, fixed=TRUE)
head(as.data.frame(flag_key))
row.names(flag) <- flag_key$Gene

#get rid of the massive original file
hm3_genes_strat <- NULL

#fix the double participant ID problem 

#flag_M <- flag %>% as.data.frame() %>% select(-c("annotation", "Gene"))
#head(as.data.frame(flag_M))

meta_seq <- meta_seq[order(meta_seq$SampleID),] #order
flag_M <- as.data.frame(flag)[,which(names(flag) %in% meta_seq$SampleID)]
flag_M <- flag_M[,order(names(flag_M))] #order
row.names(flag_M) <- flag_key$Gene
print("identical(names(flag_M), meta_seq$SampleID)")
identical(names(flag_M), meta_seq$SampleID)
row.names(meta_seq) <- meta_seq$SampleID
head(flag_M)

#with fever7d (RNA level)
fit_data = Maaslin2(
  input_data = flag_M, 
  input_metadata = meta_seq, 
  output = "maaslin2_output", 
  fixed_effects = c("RNA_seq_depth", "Bristol_Stool_Score", "Time_at_RT", "fever7d_aboveavg"),
  random_effects = c("Participant_ID"), 
  #reference=c("Sex", "m"), #tell it to use plate 1 as the ref for plate compare / male as ref for Sex
  normalization="tss", #total sum scaling
  #min_abundance=0.00001, #no abundance threshold here?
  min_prevalence=0.25, #must be present in 25% of samples
  max_significance=0.05, #calm down with the non significant associations and plots
  plot_scatter=FALSE
)
flag_fevergenes <- fit_data$results 
hits7 <- flag_fevergenes %>% filter(qval<0.05&metadata=="fever7d_aboveavg")
hits7b <- flag_fevergenes %>% filter(pval<0.05&metadata=="fever7d_aboveavg") #slightly more relaxed?
hits7b #yes, they are all Lachnos. Of course they are all Lachnos.

#export (no more annotation needed??)
#write.csv(hits7b, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_output/Uniref50_fever7d_strat_hits.csv")

#additionally write out the input table for later
#write.csv(flag_M, "/ebio/abt3_projects2/uHEAT/code/R_code/uHEAT_CleanMarkdowns_Diet3/maaslin2_input/RNAUniref50_flag_M.csv")
