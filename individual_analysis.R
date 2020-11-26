# de testing for individual datasets that were filtered before

library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)

setwd("data/")

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files (already made so load txi object)
load("~/Documents/S_L_DH/samples/tximpost_all.RData")
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)
colnames(colData)<-c("Sample", "Rep", "Treatment", "Timepoint", "Genotype", "Factor", "Condition", "Resistant")
filter<-read.csv(file="~/Documents/S_L_DH/data/all_lists_filtered.csv")
counts<-txi.kallisto.tsv$counts

colData<-colData%>% 
  mutate(T2=case_when(
    Timepoint == 6 ~ "1",
    Timepoint == 24 ~ "2",
    Timepoint == 48 ~ "3",
    Timepoint == 96 ~ "4",
  ))

list<-list()
degs_filtererd<-list()
degs_no_filter<-list()
degs_filtered_post<-list()
for(set in unique(filter$Comparison)){
  comp<-filter[(filter$Comparison==set),]
  genotype<-strsplit(set, split = "")[[1]][1]
  timepoint<-strsplit(set, split = "")[[1]][2]
  metadata<-colData[(colData$Genotype==genotype),,drop=T]
  metadata<-metadata[(metadata$T2==timepoint),,drop=T]
  metadata<-droplevels(metadata)
  samples<-as.character(metadata$Sample)
  genes<-as.character(comp$GeneID)
  df<-counts[genes,samples]
  case<-colnames(df)==metadata$Sample
  list[[length(list)+1]]<-case
  dds <- DESeqDataSetFromMatrix(round(df), metadata, ~ Treatment)
  dds<-DESeq(dds)
  r<-as.data.frame(results(dds, format='DataFrame', tidy=TRUE))
  r$comparison<-paste(set)
  degs_filtererd[[length(degs_filtererd)+1]]<-r
  df2<-counts[,samples]
  case<-colnames(df2)==metadata$Sample
  list[[length(list)+1]]<-case
  dds <- DESeqDataSetFromMatrix(round(df2), metadata, ~ Treatment)
  dds<-DESeq(dds)
  r<-as.data.frame(results(dds, format='DataFrame', tidy=TRUE))
  r$comparison<-paste(set)
  degs_no_filter[[length(degs_no_filter)+1]]<-r
  r_filtered<-subset(r, r$row %in% comp$GeneID)
  degs_filtered_post[[length(degs_filtered_post)+1]]<-r_filtered
  
}
all_filtered<-as.data.frame(do.call(rbind.data.frame, degs_filtererd))
all_no_filter<-as.data.frame(do.call(rbind.data.frame, degs_no_filter))
all_filtered_post<-as.data.frame(do.call(rbind.data.frame, degs_filtered_post))
all_filtered<-all_filtered[(all_filtered$padj<0.05),]
all_filtered<-na.omit(all_filtered)
all_no_filter<-all_no_filter[(all_no_filter$padj<0.05),]
all_no_filter<-na.omit(all_no_filter)
all_filtered_post<-all_filtered_post[(all_filtered_post$padj<0.05),]
all_filtered_post<-na.omit(all_filtered_post)


all_filtered<-all_filtered%>% 
  mutate(Timepoint=case_when(
comparison == "G1" ~ 6,
comparison == "G2" ~ 24,
comparison == "G3" ~ 48,
comparison == "G4" ~ 69,
comparison == "L1" ~ 6,
comparison == "L2" ~ 24,
comparison == "L3" ~ 48,
comparison == "L4" ~ 69,
comparison == "R1" ~ 6,
comparison == "R2" ~ 24,
comparison == "R3" ~ 48,
comparison == "R4" ~ 69,
comparison == "S1" ~ 6,
comparison == "S2" ~ 24,
comparison == "S3" ~ 48,
comparison == "S4" ~ 69
))

all_filtered<-all_filtered%>% 
  mutate(Genotype=case_when(
    comparison == "G1" ~ "Stigg",
    comparison == "G2" ~ "Stigg",
    comparison == "G3" ~ "Stigg",
    comparison == "G4" ~ "Stigg",
    comparison == "L1" ~ "Longbow",
    comparison == "L2" ~ "Longbow",
    comparison == "L3" ~ "Longbow",
    comparison == "L4" ~ "Longbow",
    comparison == "R1" ~ "Resistant",
    comparison == "R2" ~ "Resistant",
    comparison == "R3" ~ "Resistant",
    comparison == "R4" ~ "Resistant",
    comparison == "S1" ~ "Susceptible",
    comparison == "S2" ~ "Susceptible",
    comparison == "S3" ~ "Susceptible",
    comparison == "S4" ~ "Susceptible"
  ))


all_filtered_post<-all_filtered_post%>% 
  mutate(Timepoint=case_when(
    comparison == "G1" ~ 6,
    comparison == "G2" ~ 24,
    comparison == "G3" ~ 48,
    comparison == "G4" ~ 69,
    comparison == "L1" ~ 6,
    comparison == "L2" ~ 24,
    comparison == "L3" ~ 48,
    comparison == "L4" ~ 69,
    comparison == "R1" ~ 6,
    comparison == "R2" ~ 24,
    comparison == "R3" ~ 48,
    comparison == "R4" ~ 69,
    comparison == "S1" ~ 6,
    comparison == "S2" ~ 24,
    comparison == "S3" ~ 48,
    comparison == "S4" ~ 69
  ))

all_filtered_post<-all_filtered_post%>% 
  mutate(Genotype=case_when(
    comparison == "G1" ~ "Stigg",
    comparison == "G2" ~ "Stigg",
    comparison == "G3" ~ "Stigg",
    comparison == "G4" ~ "Stigg",
    comparison == "L1" ~ "Longbow",
    comparison == "L2" ~ "Longbow",
    comparison == "L3" ~ "Longbow",
    comparison == "L4" ~ "Longbow",
    comparison == "R1" ~ "Resistant",
    comparison == "R2" ~ "Resistant",
    comparison == "R3" ~ "Resistant",
    comparison == "R4" ~ "Resistant",
    comparison == "S1" ~ "Susceptible",
    comparison == "S2" ~ "Susceptible",
    comparison == "S3" ~ "Susceptible",
    comparison == "S4" ~ "Susceptible"
  ))


all_no_filter<-all_no_filter%>% 
  mutate(Timepoint=case_when(
    comparison == "G1" ~ 6,
    comparison == "G2" ~ 24,
    comparison == "G3" ~ 48,
    comparison == "G4" ~ 69,
    comparison == "L1" ~ 6,
    comparison == "L2" ~ 24,
    comparison == "L3" ~ 48,
    comparison == "L4" ~ 69,
    comparison == "R1" ~ 6,
    comparison == "R2" ~ 24,
    comparison == "R3" ~ 48,
    comparison == "R4" ~ 69,
    comparison == "S1" ~ 6,
    comparison == "S2" ~ 24,
    comparison == "S3" ~ 48,
    comparison == "S4" ~ 69
  ))

all_no_filter<-all_no_filter%>% 
  mutate(Genotype=case_when(
    comparison == "G1" ~ "Stigg",
    comparison == "G2" ~ "Stigg",
    comparison == "G3" ~ "Stigg",
    comparison == "G4" ~ "Stigg",
    comparison == "L1" ~ "Longbow",
    comparison == "L2" ~ "Longbow",
    comparison == "L3" ~ "Longbow",
    comparison == "L4" ~ "Longbow",
    comparison == "R1" ~ "Resistant",
    comparison == "R2" ~ "Resistant",
    comparison == "R3" ~ "Resistant",
    comparison == "R4" ~ "Resistant",
    comparison == "S1" ~ "Susceptible",
    comparison == "S2" ~ "Susceptible",
    comparison == "S3" ~ "Susceptible",
    comparison == "S4" ~ "Susceptible"
  ))

table(all_filtered$Genotype, all_filtered$Timepoint)
table(all_filtered_post$Genotype, all_filtered_post$Timepoint)
table(all_no_filter$Genotype, all_no_filter$Timepoint)


write.csv(all_filtered, "~/Documents/S_L_DH/data/DE_tests/separate_filtered.csv")
write.csv(all_filtered_post, "~/Documents/S_L_DH/data/DE_tests/separate_filtered_post.csv")
write.csv(all_no_filter, "~/Documents/S_L_DH/data/DE_tests/separate_no_filter.csv")



all_sig<-all_sig%>% 
  mutate(Timepoint=case_when(
    comparison == "G1" ~ 6,
    comparison == "G2" ~ 24,
    comparison == "G3" ~ 48,
    comparison == "G4" ~ 69,
    comparison == "L1" ~ 6,
    comparison == "L2" ~ 24,
    comparison == "L3" ~ 48,
    comparison == "L4" ~ 69,
    comparison == "R1" ~ 6,
    comparison == "R2" ~ 24,
    comparison == "R3" ~ 48,
    comparison == "R4" ~ 69,
    comparison == "S1" ~ 6,
    comparison == "S2" ~ 24,
    comparison == "S3" ~ 48,
    comparison == "S4" ~ 69
  ))

all_sig<-all_sig%>% 
  mutate(Genotype=case_when(
    comparison == "G1" ~ "Stigg",
    comparison == "G2" ~ "Stigg",
    comparison == "G3" ~ "Stigg",
    comparison == "G4" ~ "Stigg",
    comparison == "L1" ~ "Longbow",
    comparison == "L2" ~ "Longbow",
    comparison == "L3" ~ "Longbow",
    comparison == "L4" ~ "Longbow",
    comparison == "R1" ~ "Resistant",
    comparison == "R2" ~ "Resistant",
    comparison == "R3" ~ "Resistant",
    comparison == "R4" ~ "Resistant",
    comparison == "S1" ~ "Susceptible",
    comparison == "S2" ~ "Susceptible",
    comparison == "S3" ~ "Susceptible",
    comparison == "S4" ~ "Susceptible"
  ))
table(all_sig$Genotype, all_sig$Timepoint)

all_sig_nofilter<-all_sig
write.csv(all_sig_nofilter, "DE_tests/separate_no_filter.csv")

