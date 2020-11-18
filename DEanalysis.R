library(tximport)
library(DESeq2)
library(tidyr)
library(plyr)
library(dplyr)
library(gplots)

# samples<-as.data.frame(dir("samples/"))
# samples<-samples[(5:52),]
# samples<-as.data.frame(samples)
# files <- file.path(samples$samples, "abundance.h5")
# names(files) <- paste0(samples$samples)
# all(file.exists(files))
# txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)

load("~/Documents/S_L_DH/samples/tximpost_all.RData")

counts<-txi.kallisto.tsv$counts
abundance<-txi.kallisto.tsv$abundance
# DESeq2 differentially expression testing
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor)
dds<-DESeq(dds)

RES_S1<-as.data.frame(results(dds, contrast=c("Factor","T1S","C1S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R1<-as.data.frame(results(dds, contrast=c("Factor","T1R","C1R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_G1<-as.data.frame(results(dds, contrast=c("Factor","T1G","C1G"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_L1<-as.data.frame(results(dds, contrast=c("Factor","T1L","C1L"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R2<-as.data.frame(results(dds, contrast=c("Factor","T2R","C2R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_L2<-as.data.frame(results(dds, contrast=c("Factor","T2L","C2L"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_G2<-as.data.frame(results(dds, contrast=c("Factor","T2G","C2G"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S2<-as.data.frame(results(dds, contrast=c("Factor","T2S","C2S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S3<-as.data.frame(results(dds, contrast=c("Factor","T3S","C3S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R3<-as.data.frame(results(dds, contrast=c("Factor","T3R","C3R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_G3<-as.data.frame(results(dds, contrast=c("Factor","T3G","C3G"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_L3<-as.data.frame(results(dds, contrast=c("Factor","T3L","C3L"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S4<-as.data.frame(results(dds, contrast=c("Factor","T4S","C4S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R4<-as.data.frame(results(dds, contrast=c("Factor","T4R","C4R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_G4<-as.data.frame(results(dds, contrast=c("Factor","T4G","C4G"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_L4<-as.data.frame(results(dds, contrast=c("Factor","T4L","C4L"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))


RES_S1$Cultivar<-"Susceptible"
RES_R1$Cultivar<-"Resistant"
RES_G1$Cultivar<-"Stigg"
RES_L1$Cultivar<-"Longbow"
RES_R2$Cultivar<-"Resistant"
RES_L2$Cultivar<-"Longbow"
RES_G2$Cultivar<-"Stigg"
RES_S2$Cultivar<-"Susceptible"
RES_S3$Cultivar<-"Susceptible"
RES_R3$Cultivar<-"Resistant"
RES_G3$Cultivar<-"Stigg"
RES_L3$Cultivar<-"Longbow"
RES_S4$Cultivar<-"Susceptible"
RES_R4$Cultivar<-"Resistant"
RES_G4$Cultivar<-"Stigg"
RES_L4$Cultivar<-"Longbow"


RES_G1$Timepoint<-6
RES_G2$Timepoint<-24
RES_G3$Timepoint<-48
RES_G4$Timepoint<-96
RES_L1$Timepoint<-6
RES_L2$Timepoint<-24
RES_L3$Timepoint<-48
RES_L4$Timepoint<-96
RES_R1$Timepoint<-6
RES_R2$Timepoint<-24
RES_R3$Timepoint<-48
RES_R4$Timepoint<-96
RES_S1$Timepoint<-6
RES_S2$Timepoint<-24
RES_S3$Timepoint<-48
RES_S4$Timepoint<-96


RES_G1_filtered<-subset(RES_G1, RES_G1$row %in% G1$GeneID)
RES_G2_filtered<-subset(RES_G2, RES_G2$row %in% G2$GeneID)
RES_G3_filtered<-subset(RES_G3, RES_G3$row %in% G3$GeneID)
RES_G4_filtered<-subset(RES_G4, RES_G4$row %in% G4$GeneID)
RES_L1_filtered<-subset(RES_L1, RES_L1$row %in% L1$GeneID)
RES_L2_filtered<-subset(RES_L2, RES_L2$row %in% L2$GeneID)
RES_L3_filtered<-subset(RES_L3, RES_L3$row %in% L3$GeneID)
RES_L4_filtered<-subset(RES_L4, RES_L4$row %in% L4$GeneID)
RES_R1_filtered<-subset(RES_R1, RES_R1$row %in% R1$GeneID)
RES_R2_filtered<-subset(RES_R2, RES_R2$row %in% R2$GeneID)
RES_R3_filtered<-subset(RES_R3, RES_R3$row %in% R3$GeneID)
RES_R4_filtered<-subset(RES_R4, RES_R4$row %in% R4$GeneID)
RES_S1_filtered<-subset(RES_S1, RES_S1$row %in% S1$GeneID)
RES_S2_filtered<-subset(RES_S2, RES_S2$row %in% S2$GeneID)
RES_S3_filtered<-subset(RES_S3, RES_S3$row %in% S3$GeneID)
RES_S4_filtered<-subset(RES_S4, RES_S4$row %in% S4$GeneID)

all<-rbind(RES_G1,
           RES_G2,
           RES_G3,
           RES_G4,
           RES_L1,
           RES_L2,
           RES_L3,
           RES_L4,
           RES_R1,
           RES_R2,
           RES_R3,
           RES_R4,
           RES_S1,
           RES_S2,
           RES_S3,
           RES_S4)

all_filtered<-rbind(RES_G1_filtered,
           RES_G2_filtered,
           RES_G3_filtered,
           RES_G4_filtered,
           RES_L1_filtered,
           RES_L2_filtered,
           RES_L3_filtered,
           RES_L4_filtered,
           RES_R1_filtered,
           RES_R2_filtered,
           RES_R3_filtered,
           RES_R4_filtered,
           RES_S1_filtered,
           RES_S2_filtered,
           RES_S3_filtered,
           RES_S4_filtered)

all_sig_filtered<-all_filtered[(all_filtered$padj<0.05),]
all_sig<-all[(all$padj<0.05),]
all_sig<-na.omit(all_sig)
all_sig_filtered<-na.omit(all_sig_filtered)

table(all_sig$Cultivar, all_sig$Timepoint)
table(all_sig_filtered$Cultivar, all_sig_filtered$Timepoint)

write.csv(all_sig, "~/Documents/S_L_DH/data/all_significant_filtered.csv")

T6<-all[(all$Timepoint==6),]
T24<-all[(all$Timepoint==24),]
T48<-all[(all$Timepoint==48),]
T96<-all[(all$Timepoint==96),]


write.csv(RES_G1, file="~/Documents/S_L_DH/RES_G1.csv", row.names=F)
write.csv(RES_G2, file="~/Documents/S_L_DH/RES_G2.csv", row.names=F)
write.csv(RES_G3, file="~/Documents/S_L_DH/RES_G3.csv", row.names=F)
write.csv(RES_G4, file="~/Documents/S_L_DH/RES_G4.csv", row.names=F)
write.csv(RES_L1, file="~/Documents/S_L_DH/RES_L1.csv", row.names=F)
write.csv(RES_L2, file="~/Documents/S_L_DH/RES_L2.csv", row.names=F)
write.csv(RES_L3, file="~/Documents/S_L_DH/RES_L3.csv", row.names=F)
write.csv(RES_L4, file="~/Documents/S_L_DH/RES_L4.csv", row.names=F)
write.csv(RES_R1, file="~/Documents/S_L_DH/RES_R1.csv", row.names=F)
write.csv(RES_R2, file="~/Documents/S_L_DH/RES_R2.csv", row.names=F)
write.csv(RES_R3, file="~/Documents/S_L_DH/RES_R3.csv", row.names=F)
write.csv(RES_R4, file="~/Documents/S_L_DH/RES_R4.csv", row.names=F)
write.csv(RES_S1, file="~/Documents/S_L_DH/RES_S1.csv", row.names=F)
write.csv(RES_S2, file="~/Documents/S_L_DH/RES_S2.csv", row.names=F)
write.csv(RES_S3, file="~/Documents/S_L_DH/RES_S3.csv", row.names=F)
write.csv(RES_S4, file="~/Documents/S_L_DH/RES_S4.csv", row.names=F)

write.csv(T6, "~/Documents/S_L_DH/T6.csv")
write.csv(T24, "~/Documents/S_L_DH/T24.csv")
write.csv(T48, "~/Documents/S_L_DH/T48.csv")
write.csv(T96, "~/Documents/S_L_DH/T96.csv")


Kat_genes<-unique(all_sig$row)
kat<-subset(all, all$row %in% Kat_genes)
kat2<-kat[,c(1, 3, 8, 9)]
kat2$Factor<-paste(kat2$Cultivar, kat2$Timepoint)
kat2<-kat2[,c(1, 5, 2)]
test<-spread(kat2, key="Factor", value="log2FoldChange", fill=0)
row.names(test)<-test[,1]
test$row<-NULL
test<-as.matrix(test)

heatmap<-heatmap.2(all_sig_RS,
          dendrogram="row",
          trace="none",
          margins = c(2,1),
          cexRow = 0.1,
          cexCol = 1,
          Colv = "NA",)
names<-rownames(all_sig_RS)[heatmap$rowInd]
write.csv(names, "DEGs_filtered_gene_order_RS.csv", row.names=F)

