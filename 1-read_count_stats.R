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
allowWGCNAThreads()

# where are kalliso files?

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files (already made so load txi object)
load("~/Documents/S_L_DH/samples/tximpost_all.RData")
# contains the kallisto output directories - 1 per sample.
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)
colnames(colData)<-c("Sample", "Rep", "Treatment", "Timepoint", "Genotype", "Factor", "Condition", "Resistant")
# ==================================================================================


# check that order of samples in metadata and txi object are the same

# ==================================================================================
# read count stats chunk starts here

expressed_genes<-txi.kallisto.tsv$abundance
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
# expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),] #only use if you dont want LC genes
expressed_genes<-expressed_genes[,c(97, 1:96)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:97)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Sample")
sub<-all_wheat_genes[,c(8, 2, 3, 4)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>=0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>=0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>=0.5, 1,0)
rep_wise$Sum<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3

expressed<-rep_wise[(rep_wise$Sum >=2),]

for(i in unique(expressed$Factor)){
  data<-expressed[(expressed$Factor==i),]
  factor<-paste(i)
  write.csv(data, file=paste("~/Documents/S_L_DH/DEtesting/data/", factor, ".csv", sep=""))
  assign(factor, data)
}
{
G1<-rbind(C1G, T1G)
G1<-G1[!(duplicated(G1$GeneID)),]
G2<-rbind(C2G, T2G)
G2<-G2[!(duplicated(G2$GeneID)),]
G3<-rbind(C3G, T3G)
G3<-G3[!(duplicated(G3$GeneID)),]
G4<-rbind(C4G, T4G)
G4<-G4[!(duplicated(G4$GeneID)),]

L1<-rbind(C1L, T1L)
L2<-rbind(C2L, T2L)
L3<-rbind(C3L, T3L)
L4<-rbind(C4L, T4L)
L1<-L1[!(duplicated(L1$GeneID)),]
L2<-L2[!(duplicated(L2$GeneID)),]
L3<-L3[!(duplicated(L3$GeneID)),]
L4<-L4[!(duplicated(L4$GeneID)),]

R1<-rbind(C1R, T1R)
R2<-rbind(C2R, T2R)
R3<-rbind(C3R, T3R)
R4<-rbind(C4R, T4R)
R1<-R1[!(duplicated(R1$GeneID)),]
R2<-R2[!(duplicated(R2$GeneID)),]
R3<-R3[!(duplicated(R3$GeneID)),]
R4<-R4[!(duplicated(R4$GeneID)),]

S1<-rbind(C1S, T1S)
S2<-rbind(C2S, T2S)
S3<-rbind(C3S, T3S)
S4<-rbind(C4S, T4S)
S1<-S1[!(duplicated(S1$GeneID)),]
S2<-S2[!(duplicated(S2$GeneID)),]
S3<-S3[!(duplicated(S3$GeneID)),]
S4<-S4[!(duplicated(S4$GeneID)),]

G1$Comparison<-"G1"
G2$Comparison<-"G2"
G3$Comparison<-"G3"
G4$Comparison<-"G4"
L1$Comparison<-"L1"
L2$Comparison<-"L2"
L3$Comparison<-"L3"
L4$Comparison<-"L4"
R1$Comparison<-"R1"
R2$Comparison<-"R2"
R3$Comparison<-"R3"
R4$Comparison<-"R4"
S1$Comparison<-"S1"
S2$Comparison<-"S2"
S3$Comparison<-"S3"
S4$Comparison<-"S4"

G1<-G1[,c(2, 10)]
G2<-G2[,c(2, 10)]
G3<-G3[,c(2, 10)]
G4<-G4[,c(2, 10)]
L1<-L1[,c(2, 10)]
L2<-L2[,c(2, 10)]
L3<-L3[,c(2, 10)]
L4<-L4[,c(2, 10)]
R1<-R1[,c(2, 10)]
R2<-R2[,c(2, 10)]
R3<-R3[,c(2, 10)]
R4<-R4[,c(2, 10)]
S1<-S1[,c(2, 10)]
S2<-S2[,c(2, 10)]
S3<-S3[,c(2, 10)]
S4<-S4[,c(2, 10)]}

all_filtered_lists<-rbind(G1,
                          G2,
                          G3,
                          G4,
                          L1,
                          L2,
                          L3,
                          L4,
                          R1,
                          R2,
                          R3,
                          R4,
                          S1,
                          S2,
                          S3,
                          S4)

write.csv(all_filtered_lists, file="~/Documents/S_L_DH/DEtesting/data/all_lists_filtered.csv", row.names = F)
# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Genotype, shape=Treatment), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black"))

vst_counts<-assay(vsd)

