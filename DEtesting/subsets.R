# 
# subsets
# GU RU
# GD RD
# LU SU
# LD SD


for(genotype in unique(all_filtered$Genotype)){
  data<-all_filtered[(all_filtered$Genotype==genotype),]
  assign(paste(genotype), data)
}

GU<-Stigg[(Stigg$log2FoldChange>0),]
GD<-Stigg[(Stigg$log2FoldChange<0),]
LU<-Longbow[(Longbow$log2FoldChange>0),]
LD<-Longbow[(Longbow$log2FoldChange<0),]
RU<-Resistant[(Resistant$log2FoldChange>0),]
RD<-Resistant[(Resistant$log2FoldChange<0),]
SU<-Susceptible[(Susceptible$log2FoldChange>0),]
SD<-Susceptible[(Susceptible$log2FoldChange<0),]

GU_RU<-subset(GU, GU$row %in% RU$row)[,1]
GD_RD<-subset(GD, GD$row %in% RD$row)[,1]
LU_SU<-subset(LU, LU$row %in% SU$row)[,1]
LD_SD<-subset(LD, LD$row %in% SD$row)[,1]

GU_RU<-subset(all_filtered, all_filtered$row %in% GU_RU)
GD_RD<-subset(all_filtered, all_filtered$row %in% GD_RD)
LU_SU<-subset(all_filtered, all_filtered$row %in% LU_SU)
LD_SD<-subset(all_filtered, all_filtered$row %in% LD_SD)

write.csv(GU_RU, "GU_RU.csv")
write.csv(GD_RD, "GD_RD.csv")
write.csv(LU_SU, "LU_SU.csv")
write.csv(LD_SD, "LD_SD.csv")
