S_L_genes<-as.data.frame(c("TraesCS3B02G113500", 
                           "TraesCS4D02G015000", 
                           "TraesCS5B02G040500", 
                           "TraesCS5D02G259800LC", 
                           "TraesCS6A02G445400LC", 
                           "TraesCS7A02G488100", 
                           "TraesCS7B02G454000", 
                           "TraesCS7D02G535400LC", 
                           "TraesCSU02G452300LC", 
                           "TraesCS1D02G341000", 
                           "TraesCS3B02G238800", 
                           "TraesCS3D02G376500LC", 
                           "TraesCS4A02G256700LC", 
                           "TraesCS5A02G531900LC", 
                           "TraesCS5A02G178200", 
                           "TraesCSU02G589600LC", 
                           "TraesCS1A02G346100LC", 
                           "TraesCS5A02G555900LC", 
                           "TraesCS5B02G220300LC", 
                           "TraesCS7D02G144400", 
                           "TraesCS1A02G592500LC", 
                           "TraesCS6B02G629700LC ", 
                           "TraesCS6B02G284200", 
                           "TraesCS5D02G419100LC"))

S_L_genes$dataset<-"S_L"
colnames(S_L_genes)<-c("GeneID", "Dataset")


exp_data<-read.csv("~/Documents/S_L_DH/Limagrain/DEGs_all_data_qtls.csv")
exp_data2<-merge(exp_data, wheat_bed, by.x="row", by.y="GeneID")
write.csv(qtls_clusters, file="~/Documents/S_L_DH/Limagrain/STB_DEGs_QTLs.csv")

candidates<-subset(data,data$row %in% S_L_genes$GeneID)
candidates_for_exp<-as.character(unique(candidates$row))
candidates_for_exp<-subset(exp_data, exp_data$row %in% candidates_for_exp)

ggplot(candidates_for_exp, aes(x=as.factor(Timepoint), y=log2FoldChange, group=Genotype)) +
  geom_line(aes(colour=Genotype),size=1.5) + 
  facet_grid(row~.) + 
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange + 
                      lfcSE, colour=Genotype), width=0.3)+
  scale_colour_manual(values=c("springgreen4","cadetblue4", "chocolate3", "wheat2"))


ggplot(candidates_for_exp, aes(x=as.factor(Timepoint), y=log2FoldChange, group=Genotype)) +
  geom_col(aes(fill=Genotype), position=position_dodge(preserve="single")) + 
  facet_wrap(row~., ncol=1) + 
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange + 
                      lfcSE, group=Genotype),position=position_dodge(preserve="single", width=0.9), width=0.5) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"), axis.text.x = element_text(colour="black"))+
  xlab("Timepoint (dpi)") +
  ylab("Log(2) Fold Change")+
  scale_fill_manual(values=c("forestgreen","palegreen1", "chocolate3", "wheat2")) +
  geom_hline(aes(yintercept=1), alpha=0.5)+
  geom_hline(aes(yintercept=-1), alpha=0.5) +
  geom_hline(aes(yintercept=0))


test test

