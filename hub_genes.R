hub_genes<-c("TraesCS2D02G573200.1",
"TraesCS3A02G044400.1",
"TraesCS3A02G518500.1",
"TraesCS5B02G221600.1",
"TraesCS5B02G291900.1",
"TraesCS7B02G277000.1",
"TraesCS7B02G338900.1",
"TraesCS7D02G302000.1")

hub_exp<-subset(all_filtered, all_filtered$row %in% hub_genes)
hub_tpm<-tpm[hub_genes,]

ggplot(hub_exp, aes(x=as.factor(Timepoint), y=log2FoldChange)) + 
  geom_col(position=position_dodge(, preserve="single")) +
  facet_grid(Genotype~row)
write.csv(hub_exp, file="~/Documents/S_L_DH/data/hub_genes_expression.csv")
write.csv(hub_tpm, file="~/Documents/S_L_DH/data/hub_genes_tpm.csv")
