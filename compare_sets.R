# Compare analyses

sep_no_filter<-read.csv("~/Documents/S_L_DH/data/DE_tests/separate_no_filter.csv")
sep_filter<-read.csv("~/Documents/S_L_DH/data/DE_tests/separate_filtered.csv")
together_no_filter<-read.csv("~/Documents/S_L_DH/data/DEGs_all_data.csv")
together_filter<-read.csv("~/Documents/S_L_DH/data/DE_tests/together_no_filter.csv")

gene_list1<-as.data.frame(sep_no_filter$row)
gene_list2<-as.data.frame(sep_filter$row)
gene_list3<-as.data.frame(together_filter$row)
gene_list4<-as.data.frame(together_no_filter$row)
gene_list1$Dataset<-1
gene_list2$Dataset<-1
gene_list3$Dataset<-1
gene_list4$Dataset<-1

colnames(gene_list1)<-c("Gene", "Dataset")
colnames(gene_list2)<-c("Gene", "Dataset")
colnames(gene_list3)<-c("Gene", "Dataset")
colnames(gene_list4)<-c("Gene", "Dataset")

gene_list1<-gene_list1[!duplicated(gene_list1$Gene),]
gene_list2<-gene_list2[!duplicated(gene_list2$Gene),]
gene_list3<-gene_list3[!duplicated(gene_list3$Gene),]
gene_list4<-gene_list4[!duplicated(gene_list4$Gene),]

gene_list<-rbind(gene_list1, gene_list2, gene_list3, gene_list4)
gene_list<-as.data.frame(gene_list[!duplicated(gene_list$Gene),1])
colnames(gene_list)<-"Gene"

joined<-left_join(gene_list, gene_list1, by="Gene")
joined<-left_join(joined, gene_list2, by="Gene")
joined<-left_join(joined, gene_list3, by="Gene")
joined<-left_join(joined, gene_list4, by="Gene")
joined[is.na(joined)] <- 0

colnames(joined)<-c("Gene", "sep_no_filter",
                    "sep_filter",
                    "Together_filter",
                    "Together_no_filter")
rownames(joined)<-joined$Gene
joined$Gene<-NULL
joined<-as.matrix(joined)
