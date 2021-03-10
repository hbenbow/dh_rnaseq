library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
allowWGCNAThreads()

# Import kallisto files with txi 
samples<-as.data.frame(dir("samples/"))
files <- file.path(samples, "abundance.h5")
names(files) <- paste0(samples)
all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)

# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor)

# filter out low abundance transcripts
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= .5 ) >= 64
dds <- dds[idx,]
# This would say, e.g. filter out genes where there are fewer 
# than 32 samples with normalized counts greater than or equal to 5. Chose 32 as it is 2/3 of treated samples.

# variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst_counts<-assay(vsd)
# expression data as df, samples as rows, genes as columns
dataExp0<-as.data.frame(t(vst_counts))

pca_data<-plotPCA(vsd, intgroup=c("Genotype", "Treatment"), returnData=T)
ggplot(pca_data, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Genotype), size=4)

# check for genes and samples with too many missing values
gsg = goodSamplesGenes(dataExp0, verbose = 3)
gsg$allOK

# remove bad genes
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(dataExp0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(dataExp0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExp = dataExp0[gsg$goodSamples, gsg$goodGenes]
}
if (gsg$allOK){
  datExp = dataExp0
}

# check for outliers
sampleTree = hclust(dist(datExp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.1);
par(mar = c(4,6,4,4))
plot(sampleTree, main = "Sample clustering of expression patterns", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 240, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 240, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1 | clust==1)
datExpr2 = datExp[keepSamples, ]
nGenes = ncol(datExpr2)
nSamples = nrow(datExpr2)

# Choose a set of soft-thresholding powers
powers = c(1:10)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr2, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



softPower =4;
allowWGCNAThreads()
bwnet = blockwiseModules(datExpr2, maxBlockSize = 20000, power = softPower, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE,saveTOMs = TRUE, saveTOMFileBase = "S_L_nobulks-blockwise", verbose = 3)

 bsizeGrWindow(6,6)
# plot<-
  plotDendroAndColors(bwnet$dendrograms[[1]], moduleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], moduleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(bwnet$dendrograms[[3]], moduleColors[bwnet$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(bwnet$dendrograms[[4]], moduleColors[bwnet$blockGenes[[4]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 4",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

all<-cbind(genes, moduleColors)
for(i in paste(unique(all$moduleColor))){
  id<-paste(i)
  data<-all[(all$moduleColor == id),]
  list<-data[,1]
  colnames(list)<-NULL
  write.table(list, file=paste0(id, "exp.txt"), row.names=F, col.names=F)
  assign(id, data)
}

moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];

# Calculate eigengenes
MEList = moduleEigengenes(datExpr2, colors = moduleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr2, moduleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "module_colours_and_labels_mixed.RData")



sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="cyan";
plotMat(t(scale(datExpr2[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )


hubs<-as.data.frame(chooseTopHubInEachModule(
  datExpr2, 
  moduleColors, 
  omitColors = "grey", 
  power = 6, 
  type = "unsigned"))


# Define numbers of genes and samples
nGenes = ncol(datExpr2);
nSamples = nrow(datExpr2)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr2, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(traitData$Resistance);
names(weight) = "Resistant"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# genes with high significance to trait
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Genotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr2)
names(datExpr2)[moduleColors=="blue"]
geneInfo0 = data.frame(GeneID = names(datExpr2),
                       moduleColor = moduleColors)
modOrder = order(-abs(cor(MEs, weight, use = "p")));
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$MM.blue));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")

library(WGCNA)
allowWGCNAThreads()



# export to cytoscape
blockNumbers<-as.data.frame(cbind(colnames(datExpr2), bwnet$blocks))
blockColours<-geneInfo[,c(1,2)]
colnames(blockNumbers)<-c("GeneID", "Block")
colnames(blockColours)<-c("GeneID", "Colour")
Blocks<-merge(blockNumbers, blockColours, by="GeneID")
# for (i in paste(unique(Blocks$Colour))){
modules = "cyan"
block<-unique(as.character(Blocks[(Blocks$Colour == modules),][,2]))
genes<-Blocks[(Blocks$Block %in% block),]
probes<-genes$GeneID
datExp3<-datExpr2
datExp3<-datExp3[,probes]
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExp3, power = 9.9);
# Select module probes
colours<-as.factor(genes$Colour)
inModule = is.finite(match(colours, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = NULL,
                               nodeAttr = colours[inModule])
# }




time_hubs<-c("TraesCS7D02G347300.1",
"TraesCS6B02G036400.1",
"TraesCS4D02G295900.1",
"TraesCS7B02G401800.2",
"TraesCS7B02G418800.1")
