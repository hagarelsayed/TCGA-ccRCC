#====================================================================================
#
# ccRCC project                     Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory.  On Windows use a forward slash / instead of the usual \.
workingDir = "D:/ccRCC project/wcgna project/f_mRNA/Output";
setwd(workingDir);
# Load WGCNA package
library(WGCNA)
library(flashClust)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#Read in the ccrcc data set exp for female tumor and normal 
############################################################

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to

file.ids.pheno=pheno_updated$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno_updated$Sample.ID[index.files]

n_f = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("female"))
t_f = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("female"))
n_m = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("male"))
t_m = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("male"))
View(n_f)
View(all_n_f_samples)
all_n_f_samples = pheno_updated[n_f,]$`Sample.ID`
all_t_f_samples = pheno_updated[t_f,]$`Sample.ID`
all_n_m_samples = pheno_updated[n_m,]$`Sample.ID`
all_t_m_samples = pheno_updated[t_m,]$`Sample.ID`


normal.exp_f = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_n_f_samples]
tumor.exp_f = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_t_f_samples]
normal.exp_m = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_n_m_samples]
tumor.exp_m = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_t_m_samples]

pheno.sub_f=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_f_samples,all_t_f_samples), c("Sample.ID", "Sample.Type")]
pheno.sub_m=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_m_samples,all_t_m_samples), c("Sample.ID", "Sample.Type")]

exp.sub_f=cbind(normal.exp_f,tumor.exp_f)
exp.sub_m=cbind(normal.exp_m,tumor.exp_m)

#prpring the traits data
res.degs_m$regulation=0
res.degs_m=as.data.frame(res.degs_m)
res.degs_m[res.degs_m$log2FoldChange<0,]$regulation =-1
res.degs_m[res.degs_m$log2FoldChange>0,]$regulation =1
res.degs2_m= cbind(rownames(res.degs_m), res.degs_m$regulation)
write.table(res.degs2_m, file("res_reg_m.tsv"))

trait_f_mrna= pheno.sub_f
trait_f_mrna = as.data.frame(trait_f_mrna)

trait_f_mrna$binary =0
trait_f_mrna = as.data.frame(trait_f_mrna)
x=trait_f_mrna
rownames(x)=x[,1]
x=x[-1]


# trait_f_mrna[trai$log2FoldChange<0,]$binary = 0
# res.degs_m[res.degs_m$log2FoldChange<0,]$binary=1
# pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")
# res.degs_m[res.degs_m$log2FoldChange<0,]$regulation =-1
# 
# trait_f_mrna[trait_f_mrna$Sample.ID %in% c("Solid Tissue Normal"),]$binary= 0
# trait_f_mrna[trait_f_mrna$Sample.ID %in% c("Primary Tumor"),]$binary= 1
rm(trait_f_mrna)

#control=0
#tumor = 1

trait$sample

library(readr)
getwd()


#Do the analysis
cond1f="Solid Tissue Normal"
cond2f="Primary Tumor"

dim(exp.sub_f)
dim(pheno.sub_f)
#====================================================================================
#
# The previous part was for prpring the data 
#
#=====================================================================================

#Read in the female liver data set
ccrcc_femData = exp.sub_f
# Take a quick look at what is in the data set:
dim(ccrcc_femData);
names(ccrcc_femData);


datExpr0 = as.data.frame(t(ccrcc_femData[,]));

names(datExpr0) = rownames(ccrcc_femData);
rownames(datExpr0) = names(ccrcc_femData);

# repeat this step over till gsg$allok = true
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#end of step


sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = 3700000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 3700000, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
# dismiss = (clust==0)
# dismiss_expr=datExpr0[dismiss, ]
#rownames(dismiss_expr)
# rm(dismiss)
# rm(dismiss_expr)


datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



save(exp.sub_f, file = "exp.sub_f.Rdata")
save(exp.sub_m, file = "exp_sub_m.Rdata")

save(pheno.sub_m, file= "pheno.sub.Rdata")

save(datExpr, file="datExpr_f.Rdata")
########################################################################
#we need to make the trait data binary

# x_clust= x[keepSamples, ]

trait_f=pheno.sub_f
trait_f.df=as.data.frame(trait_f)
binary_f= within(trait_f.df, Sample.Type <- factor(Sample.Type, labels = c(1, 0)))
binary.df_f = as.data.frame(binary_f)
dim(binary.df_f)
alltraits = binary.df_f[keepSamples,]

rm(binary)
rm(binary_f)
rm(binary.df)
rm(binary.df_f)


rm(trait_f)
rm(trait_f.df)
save(alltraits, file="all traits female 204.Rdata")


save(pheno.sub_m, file="pheno.sub_m.Rdata")


a=alltraits
rownames(a)=a$Sample.ID
a = a[-1]

datTraits = a
rm(a)
rm(res.degs)
rm(dismiss_expr)
rm(dismiss)
save(datTraits, file = "datTrait_f.Rdata")
class(datTraits)
datTraits_f= as.numeric(datTraits)

collectGarbage();

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(as.numeric(unlist(datTraits)), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "Femaleccrcc-01-dataInput.RData")



# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Femaleccrcc-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
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


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "FemaleLiver-02-networkConstruction-auto.RData")
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Femaleccrcc-02-networkConstruction-auto.RData")


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.

setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames
getwd()
#The variable lnames contains the names of loaded variables.

# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
SampleType = as.data.frame(datTraits$Sample.Type);
names(SampleType) = "Sample Type"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, SampleType, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(SampleType), sep="");
names(GSPvalue) = paste("p.GS.", names(SampleType), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for SampleType",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr)[moduleColors=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")


####triall

agregated = datExpr 
agg = as.data.frame(t(agregated[,]));

names(datExpr0) = rownames(ccrcc_femData);
rownames(datExpr0) = names(ccrcc_femData);