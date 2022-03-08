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
rm(binary)
rm(binary_f)
rm(binary.df)
rm(binary.df_f)

alltraits = binary.df_f[keepSamples,]

rm(trait_f)
rm(trait_f.df)
save(alltraits, file="femaletraitafterclust.Rdata")

getwd()
save(pheno.sub_m, file="pheno.sub_m.Rdata")
