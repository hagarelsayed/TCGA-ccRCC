#====================================================================================
#
# ccRCC project                     Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory.  On Windows use a forward slash / instead of the usual \.
workingDir = "D:/ccRCC project/wcgna project/m_mRNA/output";
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
ccrcc_maleData = exp.sub_m
# Take a quick look at what is in the data set:
dim(ccrcc_maleData);
names(ccrcc_maleData);


datExpr0_m = as.data.frame(t(ccrcc_maleData[,]));

names(datExpr0_m) = rownames(ccrcc_maleData);
rownames(datExpr0_m) = names(ccrcc_maleData);

# repeat this step over till gsg$allok = true
gsg = goodSamplesGenes(datExpr0_m, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0_m)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0_m)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0_m = datExpr0_m[gsg$goodSamples, gsg$goodGenes]
}

#end of step


sampleTree = flashClust(dist(datExpr0_m), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = 4e+06, col = "red");
#abline(h = 4.25e+06, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight =  4e+06, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0_m[keepSamples, ]

sum(duplicated.data.frame(datExpr))
nGenes = ncol(datExpr_mm)
nSamples = nrow(datExpr_mm)
rm(datExpr_mm)
# dismiss = (clust==2)
# dismiss_expr=datExpr0_m[dismiss, ]
# rownames(dismiss_expr)
# rm(dismiss)
# rm(dismiss_expr)


save(datExpr_mm, file = "datExpr_m.Rdata")
save(exp.sub_m, file = "exp_sub_m.Rdata")

################################################################################

trait=pheno.sub_m
trait.df=as.data.frame(trait)
binary= within(trait_m.df, Sample.Type <- factor(Sample.Type, labels = c(1, 0)))
binary.df = as.data.frame(binary)

rm(trait)
rm(binary)
rm(binary_m)
rm(binary.df)
rm(binary.df_f)


alltraits = binary.df[keepSamples,]

save(alltraits, file= "alltraits.Rdata")

datTraits = a

library(dplyr)
dupl.df=duplicated.data.frame(d)
View(dupl.df)
no_dupl= distinct(d)

alltraits = no_dupl

table(no_dupl$Sample.ID,no_dupl$Sample.Type)

sum(dupl.df)
rm(femaleSamples)

e=alltraits
rownames(e)=e$Sample.ID
e = e[-1]
datTraits

save(datTraits, file = "datTrait_f.Rdata")

collectGarbage();







