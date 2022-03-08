###########################################
#  -Integrative bioinformatics            # 
#  - RNA-Seq                              #
#  - 2018- 12- 18                         #
#  - Copyright: Mohamed Hamed             #
###########################################
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
# 
# BiocManager::install("vsn")

library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)
library("dplyr")




##### load the mirna-Seq data #####
data.path="D:/mirna/validated M&F miRNA/gdc_download_20191121_192643.372334"
files <- list.files(path=data.path,recursive=T, pattern = "mirnas.quantification.txt")

# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]

#open a connection to your gz file and read the file
temp <- read.table(file.path(data.path,files[1]), header=T)



#create a storing object mirna.exp to save the whole read counts of each file read in an iteration
mirna.exp=cbind(temp[1],temp[2])

rownames(mirna.exp)=mirna.exp[,1]
mirna.exp=mirna.exp[-1]
colnames(mirna.exp)=c(file.id)

for(i in 2: length(files))
{
  
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  temp <- read.table(file.path(data.path,files[i]), header=T)
  
  
  ## remove the first column, bec we had it already
  temp=temp[2]
  colnames(temp)=c(file.id)
  
  mirna.exp=cbind(mirna.exp,temp)
}


pheno <-  read_delim("D:/mirna/validated M&F miRNA/gdc_sample_sheet.2019-11-21.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

clin_1 <- read_delim("D:/mirna/validated M&F miRNA/clinical.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
class(clin_1)
clin_1_df= as.data.frame(clin_1)

a = clin_1_df[,2]
b= clin_1_df[,13]


clinical = cbind(a,b)
View(clinical)

clinical.df = as.data.frame(clinical)
names(clinical.df)=c("Case ID","Gender")

cl=clinical.df
p = pheno
merged_pheno=merge(p,cl,by="Case ID")


dim(merged_pheno)


dupl.df=duplicated.data.frame(merged_pheno)
sum(dupl.df)

copy_merged= merged_pheno
no_dupl= distinct(copy_merged)
table(no_dupl$`Sample Type`,no_dupl$Gender)

getwd()
setwd("D:/mirna/validated M&F miRNA")
write.table(no_dupl, file = "updated_pheno_miRNA.txt")

pheno_updated = no_dupl
View(pheno_updated)


rm(cl)
rm(clin_1)
rm(clin_1_df)
rm(clinical)
rm(clinical.df)
rm(copy_merged)
rm(p)
rm(no_dupl)

rm(merged_pheno)
###########################################################


# rename column names : replace the spaces with dots 
pheno.names=names(pheno_updated)
names(pheno_updated)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno_updated$Sample.Type,pheno_updated$Gender)


#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 

file.ids=colnames(mirna.exp)

file.ids.pheno=pheno_updated$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mirna.exp)=pheno_updated$Sample.ID[index.files]


n_f = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("female"))
t_f = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("female"))


all_n_f_samples = pheno_updated[n_m,]$`Sample.ID`
all_t_f_samples = pheno_updated[t_m,]$`Sample.ID`


# now we will retrieve  the exp data for these 20 samples only and also the pheno data
normal.exp_f=mirna.exp[, names(mirna.exp)%in% all_n_f_samples]
tumor.exp_f =mirna.exp[, names(mirna.exp)%in% all_t_f_samples]


pheno.sub_f=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_f_samples,all_t_f_samples), c("Sample.ID", "Sample.Type")]


exp.sub_f=cbind(normal.exp_f,tumor.exp_f)
exp.sub_f=apply (exp.sub_f, 2,as.integer)

rownames(exp.sub_f)=rownames(normal.exp_f)


save(mirna.exp,pheno_updated, exp.sub_f,pheno.sub_f ,file="microRNA-seq.RDATA")


###### DO the differential EXP analysis using DeSeq2
cond1f="Solid Tissue Normal"
cond2f="Primary Tumor"

dim(exp.sub_f)
dim(pheno.sub_f)
dds_f = DESeqDataSetFromMatrix( countData = exp.sub_f , colData = pheno.sub_f , design = ~ Sample.Type)
dds.run_f = DESeq(dds_f)

### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)

res_f=results(dds.run_f)
res_f=results(dds.run_f, contrast = c("Sample.Type",cond1f ,cond2f) )


# remove nulls
res_f=res_f[complete.cases(res_f), ]
summary(res_f)

res.df_f=as.data.frame(res_f)

write.table(res.df_f, file = "res_f_mirna.txt")


plotMA(res_f, ylim=c(-1,1))

res.degs_f=res.df_f[res.df_f$padj< 0.05 & abs(res.df_f$log2FoldChange)>log2(2),]
library(xlsx)
write.table(res.degs_f, file = "res_degs_f.txt")


#for Tfmir
res.degs_f$regulation=0
res.degs_f=as.data.frame(res.degs_f)
res.degs_f[res.degs_f$log2FoldChange<0,]$regulation =-1
res.degs_f[res.degs_f$log2FoldChange>0,]$regulation =1
res.degs2_f= cbind(rownames(res.degs_f), res.degs_f$regulation)
write.table(res.degs2_f, file("res_reg_f_mirna.tsv"), quote = F, col.names = F, row.names = F)


#expression of these degs for plots
exp.degs_f= exp.sub_f[rownames(exp.sub_f) %in% rownames(res.degs_f), ]
write.table(exp.degs_f, file = "exp_dems_f_mirna.txt")

#######################################
#### get the normalized and loggedtransformed values of all exp data
#using the the variance stabilizing transformation. vsn package
ntd_f=normTransform(dds_f)
exp.sub.norm_f= assay(ntd_f)




# Make a basic volcano plot
par(mfrow=c(1,1))
with(res_f, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_f, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_f, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



########################################################################
##2) second comparison between male and control : now we will retrieve  the exp data for these
#male samples only and also the pheno data



n_m = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("male"))
t_m = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("male"))


all_n__msamples = pheno_updated[n_m,]$`Sample.ID`
all_t_m_samples = pheno_updated[t_m,]$`Sample.ID`


# now we will retrieve  the exp data for these 20 samples only and also the pheno data
normal.exp_m=mirna.exp[, names(mirna.exp)%in% all_n_m_samples]
tumor.exp_m =mirna.exp[, names(mirna.exp)%in% all_t_m_samples]


pheno.sub_m=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_m_samples,all_t_m_samples), c("Sample.ID", "Sample.Type")]


exp.sub_m=cbind(normal.exp_m,tumor.exp_m)
exp.sub_m=apply (exp.sub_m, 2,as.integer)

rownames(exp.sub_m)=rownames(normal.exp_m)


save(mirna.exp,pheno_updated, exp.sub_m,pheno.sub_m ,file="microRNA-seq.RDATA")


###### DO the differential EXP analysis using DeSeq2
cond1m="Solid Tissue Normal"
cond2m="Primary Tumor"

dim(exp.sub_m)
dim(pheno.sub_m)
dds_m = DESeqDataSetFromMatrix( countData = exp.sub_m , colData = pheno.sub_m , design = ~ Sample.Type)
dds.run_m = DESeq(dds_m)

### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)

res_m=results(dds.run_m)
res_m=results(dds.run_m, contrast = c("Sample.Type",cond1m ,cond2m) )


# remove nulls
res_m=res_m[complete.cases(res_m), ]
summary(res_m)

res.df_m=as.data.frame(res_m)

write.table(res.df_m, file = "res.txt")


plotMA(res_m, ylim=c(-1,1))
summary (res_m)

res.degs_m=res.df_m[res.df_m$padj< 0.05 & abs(res.df_m$log2FoldChange)>log2(2),]
library(xlsx)
write.table(res.degs_m, file = "res_degs_m.xlsx")


#for Tfmir
res.degs_m$regulation=0
res.degs_m=as.data.frame(res.degs_f)
res.degs_m[res.degs_m$log2FoldChange<0,]$regulation =-1
res.degs_m[res.degs_m$log2FoldChange>0,]$regulation =1
res.degs2_m= cbind(rownames(res.degs_m), res.degs_m$regulation)
write.table(res.degs2_m, file("res_reg_m.tsv"), quote = F, col.names = F, row.names = F)


#expression of these degs for plots
exp.degs_m= exp.sub_m[rownames(exp.sub_m) %in% rownames(res.degs_m), ]
write.table(exp.degs_m, file = "exp_dems_m.xlsx")

#######################################
#### get the normalized and loggedtransformed values of all exp data
#using the the variance stabilizing transformation. vsn package
ntd_m=normTransform(dds_m)
exp.sub.norm_m= assay(ntd_m)




# Make a basic volcano plot
par(mfrow=c(1,1))
with(res_m, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_m, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_m, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


