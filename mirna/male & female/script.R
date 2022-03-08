install.packages("dplyr")
library("dplyr")

library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)



clin_1 <- read_delim("D:/mirna/clinical.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
class(clin_1)

x= as.data.frame(clin_1)
View(x)
a = x[,2]
as.data.frame(clin_1)
View(a)
b= x[,13]
View(b)
class(b)
clinical = cbind(a,b)
View(clinical)


pheno <- read_delim("D:/mirna/only mirna no isoform/gdc_sample_sheet.2019-11-17.tsv", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)
View(pheno )
clinical.df = as.data.frame(clinical)

 clinical.df = as.data.frame(clinical)
 names(clinical.df)=c("Case ID","Gender")
 View(clinical.df)
 case_ids = clinical.df$`Case ID`
 case_ids_clin = pheno$`Case ID`
 index_file= match(case_ids,case_ids_clin)
 s=clinical.df
 y = pheno
 z = y$`Sample Type`[index_file] 
#g = s$Gender[index_file]
all= cbind(s,z)
View(all)
 View(y)
?merge
merged=merge(y,s,by="Case ID")

merged=merge(y,s,by="Case ID")
 View(merged)
 dim(merged)
#[1] 1232    9
 dupl=duplicated(merged$`Case ID`)  
 sum(dupl)
#[1] 716
 
 table(merged$`Sample Type`,merged$Gender)
# 
#                             female male
# Additional - New Primary      0    2
 #Primary Tumor               376  712
 #Solid Tissue Normal          38  104

 dupl=duplicated(merged$`Case ID`)  


install.packages("dplyr")
library("dplyr")
dupl.df=duplicated.data.frame(merged)
sum(dupl.df)
# [1] 616
copy_merged= merged
no_dupl= distinct(copy_merged)
table(no_dupl$`Sample Type`,no_dupl$Gender)
# 
#                           female male
# Additional - New Primary      0    1
# Primary Tumor               188  356
# Solid Tissue Normal          19   52
getwd()
setwd("D:/mirna/male & female")
write.table(no_dupl, file = "pheno_w_gender.txt")
getwd()
pheno_updated = no_dupl
View(pheno_updated)

#all.normal.samples= pheno[ pheno$Sample.Type %in% c("Solid Tissue Normal"),]$Sample.ID
#which( data$V1 > 2 | data$V2 < 4)  # this means orr not and
#my.data.frame <- data[(data$V1 > 2) & (data$V2 < 4), ]
n_f = which(pheno_updated$`Sample Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("female"))
t_f = which(pheno_updated$`Sample Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("female"))
n_m = which(pheno_updated$`Sample Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("male"))
t_m = which(pheno_updated$`Sample Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("male"))
all_n_f_samples=pheno_updated[n_f,]$`Sample ID`
View(all_n_f_samples)
all_t_f_samples = pheno_updated[t_f,]$`Sample ID`
all_n_m_samples = pheno_updated[n_m,]$`Sample ID`
all_t_m_samples = pheno_updated[t_m,]$`Sample ID`

rm(a)
rm(b)
rm(case_ids)
rm(case_ids_clin)
rm(dupl)
rm(dupl.df)
rm(lol)
rm(z)

###########################################################


##### load the mirna-Seq data #####

data.path ="D:/mirna/only mirna no isoform/gdc_download_20191117_081915.282629"
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

#pheno here is pheno_updated
table(pheno_updated$`Sample Type`,pheno_updated$Gender)

pheno.names=names(pheno_updated)
names(pheno_updated)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))

file.ids=colnames(mirna.exp)

file.ids.pheno=pheno_updated$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mirna.exp)=pheno_updated$Sample.ID[index.files]

# 1) first comparison between Female and control : now we will retrieve  the exp data for these female samples only and also the pheno data

normal.exp_f =mirna.exp[, names(mirna.exp)%in% all_n_f_samples]
tumor.exp_f =mirna.exp[, names(mirna.exp)%in% all_t_f_samples]


pheno_sub_f=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_f_samples,all_t_f_samples), c("Sample.ID", "Sample.Type")]
rm(pheno_sub)
rm(s)
rm(x)
rm(y)

exp_sub_f=cbind(normal.exp_f,tumor.exp_f)
rm(exp_sub)
exp_sub_f=apply (exp_sub_f, 2,as.integer)
rownames(exp_sub_f)=rownames(normal.exp_f)

getwd()
setwd("D:/mirna/male & female/female")
save(mirna.exp,pheno_updated, exp_sub_f,pheno_sub_f ,file="microRNA-seq.RDATA")


###### DO the differential EXP analysis using DeSeq2
cond1f="Solid Tissue Normal" 
cond2f="Primary Tumor"

dim(exp_sub_f)
dim(pheno_sub_f)

dds_f = DESeqDataSetFromMatrix( countData = exp_sub_f , colData = pheno_sub_f , design = ~ Sample.Type)
dds_run_f = DESeq(dds_f)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res_f=results(dds_run_f)
res_f=results(dds_run_f, contrast = c("Sample.Type",cond1f ,cond2f) )

# remove nulls
res_f=res_f[complete.cases(res_f), ]
summary(res_f)


res_f.df=as.data.frame(res_f)
getwd()
write.table(res_f.df, file = "res.txt")

plotMA(res, ylim=c(-1,1)) 
summary (res)

res_f.degs=res_f.df[res_f.df$padj< 0.05 & abs(res_f.df$log2FoldChange)>log2(2),]

write.table(res_f.degs, file = "res_degs.txt")

#for Tfmir
res_f.degs$regulation=0
res_f.degs=as.data.frame(res_f.degs)
res_f.degs[res_f.degs$log2FoldChange<0,]$regulation =-1
res_f.degs[res_f.degs$log2FoldChange>0,]$regulation =1
res_f.degs2= cbind(rownames(res_f.degs), res_f.degs$regulation)
write.table(res_f.degs2, file="res_f_reg.tsv", quote = F, col.names = F, row.names = F)

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




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# 2) second comparison between male and control : now we will retrieve  the exp data for these male samples only and also the pheno data
normal.exp_m =mirna.exp[, names(mirna.exp)%in% all_n_m_samples]
View(normal.exp_f)
View(normal.exp_m)
tumor.exp_m =mirna.exp[, names(mirna.exp)%in% all_t_m_samples]


normal.exp_m =mirna.exp[, names(mirna.exp)%in% all_n_m_samples]
tumor.exp_m =mirna.exp[, names(mirna.exp)%in% all_t_m_samples]


pheno_sub_m=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_m_samples,all_t_m_samples), c("Sample.ID", "Sample.Type")]

exp_sub_m=cbind(normal.exp_m,tumor.exp_m)
exp_sub_m=apply (exp_sub_m, 2,as.integer)
rownames(exp_sub_m)=rownames(normal.exp_m)

getwd()
setwd("D:/mirna/male & female/male")
save(mirna.exp,pheno_updated, exp_sub_m,pheno_sub_m ,file="microRNA-seq.RDATA")


###### DO the differential EXP analysis using DeSeq2
cond1m="Solid Tissue Normal" 
cond2m="Primary Tumor"

dim(exp_sub_m)
dim(pheno_sub_m)

dds_m = DESeqDataSetFromMatrix( countData = exp_sub_m , colData = pheno_sub_m , design = ~ Sample.Type)
dds_run_m = DESeq(dds_m)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res_m=results(dds_run_m)
res_m=results(dds_run_m, contrast = c("Sample.Type",cond1m ,cond2m) )

# remove nulls
res_m=res_m[complete.cases(res_m), ]
summary(res_m)


res_m.df=as.data.frame(res_m)
getwd()
write.table(res_m.df, file = "res_m.txt")

plotMA(res_m, ylim=c(-1,1)) 
summary (res_m)

res_m.degs=res_m.df[res_m.df$padj< 0.05 & abs(res_m.df$log2FoldChange)>log2(2),]

write.table(res_m.degs, file = "res_m_degs.txt")

#for Tfmir
res_m.degs$regulation=0
res_m.degs=as.data.frame(res_m.degs)
res_m.degs[res_m.degs$log2FoldChange<0,]$regulation =-1
res_m.degs[res_m.degs$log2FoldChange>0,]$regulation =1
res_m.degs2= cbind(rownames(res_m.degs), res_m.degs$regulation)
write.table(res_m.degs2, file="res_m_reg.tsv", quote = F, col.names = F, row.names = F)

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




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# main Question

##3) third comparison between male and female :after finishing the pathway we will compare the network of both of them