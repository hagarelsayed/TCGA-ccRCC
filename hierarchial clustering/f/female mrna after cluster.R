
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
source("https://bioconductor.org/biocLite.R")
biocLite("data.table")
BiocManager::install("AnnotationDbi")
install.packages("checkmate")
BiocManager::install("vsn")

library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)
library("dplyr")
clin_1 <- read_delim("D:/mRNA in ccRCC/clinical.tsv", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
class(clin_1)
clin_1_df= as.data.frame(clin_1)

a = clin_1_df[,2]
b= clin_1_df[,13]


clinical = cbind(a,b)
View(clinical)

clinical.df = as.data.frame(clinical)
names(clinical.df)=c("Case ID","Gender")

pheno <- read_delim("D:/mRNA in ccRCC/gdc_sample_sheet.2019-11-20.tsv", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)


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
setwd("D:/mRNA in ccRCC/m & f")
write.table(no_dupl, file = "pheno & gender.txt")

pheno_updated = no_dupl
View(pheno_updated)


rm(cl)
rm(clin_1)
rm(clin_1_df)
rm(clinical)
rm(clinical.df)
rm(p)
rm(no_dupl)
rm(copy_merged)
rm(merged_pheno)
###########################################################


##### load the mrna-Seq data #####

data.path ="D:/mRNA in ccRCC/gdc_download_20191120_094514.709651"
files <- list.files(path=data.path,recursive=T, pattern = "gz")

# read the first file for the first time
file=files[1]

file.id=strsplit(file,"/")[[1]][1]


#open a connection to your gz file and read the file
gz.con=gzfile(file.path(data.path,files[1]))
temp <- read.table(gz.con, header=F)


#create a storing object mrna.exp to save the whole read counts of each file read in an iteration
mrna.exp=temp
rownames(mrna.exp)=mrna.exp[,1]
mrna.exp=mrna.exp[-1]
colnames(mrna.exp)=c(file.id)

for(i in 2: length(files))
{
  
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  gz.con=gzfile(file.path(data.path,files[i]))
  temp <- read.table(gz.con, header=F)
  
  ## remove the first column, bec we had it already
  temp=temp[-1]
  colnames(temp)=c(file.id)
  
  mrna.exp=cbind(mrna.exp,temp)
}
View(mrna.exp)

### do the mapping of ensembel.id to gene symbol###############

# prepare the ensembel id to be like the one in the database
ensemble.id=sapply(rownames(mrna.exp), function(x) strsplit(as.character(x),"\\.")[[1]][1])
View(ensemble.id)
mrna.exp=cbind(ensemble.id,mrna.exp)

mapper<- mapIds(org.Hs.eg.db, keys=ensemble.id, column="SYMBOL",keytype="ENSEMBL", multiVals="first") #??????
mapper.df=as.data.frame(mapper)

mapper.df=cbind(rownames(mapper.df), mapper.df)
names(mapper.df)=c("ensemble.id","symbol")

mrna.exp2=merge(mrna.exp,mapper.df,by="ensemble.id",all.x=T) #????? ensembl.id of y


# drop the first column (ensemble.id)
mrna.exp2=mrna.exp2[-1]

mrna.exp2=mrna.exp2[ ! is.na(mrna.exp2$symbol),]
View(mrna.exp2)
# check duplciation of of gene symbols?  
x=duplicated(mrna.exp2$symbol)  
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
mrna.exp.data=mrna.exp2[-dim(mrna.exp2)[2]]
mrna.exp.data=apply(mrna.exp.data,2, as.numeric)

####remove  duplication by aggregation
mrna.exp.data.agg= aggregate(mrna.exp.data, list(mrna.exp2$symbol),FUN=mean)

rownames(mrna.exp.data.agg)=mrna.exp.data.agg$Group.1
mrna.exp.data.agg=mrna.exp.data.agg[-1]

file.ids=colnames(mrna.exp.data.agg)

old_exp.sub_f= exp.sub_f
exp.sub_f = agg

###### load the mrna sample sheets

###here it's the pheno_updated that we prepared earlier to add the gender colum to it 

# rename column names : replace the spaces with dots
org_pheno_updated = pheno_updated
pheno.names=names(pheno_updated)

#"TCGA-BP-4326-01A" "TCGA-BP-4770-01A"

names(pheno_updated)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))

pheno_updated

table(pheno_updated$Sample.Type,pheno_updated$Gender)


write.table(pheno_updated, file="pheno_updated.tsv", quote = F, col.names = F, row.names = F)
write.table(pheno, file="pheno_only.tsv", quote = F, col.names = F, row.names = F)

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to

file.ids.pheno=pheno_updated$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno_updated$Sample.ID[index.files]



n_f = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("female"))
t_f = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("female"))


all_n_m_samples = pheno_updated[n_m,]$`Sample.ID`
all_t_m_samples = pheno_updated[t_m,]$`Sample.ID`


#############################################################
#############################################################

# 1) first comparison between Female and control : 
#now we will retrieve  the exp data for these female samples only and also the pheno data

## this part is omitted by hajar so that the beloved noha and sahar will no get confused
#  while doing their male RNA analysis :)))

########################################################################
##2) second comparison between male and control : now we will retrieve  the exp data for these
#male samples only and also the pheno data

normal.exp_m = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_n_m_samples]
tumor.exp_m = mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% all_t_m_samples]

pheno.sub_m=pheno_updated[pheno_updated$Sample.ID %in% c(all_n_m_samples,all_t_m_samples), c("Sample.ID", "Sample.Type")]


exp.sub_m=cbind(normal.exp_m,tumor.exp_m)
exp.sub_m=apply (exp.sub_m, 2,as.integer)

rownames(exp.sub_m)=rownames(normal.exp_m)

#DO the differential EXP analysis using DeSeq2
cond1m="Solid Tissue Normal"
cond2m="Primary Tumor"

dim(exp.sub_m)
dim(pheno.sub_m)
dds_m = DESeqDataSetFromMatrix( countData = exp.sub_m , colData = pheno.sub_m , design = ~ Sample.Type)
dds.run_m = DESeq(dds_m)

### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
res_m=results(dds.run_m)
res_m=results(dds.run_m, contrast = c("Sample.Type",cond1m ,cond2m) )


# remove nulls
res_m=res_m[complete.cases(res_m), ]
summary(res_m)

res.df_m=as.data.frame(res_m)


#plotMA(res, ylim=c(-1,1))
#summary (res)

res.degs_m=res.df_m[res.df_m$padj< 0.01 & abs(res.df_m$log2FoldChange)>log2(2),]
library(xlsx)
write.xlsx(res.degs_m, file = "res_degs_m.xlsx")


#for Tfmir
res.degs_m$regulation=0
res.degs_m=as.data.frame(res.degs_m)
res.degs_m[res.degs_m$log2FoldChange<0,]$regulation =-1
res.degs_m[res.degs_m$log2FoldChange>0,]$regulation =1
res.degs2_m= cbind(rownames(res.degs_m), res.degs_m$regulation)
write.table(res.degs2_m, file("res_reg_m.tsv"))


#expression of these degs for plots
exp.degs_m= exp.sub_m[rownames(exp.sub_m) %in% rownames(res.degs_m), ]
write.xlsx(exp.degs_m, file = "exp_degs_m.xlsx")

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



