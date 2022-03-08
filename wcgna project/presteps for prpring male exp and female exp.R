
###########################################
#  -Integrative bioinformatics            # 
#  - miRNA-Seq  for male                  #
#  - 2019- 12- 18                         #
#  - Copyright: Hajar Elsayed             #
###########################################
# after installing global env to the r 
# do these presteps to get the female expression for both normal and tumor
# also for male normal and tumor

file.ids.pheno=pheno_updated$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno_updated$Sample.ID[index.files]

n_f = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("female"))
t_f = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("female"))
n_m = which(pheno_updated$`Sample.Type` %in% c("Solid Tissue Normal")  &  pheno_updated$Gender %in% c("male"))
t_m = which(pheno_updated$`Sample.Type` %in% c("Primary Tumor")  &  pheno_updated$Gender %in% c("male"))

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
exp.sub_m=apply (exp.sub_m, 2,as.integer)



library(readr)

#there is something wrong in the tables
# write.table(exp.sub_f, file="mRNAexp_f_n&t.tsv", quote = F, col.names = T, row.names = T)
# write.table(pheno, file="pheno_only.tsv", quote = F, col.names = T, row.names = T)
# write.table(pheno_updated, file="pheno_updated.tsv", quote = F, col.names = T, row.names = T)
 

#exp.sub_f=apply (exp.sub_m, 2,as.integer)
#rownames(exp.sub_m)=rownames(normal.exp_m)
#rownames(exp.sub_m)=rownames(normal.exp_m)

rm(mapper.df)
rm(mrna.exp)
rm(mrna.exp.data)
rm(mrna.exp2)
rm(temp)


rm(file)
rm(b)
rm(a)
rm(data.path)
rm(dupl)
rm(dupl.df)
rm(ensemble.id)
rm(file.id)
rm(gz.con)
rm(i)
rm(mapper)
rm(x)

y<- read_table2("D:/ccRCC project/mRNA in ccRCC/m & f/female mRNA/res_degs_f1.txt")
res.degs = y
  
res.degs$"regulation"=0
res.degs=as.data.frame(res.degs)
res.degs[res.degs$log2FoldChange<0,]$"regulation" =-1
res.degs[res.degs$log2FoldChange>0,]$"regulation" =1
res.degs2_m= cbind(rownames(res.degs_m), res.degs_m$regulation)
write.table(res.degs2_m, file("res_reg_m.tsv"))
