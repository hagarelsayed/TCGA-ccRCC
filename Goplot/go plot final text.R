
library(readr)
rm(female_david_pb)

all <- read_delim("D:/ccRCC project/david/hagar/f/all.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
res_degs_f1 <- read_table2("D:/ccRCC project/mRNA in ccRCC/m & f/female mRNA/res_degs_f1.txt")

female_david_pb <- read_delim("D:/ccRCC project/david/hagar/f/female david pb.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

library(GOplot)
library(readr)

genelist = as.data.frame(res_degs_f1)
class(genelist)
getwd()
setwd("D:/ccRCC project/Goplot")

# 
# library(xlsx)
# 
# library(readxl)
# res_degs_m <- read_excel("D:/ccRCC project/david/res_degs_m.xlsx")
# 
# 
# Male_David <- read_delim("D:/ccRCC project/david/Male David.txt", 
#                          "\t", escape_double = FALSE, trim_ws = TRUE)

david = as.data.frame(female_david_pb)
library("stringr")
#?str_split_fixed
#bp_m_split=str_split_fixed(bp_m_a$GOTERM_BP_DIRECT , "~", 2)
splitted = str_split_fixed(david$Term, "~", 2)
spl_df=as.data.frame(splitted)
colnames(spl_df) =c("ID","Term")
davidd = cbind.data.frame(spl_df,david)

# 
# imp_david = a[, c(1:3,5,7,8,12:15)];
# rm(david)
# rm(a)
# rm(spl_df)
# rm(res_degs_f1)
# rm(splitted)
# rm(female_bp_cc_ts)
# rm(female_david_pb)

genelist = res_degs_f1
a=as.data.frame(sapply(genelist, function(x) gsub("\"", "", x)))

a.names=names(a)
names(a)= as.character( sapply ( a.names, function(x) gsub("\"", "",x)))

names(a)[1]="ID"
d= a
# Gene names and their logFC 
names(d)[1] = "ID"
names(d)[2] = "logFC"


save (genelist , file = "res.degs.f22.Rdata")
write.table(genelist, file("res.degs.f22.tsv"))
rm(b)
genelist = a

EC = list (david = imp_david,genelist = genelist)

david_df1 = imp_david [1:5]
# allTraits = traitData[, -c(31, 16)];
# allTraits = allTraits[, c(2, 11:36) ];
genelist_df1 = d[, c(1,2)]
genelist = data.frame(genelist_df1$ID,genelist_df1$lfcSE)
names(genelist) = c("ID","logFC")

a = davidd
bp= davidd
bp= data.frame(bp$Category,bp$ID,bp$Term,bp$Genes,bp$Benjamini)
names(bp) = c("Category","ID","Term","Genes","adj_pval")

bp$Category= "BP"
genelist = data.frame(genelist_df1$ID,genelist_df1$logFC)
names(genelist) = c("ID","logFC")
BP = bp
circ <- circle_dat(BP,genelist)
GOBar(subset(circ, category == 'BP'))

GOBar(circ, display = 'multiple')

GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
getwd()
GOBubble(circ, labels = 3)
options(stringsAsFactors = FALSE);

