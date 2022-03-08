
library(GOplot)

library(readr)

genelist = as.data.frame(res_degs_f1)
class(genelist)
getwd()
setwd("D:/ccRCC project/Goplot")


library(xlsx)

library(readxl)
res_degs_m <- read_excel("D:/ccRCC project/david/res_degs_m.xlsx")


Male_David <- read_delim("D:/ccRCC project/david/Male David.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

david = as.data.frame(Male_David)
library("stringr")
#?str_split_fixed
#bp_m_split=str_split_fixed(bp_m_a$GOTERM_BP_DIRECT , "~", 2)
splitted = str_split_fixed(david$Term, "~", 2)
spl_df=as.data.frame(splitted)
colnames(spl_df) =c("ID","Term")
davidd = cbind.data.frame(spl_df,david)

a = davidd
imp_david = a[, c(1:3,7,8)];
rm(david)
rm(spl_df)
rm(res_degs_f1)
rm(splitted)

library(xlsx)

library(readxl)
res_degs_m <- read_excel("D:/ccRCC project/david/res_degs_m.xlsx")

