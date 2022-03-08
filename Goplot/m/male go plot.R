
library(readr)
male_bp_david <- read_delim("D:/ccRCC project/Goplot/m/male bp david.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)



# res_degs_f1 <- read_table2("D:/ccRCC project/mRNA in ccRCC/m & f/female mRNA/res_degs_f1.txt")

# female_david_pb <- read_delim("D:/ccRCC project/david/hagar/f/female david pb.txt", 
#                               "\t", escape_double = FALSE, trim_ws = TRUE)

library(GOplot)
library(readr)
res_degs <- read_delim("D:/ccRCC project/Goplot/m/res.degs.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)


library(GOplot)
library(readr)
genelist = as.data.frame(res_degs)

getwd()
setwd("D:/ccRCC project/Goplot")
                          "\t", escape_double = FALSE, trim_ws = TRUE)

david = as.data.frame(male_bp_david)
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

genelist = res_degs
a=as.data.frame(sapply(genelist, function(x) gsub("\"", "", x)))

a.names=names(a)
names(a)= as.character( sapply ( a.names, function(x) gsub("\"", "",x)))

# Gene names and their logFC 
names(d)[1] = "ID"
names(d)[2] = "logFC"


a= genelist


# genelist_df1 = a[, c(1,2)]
genelist = data.frame(a$X1,a$log2FoldChange)
names(genelist) = c("ID","logFC")


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
?GOBubble
GOBubble(circ, table.col = T)

GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  

GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  

reduced_circ <- reduce_overlap(circ, overlap = 0.75)

GOBubble(reduced_circ, labels = 2.8)

GOCircle(circ)

IDs <- c('GO:0007588', 'GO:0006814', 'GO:0006811', 'GO:0035725', 'GO:0030049', 'GO:0055085', 'GO:0055078', 'GO:0051180','GO:0034220',"GO:1902476")
GOCircle(circ, nsub = IDs)

#Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)

#b = character(head(genelist$ID))
process = bp$Term[1:7]

chord <- chord_dat(circ, genelist , process)
?chord

head(chord)
chord <- chord_dat(data = circ, genes = genelist)
chord <- chord_dat(data = circ, process = process)

chord <- chord_dat(data = circ, genes = genelist, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

set.seed(24)
GOChord(chord, gene.order = 'logFC')

GOHeat(chord[,-8], nlfc = 0)

GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

GOCluster(circ, process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
