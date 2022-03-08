library(readr)
DAVID_Kegg_pathway_male <- read_delim("D:/ccRCC project/david/DAVID-Kegg pathway - male.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
install.packages('GOplot')
library("GOplot")
library(readr)
DAVID_GO_BP_for_male <- read_delim("D:/ccRCC project/david/DAVID- GO- BP for male.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
DAVID_OMIM_disease_male <- read_delim("D:/ccRCC project/david/DAVID-OMIM disease-male.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
DAVID_GO_MF_male <- read_delim("D:/ccRCC project/david/DAVID-GO-MF-male.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
DAVID_GO_CC_male <- read_delim("D:/ccRCC project/david/DAVID-GO-CC-male.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
bp_m_a = DAVID_GO_BP_for_male
mf_m = DAVID_GO_MF_male
cc_m = DAVID_GO_CC_male
kegg_m = DAVID_Kegg_pathway_male
omim_m = DAVID_OMIM_disease_male


library("stringr")
?str_split_fixed
bp_m_split=str_split_fixed(bp_m_a$GOTERM_BP_DIRECT , "~", 2)
colnames(bp_m_split) =c("Category","Term")
bp_m = cbind.data.frame(bp_m_a,bp_m_split)

library(readr)
female_TAM_enrichment <- read_delim("C:/Users/EZZ/Downloads/TAM/female TAM enrichment.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
