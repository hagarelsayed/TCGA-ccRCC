###########################################
#  -Integrative project                   # 
#  - RNA-Seq                              #
#  - 2019- 11- 18                         #
#  - Copyright: Noha Ismail              #
###########################################
R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

# To install TCGAbiolinks
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# To set working directory
setwd("D:/Integrative project/kidney/TCGAbiolink")

# Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)



query.exp <- GDCquery(project = "TCGA-KIRC", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)
KIRC.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "KIRCExp.rda")

# get clinical data
dataClin <- GDCquery_clinic(project = "TCGA-KIRC","clinical") 

dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 
dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")


dataPrep <- TCGAanalyze_Preprocessing(object = KIRC.exp, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")                
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT") 


# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,dataSmTP],dataFilt[,dataSmNT])

View(dataDEGsFiltLevel)


#library(xlsx)
#write.xlsx(dataDEGs, file = "res_degs_LF1.xlsx")
#or
# DEGs table with expression values in normal and tumor samples


#To get biological process from debugSource

Bio_DEGs <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
                                RegulonList = rownames(dataDEGs))  

#Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)
system.time(Bio_DEGs <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))




TCGAvisualize_EAbarplot(tf = rownames(Bio_DEGs$ResBP), 
                        GOBPTab = Bio_DEGs$ResBP,
                        GOCCTab = Bio_DEGs$ResCC,
                        GOMFTab = Bio_DEGs$ResMF,
                        PathTab = Bio_DEGs$ResPat,
                        nRGTab = rownames(dataDEGs), 
                        nBar = 10)

DAVID_BP_matrix <- TCGAbiolinks:::DAVID_BP_matrix
> RegulonList <- rownames(dataDEGsFiltLevel)
> ResBP <- TCGAbiolinks:::TCGAanalyze_EA(GeneName="DEA genes Normal Vs Tumor",
                                         +                                        RegulonList = RegulonList,
                                         +                                        TableEnrichment = DAVID_BP_matrix,
                                         +                                        EAGenes = EAGenes,
                                         +                                        GOtype = "DavidBP",
                                         +                                        FDRThresh = 0.01,
                                         +                                        GeneSymbolsTable = TRUE)