library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

filename <- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)
c6 <- read.gmt(filename)
yourEntrezIdList<- c(55365,8839,1311,5919,28959,28959,894,22801,10974,
                     2308,55799,1528,4837,7439,22809,1240,85007,10418) #ENTREZID of DE genes and Annovar genes
ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.05)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes.csv")
goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont= "ALL",pAdjustMethod="BH",pvalueCutoff = 0.05,readable= TRUE)
write.csv(goEnrich,"MyGORelatedGenes.csv")
keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH", pvalueCutoff = 0.05)
write.csv(keggEnrich,"MyKEGGRelatedGenes.csv")
quit(save="no")
 