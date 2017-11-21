.libPaths(c(.libPaths(),"~/n/home13/yasinkaymaz/biotools/Rlibs/"))

library(plyr)
#library(dplyr)

args = commandArgs (T)
InputTranscriptsFile = args[1]

expdata <- read.delim(InputTranscriptsFile,row.names = 1,header = TRUE)
#expdata <- read.delim("/Users/yasinkaymaz/Documents/Harvard_Informatics/Data_Explore/test.salmon.merged.TPM.txt",row.names = 1,header = TRUE)
head(expdata)
library(stringr)
expdata$gene_id <- str_sub(expdata$X.1,18,-5)
expdata <- expdata[,-which(names(expdata) %in% c("X.1"))]

expdata.genes <- ddply(expdata,"gene_id", numcolwise(sum))
head(expdata.genes)
rownames(expdata.genes) <- expdata.genes$gene_id
expdata.genes <- expdata.genes[,-which(names(expdata.genes) %in% c("gene_id"))]

colSums(expdata.genes)
write.table(expdata.genes, file = paste(InputTranscriptsFile,"_genes",sep = ""),sep = "\t",quote = FALSE,col.names=NA)
