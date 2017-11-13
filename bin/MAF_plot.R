args = commandArgs (T)
SAMPLE_NAME = args[1]
INPUTFILE = args[2]

library(ggplot2)
library(reshape2)
library(Biobase)


mf <- read.delim(INPUTFILE,header = FALSE, row.names = 1)

mfsub<-mf[mf$V2 > 99,]$V3
pdf(paste(SAMPLE_NAME,"MAF_plot.pdf",sep = ""))
qplot(mfsub,
      geom="histogram",
      binwidth = .01,  
      main = paste("Median MAF is ",median(mfsub),sep=" "), 
      ylab = "Freq",
      xlab = "Minor Allele Frequency",  
      fill=I("blue2"), 
      col=I("red"), 
      alpha=I(1),xlim = c(0,1))
dev.off()