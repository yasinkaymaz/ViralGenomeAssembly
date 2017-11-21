#.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/3.2/"))
args = commandArgs (T)
SAMPLE_NAME = args[1]

ref1convert <- read.delim("/home/yk42w/codes/EBVseq/Type1_to_type2_reflocmatch.txt", header=TRUE)

ref2convert <- read.delim("/home/yk42w/codes/EBVseq/Type2_to_type1_reflocmatch.txt", header=TRUE)

library(ggplot2)
#setwd("/home/yasinkaymaz/Documents/EBV/Batch4/Coverage/Allreads/")
#Function to put the data in the correct format for plotting

preprocessCoverageFile <- function(file, timeSeries, gtype, ebna) {
  # Load coverage data
  data = read.delim(file, header=T)
  
  if (gtype == "NC_007605"){
    data.m <- merge(ref1convert, data, by.x='RefLocus',by.y = 'Locus')
  }
  else{
    data.m <- merge(ref2convert, data, by.x='RefLocus',by.y = 'Locus')
  }
  # Remove everything before ":"
  data.m$AlnLocus<-gsub(".*:","",data.m$AlnLocus)
  data.m$AlnLocus<-as.numeric(data.m$AlnLocus)
  data.m<-subset(data.m, select=c(AlnLocus, Total_Depth))
  data.m$Genomes<-c(as.character(timeSeries))
  attach(data.m)
  if (ebna == "EBNA2"){
    data.m <- data.m[ which(AlnLocus > 36000 & AlnLocus < 38000),]
  }
  else{
    data.m <- data.m[ which(AlnLocus > 80000 & AlnLocus < 91000),]
  }
}

#Function to plot the coverage
plotCoverage <- function(coverage, title) {
  ggplot(data=coverage, aes(x=AlnLocus, y=Total_Depth, colour=Genomes)) +
    geom_line() +
    ylab("Depth of coverage") +
    xlab("Loci") +
    labs(title=title)
}


# Example of using the above functions.

outputPDFFile = "Coverage_of_depth_EBV_onlyEBNAs.pdf"
# Load coverage data
coverage.t1.ebna2 = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_007605.sorted_DepthofCoverage",sep=""), "Type 1","NC_007605","EBNA2")
coverage.t2.ebna2 = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_009334.sorted_DepthofCoverage",sep=""), "Type 2","NC_009334","EBNA2")

coverage.t1.ebna3 = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_007605.sorted_DepthofCoverage",sep=""), "Type 1","NC_007605","EBNA3s")
coverage.t2.ebna3 = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_009334.sorted_DepthofCoverage",sep=""), "Type 2","NC_009334","EBNA3s")


# Create dataset with all data together.
bothCoverage.ebna2 <- rbind(coverage.t1.ebna2,coverage.t2.ebna2)
bothCoverage.ebna3 <- rbind(coverage.t1.ebna3,coverage.t2.ebna3)

# Plot the data to a pdf
pdf(outputPDFFile, width=30)

if (sum(coverage.t1.ebna2$Total_Depth) > sum(coverage.t2.ebna2$Total_Depth)){
  p<-t.test(coverage.t1.ebna2$Total_Depth,coverage.t2.ebna2$Total_Depth)$p.val
  print(paste("Type 1 --- EBNA2, pval=",p,   sep = " "))
}else{
  p<-t.test(coverage.t1.ebna2$Total_Depth,coverage.t2.ebna2$Total_Depth)$p.val
  print(paste("Type 2 --- EBNA2, pval=",p,   sep = " "))
}
plotCoverage(bothCoverage.ebna2, paste("Coverage comparison of EBNA2 for Genome Typing",SAMPLE_NAME,sep = " "))

if (sum(coverage.t1.ebna3$Total_Depth) > sum(coverage.t2.ebna3$Total_Depth)){
  p<-t.test(coverage.t1.ebna3$Total_Depth,coverage.t2.ebna3$Total_Depth)$p.val
  print(paste("Type 1 --- EBNA3, pval=",p,   sep = " "))
}else{
  p<-t.test(coverage.t1.ebna3$Total_Depth,coverage.t2.ebna3$Total_Depth)$p.val
  print(paste("Type 2 --- EBNA3, pval=",p,   sep = " "))
}
plotCoverage(bothCoverage.ebna3, paste("Coverage comparison of EBNA3A/B/C for Genome Typing",SAMPLE_NAME,sep = " "))

dev.off()
