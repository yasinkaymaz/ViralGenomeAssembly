#.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/3.2/"))
args = commandArgs (T)
SAMPLE_NAME = args[1]

ref1convert <- read.delim("/home/yk42w/codes/EBVseq/Type1_to_type2_reflocmatch.txt", header=TRUE)

ref2convert <- read.delim("/home/yk42w/codes/EBVseq/Type2_to_type1_reflocmatch.txt", header=TRUE)

library(ggplot2)
#setwd("/home/yasinkaymaz/Documents/EBV/Batch4/Coverage/Allreads/")
#Function to put the data in the correct format for plotting

preprocessCoverageFile <- function(file, timeSeries, gtype) {
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
  data.m
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

outputPDFFile = "Coverage_of_depth_EBV.pdf"
# Load coverage data
coverageControl = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_007605.sorted_DepthofCoverage",sep=""), "Type 1","NC_007605")
coverageSample1 = preprocessCoverageFile(paste(SAMPLE_NAME,"_alignment_NC_009334.sorted_DepthofCoverage",sep=""), "Type 2","NC_009334")

head(coverageControl)
mean(coverageControl$Total_Depth)
mean(coverageSample1$Total_Depth)

# Create dataset with all data together.
allCoverage1<-rbind(coverageControl,coverageSample1)

# Plot the data to a pdf
pdf(outputPDFFile, width=30)
plotCoverage(allCoverage1, paste("Coverage comparison for Genome Typing",SAMPLE_NAME,sep = " "))

dev.off()
