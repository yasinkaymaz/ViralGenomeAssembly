.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/3.2/"))

library(ggplot2)
setwd("/home/yasinkaymaz/Documents/EBV/Batch4/Coverage/Allreads/")
#Function to put the data in the correct format for plotting
preprocessCoverageFile <- function(file, timeSeries) {
  # Load coverage data
  data = read.delim(file, header=T)
  # Remove everything before ":"
  data$Locus<-gsub(".*:","",data$Locus)
  data$Locus<-as.numeric(data$Locus)
  data<-subset(data, select=c(Locus, Total_Depth))
  data$Amplification<-c(as.character(timeSeries))
  data
}

#Function to plot the coverage
plotCoverage <- function(coverage, title) {
  ggplot(data=coverage, aes(x=Locus, y=Total_Depth, colour=Amplification)) +
    geom_line() +
    ylab("Depth of coverage") +
    xlab("Loci") +
    labs(title=title)
}


# Example of using the above functions.

outputPDFFile = "Coverage_of_depth_EBV_PCRsWGA_AllReads.pdf"
# Load coverage data
coverageControl = preprocessCoverageFile("Raji.gatk_recal_DepthofCoverage", "Control Raji")
coverageSample1 = preprocessCoverageFile("Raji_5cyc.gatk_recal_DepthofCoverage", "Raji - 5 cyc PCR -sWGA")
coverageSample2 = preprocessCoverageFile("Raji_10cyc.gatk_recal_DepthofCoverage", "Raji - 10 cyc PCR -sWGA")
coverageSample3 = preprocessCoverageFile("Raji_15cyc.gatk_recal_DepthofCoverage", "Raji - 15 cyc PCR -sWGA")
coverageSample4 = preprocessCoverageFile("Raji_Fw_only.gatk_recal_DepthofCoverage", "Raji - Fw only PCR -sWGA")
coverageSample5 = preprocessCoverageFile("Raji_Fw_only.gatk_recal_DepthofCoverage", "Raji - Rev only PCR -sWGA")

head(coverageControl)
mean(coverageControl$Total_Depth)
mean(coverageSample3$Total_Depth)

# Create dataset with all data together.
allCoverage1<-rbind(coverageControl,coverageSample1)
allCoverage2<-rbind(coverageControl,coverageSample2)
allCoverage3<-rbind(coverageControl,coverageSample3)
allCoverage4<-rbind(coverageControl,coverageSample4)
allCoverage5<-rbind(coverageControl,coverageSample5)

# Plot the data to a pdf
#pdf(outputPDFFile, paper="a4r", width=10)
pdf(outputPDFFile, width=30)
#plotCoverage(coverageControl, "Coverage plot for control Raji")
#plotCoverage(coverageSample1, "Coverage plot for 5cyc-PCR-sWGA")
plotCoverage(allCoverage1, "Coverage plot for all 5cyc-PCR-sWGA")
plotCoverage(allCoverage2, "Coverage plot for all 10cyc-PCR-sWGA")
plotCoverage(allCoverage3, "Coverage plot for all 15cyc-PCR-sWGA")
plotCoverage(allCoverage4, "Coverage plot for all Fw-only-PCR-sWGA")
plotCoverage(allCoverage5, "Coverage plot for all Rev-only-PCR-sWGA")

dev.off()

data = read.table("5libraries.DepthofCoverage",header = TRUE)
head(data)
cor(data$Depth_for_Raji, data$Depth_for_Raji_5cyc_sWGA)
cor(data$Depth_for_Raji, data$Depth_for_Raji_10cyc_sWGA)
cor(data$Depth_for_Raji, data$Depth_for_Raji_15cyc_sWGA)
cor(data$Depth_for_Raji, data$Depth_for_Raji_Fw_sWGA)
cor(data$Depth_for_Raji, data$Depth_for_Raji_Rev_sWGA)

genes <- read.table("depthwithgenes.sample_interval_summary", header = TRUE,row.names = 1)

meancov <- data.frame(
      Samples=c("Raji","Raji_5cyc","Raji_10cyc","Raji_15cyc","Raji_Fw","Raji_Rev"),
      rbind(genes$Raji_mean_cvg,
      genes$Raji_5cyc_sWGA_mean_cvg,
      genes$Raji_15cyc_sWGA_mean_cvg,
      genes$Raji_10cyc_sWGA_mean_cvg,
      genes$Raji_Fw_sWGA_mean_cvg,
      genes$Raji_Rev_sWGA_mean_cvg))


barplot(meancov,beside = TRUE)

ggplot(data=meancov, aes(colour=Samples)) +
  geom_bar() +
  ylab("Depth of coverage") +
  xlab("Loci") +
  labs(title=title)


