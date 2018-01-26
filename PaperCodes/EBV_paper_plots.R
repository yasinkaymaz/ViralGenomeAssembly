library(cowplot)# cowplot enables side-by-side ggplots
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Figure 1-B

data <- read.delim("~/Documents/EBV/fixed_assembly_genome_CoveredLen.txt",header = TRUE,row.names = 3)
data$samplename <- row.names(data)

p <- ggplot(data, aes(x=reorder(samplename,lineorder), y=CoverageLen/1000, fill=SampleType, group=SampleType ))+
  geom_bar(stat="identity")+
  theme_bw()+
  ylab("Covered Genome Length (kb)")+
  xlab("Sequenced EBV Genomes")+
  theme(axis.text.x=element_blank())+
  geom_point(aes(x=reorder(samplename,lineorder), y=N50/1000))+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))+
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/genome_CoveredLen.pdf",width = 15,height = 6)
p
dev.off()

#Figure 2-A

data <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/ControlGenomesDistance.txt",row.names = 1,header = TRUE)
data$Genomes <- gsub("_genome_fixed_assembly_Fixed","",rownames(data))
datasub <- data[,c("T1_MMR", "T2_MMR","Jijoye_Pct")]
mdatasub <- melt(datasub,id="Jijoye_Pct")

err.j <- ggplot(data=mdatasub,aes(x=100*(value), y=Jijoye_Pct, color=variable))+geom_point(size=2)+
  xlab("Distance (% missmatch)") + ylab("Jijoye % in mixture")+
  scale_color_manual(labels = c("Type 1 (NC_007605)", "Type 2 (NC_009334)"), values = c("blue", "red"))+
  theme_bw()+
  guides(color=guide_legend("Distance to"))
#+geom_smooth(data = mdatasub, method = "auto",aes(fill=variable), formula = y ~ log10(x))

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/ControlGenomesDistance.pdf",width = 7.5,height = 5)
err.j
dev.off()

df <- data.frame(Daudi=seq(100,0),
                 Jijoye=seq(0,100))
mdf <- melt(df)
mdf$sampleid <- rep(1:101,2)


pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/ControledMixtures.pdf",width = 7.5,height = 5)
ggplot(mdf, aes(sampleid, value, color=variable,fill=variable))+geom_bar(stat="identity")+coord_flip()+
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.background =element_blank())
dev.off()

#Figure 2-B
viral <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/ViralLoadvsType.txt",header = TRUE, row.names = 1)
head(viral)

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/ViralLoadvsType.pdf")
ggplot(data=viral)+
  geom_jitter(aes(colour=genome,genome,viralload),size=3)+
  geom_boxplot(aes(colour=viral$genome,viral$genome,viral$viralload),outlier.colour = NA, alpha = 0.4)+
  theme(text = element_text(size=20))+
  labs(colour="Viral Type",x='',y="Viral Load (copy/ng DNA)")
dev.off()



#Figure 2-C
data <- read.delim("~/Documents/EBV/type.genomes.txt",header = TRUE,row.names = 1)

p <- ggplot(data, aes(x=as.factor(SampleType),group=as.factor(EBVtype),fill=as.factor(EBVtype)))+
  geom_bar(position = "dodge")+
  xlab("Case/Control")+
  ylab("Number of Individuals")+
  scale_fill_manual("EBV Type",values = c("1"="blue","2"="red"))

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EBVtypeTest.pdf")
p
dev.off()


data <- read.delim("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/dnds_2.txt", row.names = 1, header = TRUE)
head(data)
subdata <- data[,c("gene","SynPerGenomePerKb","NonsynPerGenomePerKb")]
#subdata$genes <- rownames(data)
subdata <- melt(subdata)
subdata$lineorder <- rep(seq(1,87),2)

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/Syn-NonSyn_Profiles.pdf",width = 18,height = 5)
ggplot(subdata, aes(x=reorder(gene,lineorder), value,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  ylab("Number of variant per Genome per Kb")+
  xlab("Viral Genes (in genomic order)")
dev.off() 



i=0
plots <- list() 
for (genome in c("Jijoye","Namalwa","eBL-Tumor-0033", "eBL-Plasma-0049","eBL-Tumor-0012","LN827563.2_sLCL-1.18")){
  print(genome)
  i=i+1
  data1 <- read.delim(paste("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/",genome,"_Mismatch_positions_with_type1_smooth.bed",sep = ""),header = TRUE)
  data2 <- read.delim(paste("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/",genome,"_Mismatch_positions_with_type2_smooth.bed",sep = ""),header = TRUE)
  data <- data.frame(Pos=data1$start,Sim2Type1=data1$Sim2Type1,Sim2Type2=data2$Sim2Type2)
  mdata <- melt(data, id=c("Pos"))
  
  p1 <- ggplot(mdata, aes(Pos/1000, y=value, colour=variable))+geom_line(size=1)+
    scale_colour_manual(labels = c("Similarity to Type 1", "Similarity to Type 2"), values=c("blue", "red"))+
    xlab("Genomic Postion (kb)")+
    ylab("Percent Similarity")+
    scale_x_continuous(breaks = round(seq(min(mdata$Pos/1000), max(mdata$Pos/1000), by = 20) ),expand = c(0, 0))+
    theme(axis.text.x = element_text(angle = 0, hjust = .5,vjust=0.5))

plots[[i]] <- p1
}

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/HybridGenomes-3.pdf",width = 18,height = 16)
multiplot(plotlist = plots)
dev.off()

#Hybrid genome similarity for control. Supp Figure

i=0
plots <- list() 
for (genome in c("Jijoye","Daudi","Namalwa","Raji","eBL-Tumor-0033", "eBL-Plasma-0049","eBL-Tumor-0012","LN827563.2_sLCL-1.18")){
  print(genome)
  i=i+1
  data1 <- read.delim(paste("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/",genome,"_Mismatch_positions_with_type1_smooth.bed",sep = ""),header = TRUE)
  data2 <- read.delim(paste("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/",genome,"_Mismatch_positions_with_type2_smooth.bed",sep = ""),header = TRUE)
  data <- data.frame(Pos=data1$start,Sim2Type1=data1$Sim2Type1,Sim2Type2=data2$Sim2Type2)
  mdata <- melt(data, id=c("Pos"))
  
  p1 <- ggplot(mdata, aes(Pos/1000, y=value, colour=variable))+geom_line(size=1)+
    scale_colour_manual(labels = c("Similarity to Type 1", "Similarity to Type 2"), values=c("blue", "red"))+
    xlab("Genomic Postion (kb)")+
    ylab("Percent Similarity")+
    scale_x_continuous(breaks = round(seq(min(mdata$Pos/1000), max(mdata$Pos/1000), by = 20) ),expand = c(0, 0))+
    theme(axis.text.x = element_text(angle = 0, hjust = .5,vjust=0.5))
  
  plots[[i]] <- p1
}

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/HybridGenomes-Supplementary.pdf",width = 18,height = 16)
multiplot(plotlist = plots)
dev.off()


#Deleted Genome CoveragePlot
i=0
plots <- list()
for (genome in c("eBL-Tumor-0031", "Raji_CellLine_longRead","Daudi_CellLine")){
  i=i+1
  print(genome)
  cov <- read.delim(paste("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/",genome,"_alignment_NC_007605.sorted_DepthofCoverage",sep=""),row.names = 1, header=TRUE)
  rownames(cov) <- gsub("NC_007605:","",rownames(cov))
  d <- data.frame(Pos=NULL,Cov=NULL)
  wind=1000
  slid=500

  for (x in seq(0,max(as.numeric(rownames(cov)))-wind,slid ) ){
    d <- rbind(d, data.frame(Pos=x, Cov=sum(cov[x:x+wind, ]$Total_Depth) ))
  }

  p1 <- ggplot(d, aes(x=Pos/wind,y=Cov))+
    geom_area(size=1,fill="turquoise2",color="brown")+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(breaks = round(seq(min(d$Pos/wind), max(d$Pos/wind), by = 20) ),expand = c(0, 0))+
    theme(axis.text.x = element_text(angle = 0, hjust = .5,vjust=0.5))+
    xlab("Genomic Postion (kb)")+
    ylab("Depth of Coverage")

  plots[[i]] <- p1

}
#geom_line(size=1,color="brown")

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/DeletedGenomes2.pdf",width = 18,height = 12)
multiplot(plotlist = plots)
dev.off()

#PCA with MSA
#library(bios2mds)
#aln <- import.fasta(system.file("/Users/yasinkaymaz/Documents/EBV/ourbatch/CorrectFastas/my_ICed_A73.aln.filtered.fasta", package = "bios2mds"))




df <- data.frame(date=c("12Ar","13Ar","13Ar"),
                 Mtype=c("m1","m2","m3"),
                 Count=c(3,1,1))
df
library(reshape)

udf <- untable(df[,c(1,2)], num=df[,3])

tdf <- as.data.frame(table(udf[,c(1,2)]))
tdf


for (day in seq(1,28,7) ){
  print(day)
}


#Figure 5C

library(devtools)
#install_github("drveera/ggman")
library(ggman)

#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.bed",header = TRUE, row.names = 4)
#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.1Milperm.bed",header = TRUE, row.names = 4)
#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.1Milperm-strata.bed",header = TRUE, row.names = 4)
assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.combined1.bed",header = TRUE, row.names = 4)

assoc <- assoc[assoc$P < 0.75,]
assoc.ebna3c <- assoc[ which(89135 > assoc$position & assoc$position > 86083),]$position
assoc[ -log10(assoc$P) > 3,]
write.table(assoc[ -log10(assoc$P) > 2,], file="~/Dropbox/Papers/EBV_project/workspace/data/SignificantVariatAssociations.txt",sep="\t")

p1 <- ggman(assoc, snp = "position", bp = "position", chrom = "chrom", pvalue = "P",relative.positions = TRUE,pointSize = 1.5)

p1 <- ggmanHighlight(p1,highlight = assoc.ebna3c)+ theme_bw()

ggsave(plot=p1,width = 10,height = 6, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/assoc3.pdf", useDingbats=FALSE )


ggman(assoc, snp = "snp", bp = "position", chrom = "chrom", pvalue = "OR",relative.positions = TRUE,pointSize = 0.5, logTransform = FALSE, ymax = 100)
