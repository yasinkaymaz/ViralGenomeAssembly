library(cowplot)# cowplot enables side-by-side ggplots
library(reshape2)
library(dplyr)
library(ggplot2)
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
