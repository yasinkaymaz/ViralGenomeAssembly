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

err.j <- ggplot(data=mdatasub,aes(x=100*(value), y=Jijoye_Pct, color=variable))+geom_point(size=3)+
  xlab("Distance (% missmatch)") + ylab("Jijoye % in mixture")+
  scale_color_manual(labels = c("Type 1 (NC_007605)", "Type 2 (NC_009334)"), values = c("blue", "red"))+
  theme_bw()+
  guides(color=guide_legend("Distance to"))
#+geom_smooth(data = mdatasub, method = "auto",aes(fill=variable), formula = y ~ log10(x))

ggsave(plot=err.j,width = 11,height = 9, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/ControlGenomesDistance2.pdf", useDingbats=FALSE )

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




vlplot <- ggplot(data=viral,aes(colour=genome,genome,viralload))+  
          geom_violin(position=position_dodge(1))+
          geom_dotplot(aes(fill=genome),binaxis='y', stackdir='center', dotsize=1,position=position_dodge(1))+
          theme(text = element_text(size=20))+
          labs(colour="Viral Type",x='',y="Viral Load (copy/ng DNA)")

ggsave(plot=vlplot,width = 12,height = 10, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/ViralLoadvsType2.pdf", useDingbats=FALSE )



#Figure 2-C
data <- read.delim("~/Documents/EBV/type.genomes.txt",header = TRUE,row.names = 1)
names(data)
datasub <- data[,c("EBVtype", "SampleType")]
datasub <- rbind(datasub, data.frame(EBVtype =c("1","2","1"), SampleType=c("BL","BL","BL")))

p <- ggplot(datasub, aes(x=as.factor(SampleType),group=as.factor(EBVtype),fill=as.factor(EBVtype)))+
  geom_bar(position = "dodge")+
  xlab("Case/Control")+
  ylab("Number of Individuals")+
  scale_fill_manual("EBV Type",values = c("1"="blue","2"="red"))

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EBVtypeTest_v2.pdf")
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
  
  p1 <- ggplot(mdata, aes(Pos/1000, y=value, colour=variable))+geom_line(size=1,show.legend = FALSE)+
  #p1 <- ggplot(mdata, aes(Pos/1000, y=value, colour=variable))+geom_area(size=1)+
    scale_colour_manual(labels = c("Similarity to Type 1", "Similarity to Type 2"), values=c("blue", "red"))+
    xlab("Genomic Postion (kb)")+
    ylab("Percent Similarity")+
    scale_x_continuous(breaks = round(seq(min(mdata$Pos/1000), max(mdata$Pos/1000), by = 20) ),expand = c(0, 0))+
    theme(axis.text.x = element_text(angle = 0, hjust = .5,vjust=0.5))
  print(p1)
plots[[i]] <- p1
}

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/HybridGenomes-6.pdf",width = 18,height = 10)
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
plots.d <- list()
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

  plots.d[[i]] <- p1

}
#geom_line(size=1,color="brown")

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/DeletedGenomes2-2.pdf",width = 18,height = 6)
multiplot(plotlist = plots.d)
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

#library(devtools)
#install_github("drveera/ggman")
#library(ggman)

#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.bed",header = TRUE, row.names = 4)
#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.1Milperm.bed",header = TRUE, row.names = 4)
#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.filtered.1Milperm-strata.bed",header = TRUE, row.names = 4)
#assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/assocP.combined1.bed",header = TRUE, row.names = 4)
#assoc <-read.delim("/Users/yasinkaymaz/Documents/EBV/ourbatch2/Secondround/test/Plink.1M.nullFixed.sub",header = TRUE)
assoc <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/AssocTest.03-27-2018_v3.txt", header = TRUE)
head(assoc)
dim(assoc)
#assoc <- assoc[assoc$P < 0.75,]
#assoc <- assoc[assoc$P < 1.0,]
#assoc <- droplevels(assoc[assoc$FILTER == "PASS",])

assoc.sig <- assoc[ which(assoc$P < 0.01),]$position

assoc.save <- droplevels(assoc[ -log10(assoc$P+0.0000001) > 2,c("chrom", "position", "REF", "ALT", "MINA","MINU", "OBSA","OBSU","P","OR","VarAAposition",  "VarType", "RefCodon", "VarCodon", "RefAA", "VarAA",  "Gene")])
write.table(assoc.save, file="~/Dropbox/Papers/EBV_project/workspace/data/SignificantVariatAssociations.txt",sep="\t",quote = FALSE,row.names = FALSE)

# p1 <- ggman(assoc,sigLine = 2, snp = "position", bp = "position", chrom = "chrom", pvalue = "P",relative.positions = TRUE,pointSize = 1.5)
# 
# #p1 <- ggmanHighlight(p1,highlight = assoc.ebna3c)+ theme_bw()
# p1 <- ggmanHighlight(p1,highlight = assoc.sig)+ theme_bw()
# 
# ggsave(plot=p1,width = 12,height = 6, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/assoc5.pdf", useDingbats=FALSE )

head(assoc)
p1 <- ggplot(assoc, aes(x=position/1000, y=-log10(P)))+ geom_point(size=4, color="grey")+
  geom_hline(aes(yintercept= as.numeric(2)),colour = "red", linetype="dashed", size = 1) + theme_bw()+
  geom_point(data=assoc[ -log10(assoc$P+0.0000001) > 2,],aes(x=position/1000, y=-log10(P)),size=4, color="red")+
  ylab("−log10 (P value)")+
  xlab("Genomic Postion (kb)")+
  ylim(0,5)+
  theme(axis.text.x = element_text(color = "black",face = "bold", size=18), 
        axis.text.y = element_text(color = "black",face = "bold",size=18),
        axis.title.y = element_text(color = "black",face = "bold",size=20),
        axis.title.x = element_text(color = "black",face = "bold",size=20))

ggsave(plot=p1,width = 18,height = 9, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/assoc5.pdf", useDingbats=FALSE )


#Figure for Pre-amp dNTP optimization:
dntp <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/dNTP_optimization.txt", header = TRUE)
dntp
mdntp <- melt(dntp)

ggplot(mdntp, aes(x=variable, y=value,group=as.factor(Input_EBV_copy.uL)))+geom_bar(stat="identity",position = "dodge")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        plot.background =element_blank())+theme_bw()+scale_fill_gradient()



pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/dNTPOptimization0.pdf",width = 8,height = 4)
ggplot(mdntp, aes(x=variable, y=value, group=as.factor(Input_EBV_copy.uL),fill=as.factor(Input_EBV_copy.uL)))+geom_bar(stat="identity",position = "dodge")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        plot.background =element_blank())+theme_bw()
dev.off()

#Figure for Incubation optimization:
incdata <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/Incubation_optimization.txt",header = TRUE)
incdata

pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/incubationOptimization0.pdf",width = 8,height = 4)
ggplot(incdata, aes(x=Incubation, y=EBVincrease,group=as.factor(Denaturation.Buffer),fill=as.factor(Denaturation.Buffer)))+
  geom_bar(stat="identity",position = "dodge")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        plot.background =element_blank())+theme_bw()+scale_fill_manual(values=c("#E69F00", "#56B4E9"))
dev.off()

#cat my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified-PART1_filtered-OnlyNonSynVariants_Effects.txt |tr -s " "|tr " " '\t' > 1.txt
vareff <- read.delim("/Users/yasinkaymaz/Documents/EBV/ourbatch2/Secondround/test/Filtered-OnlyNonSynVariants_Effects_Fixed.txt",header = FALSE)
head(vareff)
vareff <- droplevels(vareff[!(vareff$V11 == "Z"),])
mat <- confusionMatrix(vareff$V10,vareff$V11)

library(corrplot)
#pheatmap(mat$table,show_rownames = TRUE,cluster_rows = FALSE,cluster_cols=FALSE, legend = FALSE,	cellwidth = 15, cellheight = 15,
#         border_color='white',display_numbers = mat$table)


svg(filename = "/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/nonSynAAtable.svg", width = 10,height = 10)
#corrplot(mat$table, method="circle",is.corr = FALSE, type="upper", order="hclust")
corrplot(mat$table, method="circle",is.corr = FALSE,tl.srt = 0)
dev.off()


#Gene Function Table:
gf <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/Gene.Functions.txt",header = TRUE)
head(gf)

gf[!(gf$Gene %in% assoc$Gene),]
rownames(gf) <- gf$Gene

data <- read.delim("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/dnds_2.txt", row.names = 1, header = TRUE)
head(data)
subdata <- data[,c("gene","SynPerGenomePerKb","NonsynPerGenomePerKb")]
head(subdata)

synCat <- NULL
for (c in levels(gf$Category)){
  print(c)
  cgenes <- rownames(gf[gf$Category == c,])
  print(subdata[cgenes,])
  categ <- data.frame(Category=c,
                      Syn=mean(subdata[cgenes,"SynPerGenomePerKb"]),
                      Nyn=mean(subdata[cgenes,"NonsynPerGenomePerKb"]))
  print(categ)
  synCat <- rbind(synCat, categ)
  
}
rownames(synCat) <- synCat$Category 
#synCat <- synCat[,-c(1)]
msynCat <- melt(synCat)
synplot <- ggplot(msynCat,aes(x=Category, y=value, fill=variable))+geom_bar(stat="identity", position=position_dodge())+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  ylab("Average Variant (per Gene)")+
  xlab("Functional Categories")

ggsave(plot=synplot,width = 8,height = 6, dpi=200, filename="/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/synCategories.pdf", useDingbats=FALSE )

#Epitopes table:
epp <- read.table("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/Epitopes/Epitopes_Table.txt", header = 1)
head(epp)


adist("RPFFHPVGEADYFEY", epp[(epp$EpitopeSequence == "RPFFHPVGEADYFEY"),]$SampleEpitopeSequence)
stringdist::stringsim("RPFFHPVGEADYFEY", as.character(epp[(epp$EpitopeSequence == "RPFFHPVGEADYFEY"),]$SampleEpitopeSequence))

refpep <- "QVADVVRAPGVPAMQPQYF"
testpeps <- c("PWKTLGHYNDOQLTLMNTR","TTTQVAXXXXXXXXXXXX","QVAXXXXXXXXXXXXXXXX","QVADVVRAPXXXXXXXQXX","QVADVVRAPXXXXXXXXXX","QVADVVRAPGVPAMXXXXX","QVADVVRAPGVPAMQPXXX")

for (pep in testpeps){
  newpep <- gsub("X","",pep)
  print(paste(newpep, stringdist::stringsim(refpep, as.character(newpep),method='hamming'), sep=" "))
}


#distmat <- adist(refepp$EpitopeSequence, refepp$SampleEpitopeSequence)
#pheatmap(distmat,cluster_rows = FALSE,cluster_cols = FALSE)
#corrplot(distmat, method="circle",is.corr = FALSE,tl.srt = 0)

#xxx<-head(epp,200)
xxx <- epp

LenPolymorhicPeps <- c("VQPPQLTLQV","RKCCRAKFKQLLQHYR","PHDITYPYTARNIRDAACRAV","APGPGPQPLRESIVCYFM","LDLFGQLTPHTKAVYQPRGA")
RepeatSamples <- c("Raji_new_Rep","eBL-Tumor-0022_Rep", "eBL-Tumor-0024_Rep", "eBL-Tumor-0026_Rep", "eBL-Tumor-0037_Rep","eBL-Tumor-0038_Rep")
BadEpitopes <- c("ILRQLLTGGVKKGRPSLKLQ","QAPTEYTRERRGVGPMPPT","IVTDFSVIK","RLRAEAQVK","SVRDRLARL","GPWVPEQWMFQGAPPSQGTP","QVADVVRAPGVPAMQPQYF","HRCQAIRKK","FIEFVGWLCKKDHTHIREWF","TYSAGIVQI","AVFDRKSDAK","VEITPYKPTW","EDLPCIVSRGGPKVKRPPIF","LEKARGSTY","KRPPIFIRRL","FLRGRAYGL","VSFIEFVGW","RRFPLDLR","RYSIFFDY","VYFQDVFGTMWCHHA","AYSSWMYSY","RQAIRDRRRNPASRR","RRARSLSAERY")
BadSamples <- c("eBL-Tumor-0003", "eBL-Plasma-0009")

xxx <- xxx[!(xxx$EpitopeSequence %in% LenPolymorhicPeps),]
xxx <- xxx[!(xxx$EpitopeSequence %in% BadEpitopes),]

xxx <- xxx[!(xxx$SampleName %in% RepeatSamples),]
xxx <- xxx[!(xxx$SampleName %in% BadSamples),]

#Calculate the similarity between Sample and reference Epitope sequences
Simi=NULL
for (i in 1:length(xxx[,1])){
  #samplepep <- gsub("X","",as.character(xxx[i,]$SampleEpitopeSequence))
  #Simi <- rbind(Simi,stringdist::stringsim(as.character(xxx[i,]$EpitopeSequence), as.character(samplepep),method='hamming'))
  Simi <- rbind(Simi,stringdist::stringsim(as.character(xxx[i,]$EpitopeSequence), as.character(xxx[i,]$SampleEpitopeSequence)))
}
xxx$Similarity <- Simi

refepp<- xxx[(xxx$SampleID == "Type1_Ref"),]

write.table(xxx, file="~/Dropbox/Papers/EBV_project/workspace/data/EpitopesTable_withSimilarities.txt",sep="\t",quote = FALSE,row.names = FALSE)
simMat <- reshape2::dcast(xxx, EpitopeSequence ~ SampleName )


#simMat[is.na(simMat)] <- 0
rownames(simMat) <- simMat$EpitopeSequence
simMat <- simMat[,-c(1)]

colannot <- NULL
for (sample in colnames(simMat)){
  colannot <- rbind(colannot,unique(epp[(epp$SampleName == sample),c("SampleID","ViralSubtype", "SampleGroup")]))
}

rownames(colannot) <- colannot$SampleID
colannot <- colannot[,-c(1)]

rowannot <- refepp[,c("AntigenName","EpitopeType")]
rownames(rowannot) <- refepp$EpitopeSequence

SampleGroup        <- c("purple", "orange", "green")
names(SampleGroup) <- c("eBL", "HealthyControl","CultureBL")

ViralSubtype <- c("red", "blue")
names(ViralSubtype) <- c("Type1", "Type2")

EpitopeType <- c("cyan", "orchid1")
names(EpitopeType) <- c("CD4_T_cell_epitope","CD8_T_cell_epitope")

AntigenName <- c("hotpink","brown","limegreen","navy","firebrick1","gold","dodgerblue","gray87")
names(AntigenName) <- c("BLLF1", "BZLF1", "EBNA1","EBNA2", "EBNA3A", "EBNA3B", "EBNA3C", "LMP1")
simMat <- cbind(simMat, 
                Daudi_new =simMat$Daudi_new)
colnames(simMat)
simMat <- simMat[,-c(1)]

reorderedCols <- rownames(colannot[with(colannot, order(SampleGroup, ViralSubtype)),])
simMat <- simMat[,reorderedCols]

colannot_colors <- list(SampleGroup = SampleGroup,
                        ViralSubtype = ViralSubtype,
                        EpitopeType=EpitopeType,
                        AntigenName=AntigenName)

out <- pheatmap(simMat,border_color='white',cellwidth = 8,cellheight = 10,cluster_rows = FALSE, cluster_cols = FALSE,annotation_row = rowannot,annotation_col = colannot,show_colnames = TRUE,annotation_colors = colannot_colors)


pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EpitopePool-Hamming_TrimX2.pdf",width = 25,height = 25,onefile = FALSE)
out
dev.off()


Conservation=NULL
for (pep in rownames(rowannot)){
line <- simMat[pep,]
covredpeps <- length(line[line == 1])*100/length(line[line != 0])
Conservation <- rbind(Conservation, pep=covredpeps)
}
rownames(Conservation) <- rownames(rowannot)
colnames(Conservation) <- "ConservationLevel"
Conservation <- as.data.frame(Conservation)
Conservation$Epitopes <- rownames(Conservation)
Conservation <- Conservation[rownames(simMat[out$tree_row[["order"]],]),]
Conservation$Epitopes <- factor(Conservation$Epitopes, levels=rownames(simMat[out$tree_row[["order"]],]))

library(ggplot2)
pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EpitopePool-Hamming_TrimX2_cons.pdf",width = 25,height = 5,onefile = FALSE)
ggplot(Conservation, aes(x=Epitopes,y=ConservationLevel))+geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  scale_x_discrete(position = "top")+  scale_y_reverse() +
  theme(axis.text.y = element_text(angle = 90, hjust = 1,vjust=0.5))
dev.off()


#Alternative epitopes plot
library(pheatmap)
out <- pheatmap(simMat[rownames(rowannot),],border_color='white',cellwidth = 8,cellheight = 10,cluster_rows = FALSE, cluster_cols = FALSE,annotation_row = rowannot,annotation_col = colannot,show_colnames = TRUE,annotation_colors = colannot_colors)

head(colannot)

t1.genomes <- rownames(colannot[which(colannot$ViralSubtype == "Type1"),])
t2.genomes <- rownames(colannot[which(colannot$ViralSubtype == "Type2"),])
ebl.genomes <- rownames(colannot[which(colannot$SampleGroup == "eBL"),])
hc.genomes <- rownames(colannot[which(colannot$SampleGroup == "HealthyControl"),])

Conservations <- data.frame(Epitopes = rownames(rowannot), rowannot)
simMat <- simMat[rownames(rowannot),]
Conservations$Type1 <- apply(simMat[,t1.genomes], 1, function(x) length(x[x==1])*100/length(x[x !=0])  )
Conservations$Type2 <- apply(simMat[,t2.genomes], 1, function(x) length(x[x==1])*100/length(x[x !=0])  )
Conservations$eBL <- apply(simMat[,ebl.genomes], 1, function(x) length(x[x==1])*100/length(x[x !=0])  )
Conservations$HC <- apply(simMat[,hc.genomes], 1, function(x) length(x[x==1])*100/length(x[x !=0])  )

Conservations$Epitopes <- factor(Conservations$Epitopes, levels=rownames(rowannot))

Conservations <- reshape::melt(Conservations)
tmp <- simMat[,t2.genomes]

library(ggplot2)
pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EpitopePool-Hamming_TrimX2_cons_all.pdf",width = 40,height = 5,onefile = FALSE)
  ggplot(Conservations, aes(x=Epitopes,y=value))+
    geom_bar(aes(fill=variable),stat="identity")+
    #facet_wrap(~variable)+
    facet_grid(~variable)+
    scale_fill_manual(values =c("red", "blue","purple", "orange"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
    scale_x_discrete(position = "top")+  scale_y_reverse() +
    theme(axis.text.y = element_text(angle = 90, hjust = 1,vjust=0.5))
dev.off()
  


pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/EpitopePool-Hamming_TrimX2_test.pdf",width = 25,height = 25,onefile = FALSE)
out
dev.off()





data <- read.table("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/data/sample_characteristics.txt",header=TRUE,row.names = 1)
ggplot(data,aes(x=reorder(SampleID,Age),y=Age,color=Viral_Type))+geom_point(aes(shape=Diagnose),size=3)

data.e <- droplevels(data[which(data$Diagnose == "eBL"),])
ggplot(data.e,aes(x=reorder(SampleID,Age),y=Age,color=Viral_Type))+
  geom_point(aes(shape=Gender),size=3)+ 
  geom_vline(xintercept = round(length(rownames(data.e))/2), linetype="dotted")+
  geom_hline(yintercept = median(data.e$Age), linetype="dotted")+ 
  scale_y_continuous(breaks=seq(0,15,1))+
  theme(axis.text.x=element_blank())+
  xlab("Individuals")

data.h <- droplevels(data[which(data$Diagnose == "HealthyControl"),])
ggplot(data.h,aes(x=reorder(SampleID,Age),y=Age,color=Viral_Type))+
  geom_point(aes(shape=Gender),size=3)+ 
  geom_vline(xintercept = round(length(rownames(data.h))/2), linetype="dotted")+
  geom_hline(yintercept = median(data.h$Age), linetype="dotted")+ 
  scale_y_continuous(breaks=seq(0,15,1))+
  theme(axis.text.x=element_blank())+
  xlab("Individuals")


#PCoA
library(vcfR)
library(dartR)
require(gridExtra)
library(ggplot2)
library(plotly)
library(tidyverse)
source(here::here("PapeRCodes/functions.R"))
pops <- read.delim(here::here("workspace/data/pop_info.txt"),row.names = 1, header = T)
vcf.ebv <- read.vcfR(here::here("workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_repeatFiltered.vcf"))
gl <- vcfR2genlight(vcf.ebv)
gl <- gl.filter.callrate(gl, method = "loc",threshold = 0.25)
gl <- gl.filter.callrate(gl, method = "ind",threshold = 0.25)
gl$pop <- pops[indNames(gl),2]
type <- pops[indNames(gl),3]
samplelables <- pops[indNames(gl),4]
pc <- gl.pcoa(gl,nfactors = 20)

#Run PCOA analysis for all variants


source(here::here("PapeRCodes/functions.R"))
gl.pcoa.plot(pc, gl, type = samplelables, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "")


pcap <- gl.pcoa.plot(pc, gl, type = samplelables, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "")
pcap <- pcap + xlim(c(-25,25)) + ylim(c(-12.5,12.5))
ggsave(plot=pcap,width = 8,height = 6, dpi=200, filename="~/Dropbox/Papers/EBV_project/workspace/Figure_temps/PCoA1.pdf", useDingbats=FALSE )


pcap <- gl.pcoa.plot(pc, gl, type = type, xaxis=2, yaxis=3, ellipse = F, labels = "pop", title = "")
pcap <- pcap + xlim(c(-15,15)) + ylim(c(-15,15))
ggsave(plot=pcap,width = 8,height = 6, dpi=200, filename="~/Dropbox/Papers/EBV_project/workspace/Figure_temps/PCoA2.pdf", useDingbats=FALSE )


p34 <- gl.pcoa.plot(pc, gl, type = samplelables, xaxis=3, yaxis=4, ellipse = F, labels = "pop", title = "for All variants")
p56 <- gl.pcoa.plot(pc, gl, type = samplelables, xaxis=5, yaxis=6, ellipse = F, labels = "pop", title = "for All variants")
p78 <- gl.pcoa.plot(pc, gl, type = samplelables, xaxis=7, yaxis=8, ellipse = F, labels = "pop", title = "for All variants")

ggsave(plot=grid.arrange(p12,p23,p45, ncol=3),width = 18,height = 5, dpi=200, filename="~/Dropbox/Papers/EBV_project/workspace/Figure_temps/PCoA2.supp1.pdf", useDingbats=FALSE )

#Samples processing table:
data <- read.delim("~/Dropbox/Papers/EBV_project/workspace/data/SamplesProcessing.txt", header=T)
head(data)
table(data[,c("Amplification","Sample.Source")])


data <- read.delim("~/Dropbox/EBV_Capture/MyBait/Mild_Balanced_Boosted_Probes.txt", header = F, row.names = 2)
head(data)
hist(data$V4)
