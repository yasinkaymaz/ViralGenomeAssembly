#This is a script in which I explore the viral genomic variants and population structures.
library(vcfR)
library(dartR)
require(gridExtra)
library(ggplot2)
library(plotly)
library(tidyverse)
source(here::here("PapeRCodes/functions.R"))


pops <- read.delim(here::here("workspace/data/pop7.txt"),row.names = 1, header = T)

summary(pops)

vcf.ebv <- read.vcfR(here::here(("workspace/data/my_ICed_sequence.aln.Sub_for_Assoc_modified.vcf")))
gl <- vcfR2genlight(vcf.ebv)
gl <- gl.filter.callrate(gl, method = "ind",threshold = .25, recalc = TRUE)
gl$pop <- pops[indNames(gl),2]
type <- pops[indNames(gl),1]

pc <- gl.pcoa(gl,nfactors = 20)

df <- as.data.frame(pc$loadings)
ggplot(df, aes(x=as.numeric(gl@position),y=Axis1, group = 1))+
  geom_point(size=.1)+
  labs(x="genomic positions")+
  ggtitle("PCA Loadings from Axis 1")

p <- gl.pcoa.plot(pc, gl, type = type, xaxis=1, yaxis=3, ellipse = F, labels = "pop", title = "for All variants")
ggplotly(p)

gl
gl.filt <- gl.filter.callrate(gl, method = "loc",threshold = 0.25)
gl.filt <- gl.filter.callrate(gl.filt, method = "ind",threshold = 0.25)



Fstats.ebvtype <- NAM::Fst(as.matrix(gl.filt),fam=as.numeric(pops[indNames(gl.filt),]$EBV_Type) )
Fstats.phenotype <- NAM::Fst(as.matrix(gl.filt),fam=as.numeric(pops[indNames(gl.filt),]$Condition) )
plot(Fstats.ebvtype)
plot(Fstats.phenotype)

fstdata <- data.frame(loc=as.numeric(gl.filt@position), Fst.EBVTypes=Fstats.ebvtype$fst, Fst.pheno = Fstats.phenotype$fst )
ggplot(fstdata, aes(loc, Fst.EBVTypes))+ geom_point(color="red",size=0.2) + ylim(0,1) 
ggplot(fstdata, aes(loc, Fst.pheno))+ geom_point(color="blue",size=0.2) + ylim(0,1)






vcf.ebv.m <- read.vcfR(here::here("workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_filtered.vcf"))
gl.m <- vcfR2genlight(vcf.ebv.m)
gl.m <- gl.filter.callrate(gl.m, method = "ind",threshold = .25, recalc = TRUE)
gl.m$pop <- pops[indNames(gl.m),2]
type.m <- pops[indNames(gl.m),1]

#Run PCoA
pc.m <- gl.pcoa(gl.m)



pdf("workspace/Figure_temps/PCoA.2.pdf",width = 30,height = 8)
p12 <- gl.pcoa.plot(pc, gl, type = type, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "for All variants")
p23 <- gl.pcoa.plot(pc, gl, type = type, xaxis=2, yaxis=3, ellipse = F, labels = "pop", title = "for All variants")
p45 <- gl.pcoa.plot(pc, gl, type = type, xaxis=4, yaxis=5, ellipse = F, labels = "pop", title = "for All variants")
grid.arrange(p12,p23,p45, ncol=3)

p12 <- gl.pcoa.plot(pc.m, gl.m, type = type.m, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "for non-type specific variants (EBNA regions excluded).")
p23 <- gl.pcoa.plot(pc.m, gl.m, type = type.m, xaxis=2, yaxis=3, ellipse = F, labels = "pop", title = "for non-type specific variants (EBNA regions excluded).")
p45 <- gl.pcoa.plot(pc.m, gl.m, type = type.m, xaxis=4, yaxis=5, ellipse = F, labels = "pop", title = "for non-type specific variants (EBNA regions excluded).")
grid.arrange(p12,p23,p45, ncol=3)

dev.off()


pops <- read.delim("workspace/data/pop7.txt",row.names = 1, header = T)

vcf.ebv <- read.vcfR("workspace/data/my_ICed_sequence.aln.Sub_for_Assoc_modified.vcf")
gl <- vcfR2genlight(vcf.ebv)
gl$pop <- pops[indNames(gl),2]
table(pop(gl))

strata(gl) <- as.data.frame(pops[indNames(gl),1])

pc <- gl.pcoa(gl,nfactors = 20)

gl.pcoa.plot(pc, gl, xaxis=1, yaxis=2, ellipse = F, labels = "pop",title ="", type = type)


gl.filt <- gl.filter.callrate(gl, method = "ind",threshold = .25, recalc = TRUE)
gl.filt$pop <- pops[indNames(gl.filt),2]
#gl.filt <- gl.filter.callrate(gl.filt, method = "loc",threshold = 0.50, recalc = TRUE)
gl.filt

pc <- gl.pcoa(gl.filt)
gl.pcoa.plot(pc, gl, xaxis=1, yaxis=2, ellipse = F, labels = "pop")
p <- gl.pcoa.plot(pc, gl.filt, xaxis=1, yaxis=2, ellipse = F, labels = "interactive")
ggplotly(p)

# library(StAMPP)
# pwfst <-stamppFst(gl.filt, nboots=100, percent=95, nclusters=1)
d <- round(dim(as.matrix(gl))[2]/100)
for (i in 1:d){
  print((i+1)*100)
  
  }


Fstats.ebvtype <- NAM::Fst(as.matrix(gl),fam=as.numeric(pops$EBV_Type) )
Fstats.phenotype <- NAM::Fst(as.matrix(gl),fam=as.numeric(pops$Condition) )
plot(Fstats.ebvtype)
plot(Fstats.phenotype)
fstdata <- data.frame(loc=as.numeric(gl@position), Fst.EBVTypes=Fstats.ebvtype$fst, Fst.pheno = Fstats.phenotype$fst )
head(fstdata)
ggplot(fstdata, aes(loc, Fst.EBVTypes))+ geom_point(size=0.2) + ylim(0,1) 
ggplot(fstdata, aes(loc, Fst.pheno))+ geom_point(size=0.2) + ylim(0,1)

#Subset types
#Type 1
pop.t1 <- pops[which(pops$EBV_Type == "type_1"),]
type_1.genomes <- rownames(pop.t1)

Fstats.t1 <- NAM::Fst(as.matrix(gl)[type_1.genomes,],fam=as.numeric(pop.t1$Condition) )
fstdata.t1 <- data.frame(loc=as.numeric(gl@position), Fst=Fstats.t1$fst)
ggplot(fstdata.t1, aes(loc, Fst))+ geom_point(size=0.2) + ylim(0,1) + ggtitle("Differentiation between eBL and Healthy Controls within Type 1")

#Type 2
pop.t2 <- pops[which(pops$EBV_Type == "type_2"),]
type_2.genomes <- rownames(pop.t2)

Fstats.t2 <- NAM::Fst(as.matrix(gl)[type_2.genomes,],fam=as.numeric(pop.t2$Condition) )
fstdata.t2 <- data.frame(loc=as.numeric(gl@position), Fst=Fstats.t2$fst)
ggplot(fstdata.t2, aes(loc, Fst))+ geom_point(size=0.2) + ylim(0,1)+ ggtitle("Differentiation between eBL and Healthy Controls within Type 2")


bindata <- read.delim("workspace/data/out.windowed.weir.fst",header = T)
head(bindata)
ggplot(bindata, aes(BIN_START, WEIGHTED_FST))+ geom_line(size=0.2) + ylim(0,1)



#REdo
source(here::here("PapeRCodes/functions.R"))

vcf.ebv <- read.vcfR("/Users/yasinkaymaz/Documents/EBV/ourbatch2/Secondround/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified.vcf")
gl <- vcfR2genlight(vcf.ebv)
gl <- gl.filter.maf(gl,v = 3)

gl <- gl.filter.callrate(gl, method = "loc",threshold = 0.25, recalc = TRUE,v=3)
gl <- gl.filter.callrate(gl, method = "ind",threshold = 0.25, recalc = TRUE,v=3)
gl

gl$pop <- pops[indNames(gl),2]
type <- pops[indNames(gl),3]
#Run PCOA analysis for all variants
pc <- gl.pcoa(gl,nfactors = 20)


data <- read.delim(here::here("workspace/data/pairwisedistances.ordered.txt"),header=F)
data1 <- read.delim(here::here("workspace/data/Type1.pairwisedistances.ordered.txt"),header=F)
data2 <- read.delim(here::here("workspace/data/Type2.pairwisedistances.ordered.txt"),header=F)

data1$Group <- "Type 1"
data2$Group <- "Type 2"

d1d2 <- rbind(data1,data2)

ggplot(d1d2, aes(x=reorder(V2,V1), y=V3, color=Group)) + 
  geom_errorbar(aes(ymin=V3-V4, ymax=V3+V4), width=.1) +
  geom_point()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("EBV Genes")+
  ylab("d (Kimura 2 Parameter)")+
  ggtitle("Gene by gene genetic distances of EBV genomes")


data <- cbind(data1,data2)
head(data)
colnames(data) <- c("col1","gene1", "d_t1", "err_t1", "group1","col2","gene2", "d_t2", "err_t2", "group2")
data$delta_d <- data$d_t1 - data$d_t2

ggplot(data, aes(x=reorder(gene1,col1), y=delta_d)) +
  geom_point(shape=23, fill="orange")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("EBV Genes")+
  ylab("Delta d (Type1 - Type2")+
  ylim(-0.0025, 0.02)+
  ggtitle("Difference between the divergence of two types")


source(here::here("PapeRCodes/functions.R"))

p <- gl.pcoa.plot(pc, gl, type = type, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "for All variants")
ggplotly(p)


gl.t1 <- gl.keep.pop(gl, pop.list=c("Type1"),recalc = TRUE)
Fstats.ebvtype <- NAM::Fst(as.matrix(gl.t1),fam=as.numeric(pops[indNames(gl.t1),]$EBV_Type) )
Fstats.phenotype <- NAM::Fst(as.matrix(gl.t1),fam=as.numeric(pops[indNames(gl.t1),]$Condition) )
fstdata <- data.frame(loc=as.numeric(gl.t1@position), Fst.EBVTypes=Fstats.ebvtype$fst, Fst.pheno = Fstats.phenotype$fst )

gl.t2 <- gl.keep.pop(gl, pop.list=c("Type2"),recalc = TRUE)
Fstats.ebvtype <- NAM::Fst(as.matrix(gl.t2),fam=as.numeric(pops[indNames(gl.t2),]$EBV_Type) )
Fstats.phenotype <- NAM::Fst(as.matrix(gl.t2),fam=as.numeric(pops[indNames(gl.t2),]$Condition) )
fstdata <- data.frame(loc=as.numeric(gl.t2@position), Fst.EBVTypes=Fstats.ebvtype$fst, Fst.pheno = Fstats.phenotype$fst )


window = 500; newfstdata <- NULL; fstdata.loc <- NULL; j=0

for(i in 0:max(unique(round(fstdata$loc/window))) ){
  fstdata.loc <- fstdata[which((fstdata$loc >= j*window ) & (fstdata$loc < i*window+window) ),]
  if(dim(fstdata.loc)[1] == 0){
    emptydf <- data.frame(c(i*window, i*window+window),c(0,0),c(0,0),c(0,0) )
    colnames(emptydf) <- colnames(newfstdata)
    newfstdata <- rbind(newfstdata,emptydf)    
  }else{
    fstdata.loc$MeanFst <- mean(fstdata.loc$Fst.pheno)
    newfstdata <- rbind(newfstdata,fstdata.loc)
  }
  j <- i+1
}

#newfstdata$Gene <- ifelse( ((newfstdata$loc >= 80000 ) & (newfstdata$loc < 90000)),"EBNA","NA")

vcf.ebv <- read.vcfR(here::here("workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_repeatFiltered.vcf"))

gl <- vcfR2genlight(vcf.ebv)
gl$pop <- pops[indNames(gl),2]

pop(gl)
gl.t2 <- gl.keep.pop(gl.t2, pop.list=c("Type2"),recalc = TRUE)
gl.t2 <- gl.filter.callrate(gl.t2, method = "loc",threshold = 0.25)
gl.t2 <- gl.filter.callrate(gl.t2, method = "ind",threshold = 0.25)

#Import genome type and patient info
gl.t2$pop <- pops[indNames(gl.t2),2]
type <- pops[indNames(gl.t2),3]
samplelabels <- pops[indNames(gl.t2),4]
method <- pops[indNames(gl.t2),1]

#Run PCOA analysis for all variants
pc <- gl.pcoa(gl.t2,nfactors = 20)
p <- gl.pcoa.plot(pc, gl.t2, type = samplelabels, xaxis=1, yaxis=2, ellipse = F, labels = "pop", title = "for All variants")
ggplotly(p)


#Average window of PCA loadings
rect <- data.frame(xmin=c(36216,79955,83065,86083,95662), xmax=c(37679,82877,85959,89135,97587), ymin=-Inf, ymax=Inf,gene=c("EBNA2","EBNA3A","EBNA3B","EBNA3C","EBNA1"))

gl$position

gl <- vcfR2genlight(vcf.ebv)

#gl <- gl.filter.maf(gl,v = 3)

gl <- gl.filter.callrate(gl, method = "loc",threshold = 0.25, recalc = TRUE,v=3)
gl <- gl.filter.callrate(gl, method = "ind",threshold = 0.25, recalc = TRUE,v=3)
gl

gl$pop <- pops[indNames(gl),2]
type <- pops[indNames(gl),3]
#Run PCOA analysis for all variants
pc <- gl.pcoa(gl,nfactors = 20)


iii=0
plots <- list() 
pcaplots <- list()
for (Axis in c("Axis1","Axis2","Axis3","Axis4","Axis5","Axis6")){
  #print(genome)
  iii=iii+1

  df <- data.frame(loc = gl$position, Load=abs(pc$loadings[,Axis])) 
  head(df)
  dim(df)
  
  window = 1000; newdf <- NULL; df.loc <- NULL; j=0
  for(i in 0:max(unique(round(df$loc/window))) ){
    df.loc <- df[which((df$loc >= j*window ) & (df$loc < i*window+window) ),]
    print(df.loc)
    if(dim(df.loc)[1] == 0){
      emptydf <- data.frame(c(i*window, i*window+window),c(0,0),c(0,0) )
      colnames(emptydf) <- colnames(newdf)
      newdf <- rbind(newdf,emptydf)    
    }else{
      df.loc$MeanLoad <- mean(df.loc$Load)
      newdf <- rbind(newdf,df.loc)
    }
    j <- i+1
  }
  
  p1 <- ggplot(newdf, aes(loc/1000, MeanLoad))+ geom_line(color="darkgreen",size=0.8) +
 #  ylim(0,0.03)+
    geom_rect(data=rect/1000, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              color="transparent", fill="orange", alpha=0.3,
              inherit.aes = FALSE, show.legend = T)+
    labs(x="Genomic Position (Kb)",y=paste("mean PC",Axis,"Loadings",sep=" ") )+
    geom_text(data=rect, aes(x=xmin/1000+((xmax-xmin)/1000)/2, y=0.025, label=gene), size=4,angle=45)

  pca.p <- gl.pcoa.plot(pc, gl, type = type, xaxis=iii, yaxis=iii+1, ellipse = F, labels = "pop", title = "for All variants")
  #print(p1)
  plots[[iii]] <- p1
  pcaplots[[iii]] <- pca.p
  
  
}


pdf("/Users/yasinkaymaz/Dropbox/Papers/EBV_project/workspace/Figure_temps/MeanPCALoadings5_type2.pdf",width = 18,height = 30)
multiplot(plotlist = plots)
multiplot(plotlist = pcaplots)
dev.off()



source("~/Dropbox/codes/ViralGenomeAssembly/PapeRCodes/functions.R")

pops <- read.delim("~/Dropbox/codes/ViralGenomeAssembly/workspace/data/pop_info.txt",row.names = 1, header = T)
summary(pops)
vcf.ebv <- read.vcfR("~/Dropbox/codes/ViralGenomeAssembly/workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_repeatFiltered.vcf")

gl <- vcfR2genlight(vcf.ebv)


#Import genome type and patient info
gl$pop <- pops[indNames(gl),2]

gl <- gl.keep.pop(gl, pop.list=c("Type2"),recalc = TRUE)

gl <- gl.filter.callrate(gl, method = "loc", threshold = 0.05)
gl <- gl.filter.callrate(gl, method = "ind", threshold = 0.25)

type <- pops[indNames(gl),3]
samplelabels <- pops[indNames(gl),4]
method <- pops[indNames(gl),1]
#Run PCOA analysis for all variants
pc <- gl.pcoa(gl,nfactors = 20)

p <- gl.pcoa.plot(pc, gl, type = type, xaxis=4, yaxis=5, ellipse = F, labels = "pop", title = "for All variants")
ggplotly(p)



library(SNPRelate)
vcf.fn <- system.file("~/Dropbox/codes/ViralGenomeAssembly/workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_repeatFiltered.vcf", package="SNPRelate")

snpgdsVCF2GDS(vcf.fn, "my.gds")


