#This is a script in which I explore the viral genomic variants and population structures.
library(vcfR)
library(dartR)
require(gridExtra)
library(ggplot2)
library(plotly)
pops <- read.delim("workspace/data/pop7.txt",row.names = 1, header = T)

vcf.ebv <- read.vcfR("workspace/data/my_ICed_sequence.aln.Sub_for_Assoc_modified.vcf")
gl <- vcfR2genlight(vcf.ebv)
gl <- gl.filter.callrate(gl, method = "ind",threshold = .25, recalc = TRUE)
gl$pop <- pops[indNames(gl),2]
type <- pops[indNames(gl),1]

vcf.ebv.m <- read.vcfR("workspace/data/my_ICed_my_edited_sequence.aln.IDreplaced.OnlyIndividuals-MinusCellLines_sub_modified_filtered.vcf")
gl.m <- vcfR2genlight(vcf.ebv.m)
gl.m <- gl.filter.callrate(gl.m, method = "ind",threshold = .25, recalc = TRUE)
gl.m$pop <- pops[indNames(gl.m),2]
type.m <- pops[indNames(gl.m),1]

#Run PCoA
pc <- gl.pcoa(gl)
pc.m <- gl.pcoa(gl.m)

library(tidyverse)

pc$loadings %>% as.tibble() %>% ggplot(aes())
df <- as.data.frame(pc$loadings)
dim(df)
ggplot(df, aes(x=as.numeric(gl@position),y=Axis2, group = 1))+geom_point(size=.1)

source("~/Dropbox/codes/ViralGenomeAssembly/PapeRCodes/functions.R")
p <- gl.pcoa.plot(pc, gl, type = type, xaxis=1, yaxis=3, ellipse = F, labels = "pop", title = "for All variants")
ggplotly(p)

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

gl.pcoa.plot(pc, gl, xaxis=1, yaxis=2, ellipse = F, labels = "pop")


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
plot(Fstats)
fstdata <- data.frame(loc=as.numeric(gl@position), Fst.EBVTypes=Fstats$fst, Fst.pheno = Fstats.phenotype$fst )
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

