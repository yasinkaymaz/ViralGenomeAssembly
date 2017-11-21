

# Requires 3.2.0 or above!!!

args = commandArgs (T)
SAMPLE_NAME = args[1]
INPUTFILE = args[2]


library(ggplot2)
library(reshape2)
library(Biobase)

haplofreqs <- read.delim(INPUTFILE,header=FALSE,row.names = 1)

out <- strsplit(as.character(haplofreqs$V4),',')
#out <- strsplit(as.character(haplofreqs[which(haplofreqs$V4 != 100),]$V4),',')

#out <- as.data.frame(lapply(out, FUN=function(x) sort(as.numeric(c(unlist(x), rep(0, max(lengths(out))-length(x)))),decreasing= TRUE)))
out <- as.data.frame(lapply(out, FUN=function(x) sort(as.numeric(c(unlist(x), rep(0, max(lengths(out))-length(x)))),decreasing= FALSE)))

col.n <- c()
i=dim(out)[1]
while (i > 0 ){
  col.n <- c(col.n , paste("hap",i,sep = ""))
  i=i-1
  }

out$types <- col.n

out.m <- melt(out)
colnames(out.m) <- c("types","loc", "hapfreq")
head(out.m)

p <- ggplot(out.m,aes(x=as.numeric(loc),y=hapfreq))+
geom_area(aes(colour=types,fill=types),data=out.m,stat="identity",position = "stack")+
ylab("Haplotype Frequency") +
  xlab("Variant Loci") +
  labs(title=SAMPLE_NAME)

pdf(paste(SAMPLE_NAME,"_HaploFreq.pdf",sep=""),width = 20,height = 5)
print(p)
dev.off()

