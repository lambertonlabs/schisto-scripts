############################################################
## In silico digest 
## adapted from arne jacobs
## christina.faust@gmail.com
## 10 sept 2018
############################################################

 setwd("~/Genomes/Schistosoma_mansoni")

# uncomment to download packages the first time
# source("http://bioconductor.org/biocLite.R") 
# biocLite("ShortRead")

library(SimRAD)
library(ShortRead)

## Read reference genome
ref = ref.DNAseq("schistosoma_mansoni.PRJEA36577.WBPS11.genomic.fa", prop.contigs = 1) #replace with fasta file of whole genome
#ref =  sim.DNAseq(size=10000000, GCfreq=0.433) # simulated genome

######
## defining enzymes for ddRAD 

## PstI-HF enzyme (6 base pairs)
pst_5p1 = "CTGCA"
pst_3p1 = "G"

#MspI enzyme
mspI_5p2 = "C"
mspI_3p2 = "CGG"

#####
#Simulating the digest 
simseq_dig = insilico.digest(ref, pst_5p1, pst_3p1,
                             mspI_5p2, mspI_3p2, verbose=TRUE)

#Simulating the flanking by two different restriction enzymes as it is the case for ddRAD
simseq_sel = adapt.select(simseq_dig, type = "AB+BA", pst_5p1, pst_3p1, mspI_5p2, mspI_3p2)

#Simulating the size selection step (change min and max size accordingly, try different size ranges):
size_simseq = size.select(simseq_sel, min.size = 0, max.size = 100000000, graph=F, verbose=TRUE)
#size_simseq_150_600 = size.select(simseq_sel, min.size = 150, max.size = 600, graph=TRUE, verbose=TRUE)

size_simseq

binneddata <- hist(size_simseq@ranges@width, breaks=seq(0,10000000,by=50), xlim=c(0, 10000))

bins = as.integer(binneddata$breaks)
binsnew = bins[2:200001]

selectiondf = cbind(binsnew, binneddata$counts)
colnames(selectiondf) <- c("fragment length", "count")
head(selectiondf)

write.csv(selectiondf, 'selectiondata.csv')

#bk<-hist((size_simseq_300_600@ranges@width), breaks=length(size_simseq_300_600)/20, plot=FALSE)$breaks
#hist(width(size_simseq_300_600),border="grey75", col="grey75", breaks=bk, main="", xlab="Locus size (bp)", ylab="Number of loci", xlim=c(0, 1000), ylim=c(0, 10000))

#hist(width(size_simseq_300_600),border="red", col="red", add=TRUE, breaks=bk)
#text(mean(c(min.size, max.size)), max(hist(width(size_simseq_300_600), breaks=bk, plot=FALSE)$counts), pos=4, labels=paste(length(size_simseq_300_600), " loci between ", min.size, " and ", max.size, " bp", sep=""), col="red", cex=0.9, font=2)

#writeFasta(simseq_sel, file="simseq_allfragments.fa", mode ="w")

#boxplot(list(width(simseq_sel), width(size_simseq_300_600)))

#hist(unselectedfraglen,
#     main="Unselected fragments",
#     xlim=c(0, 1000000),
#)
#
#unselectedfraglen = simseq_sel@ranges@width
#
#
#selectedfraglen = size_simseq_300_600@ranges@width
#hist(selectedfraglen,
#     xlim=c(0, 10000),
#     breaks=100
#)
