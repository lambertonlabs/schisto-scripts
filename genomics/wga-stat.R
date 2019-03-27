#r script for performing basic sumamry stats on schisto WGS
#if possible, use VCFtools and then ggplot to generate your nice manhattan plots
#adegenet doesn't tolerate large SNP sets very well, and popgenome a) won't analyse multiple chromosomes and b) is very finicky about inputs


#install external pkgs
install.packages("PopGenome")
install.packages("adegenet")
install.packages("vcfR")
install.packages("poppr")

#load basically all popular popgen modules
library(PopGenome)
library(adegenet)
library(vcfR)
library(poppr)

#load data and convert to object

#setwd to dir containing vcfs
setwd("~/Projects/Schisto/GENOMIC/VCF/popgenome/")

snp <- read.vcfR("CF.vcf")
x <- vcfR2genind(snp)



