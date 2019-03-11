##################################################
# RDA-based genome scan approach climate dataset #
##################################################

# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

library(adegenet)
library(robust)
library(qvalue)
library(ggplot2)
library(psych)
library(vegan)
library(adegenet)
library(randomcoloR)

file_path <- "/Users/umer/Documents/Results_Rad/World_Clim/"
setwd(file_path)

# To extract from plink.map file
#extract.PLINKmap("plink.map", x = NULL)

# after conversion with plink from .ped and .map (both obtained from populations) read .raw file
plink <- read.PLINK("plink.raw", map.file = NULL, quiet = FALSE, chunkSize = 1000,
                    parallel = require("parallel"), n.cores = NULL)
#/// GENLIGHT OBJECT /////////

# Turn into a matrix
gen <- as.matrix(plink)

# Remove NAs
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#sum(is.na(gen.imp))

# Get environmental data downloaded from WorldClim
env <- read.csv("environmantal_variables_20181106_4.5km.csv")
env$Population <- factor(env$Population, levels = c("Fortrose", "Auchmithie", "Crail", "Kildonan Castle", "Tynemouth",
                                                    "Whitby","Tenby","Llantwit Major","St Aldehams Head","West Looe",
                                                    "Fowey","Prussia cove","Cudillero","Playa de Xago","Cabo de Penas","Playas de Viodo","Playa Pedrero",
                                                    "San Juan de Gaxtelugatxe","Getaria","San Sebastian"))
# Rename environmental variables 
names(env)[names(env) == 'BIO1'] <-'Annual_Mean_Temperature'
names(env)[names(env) == 'BIO2'] <-'Mean_Diurnal_Range'
names(env)[names(env) == 'BIO3'] <-'Isothermality' 
names(env)[names(env) == 'BIO4'] <-'Temperature_Seasonality'
names(env)[names(env) == 'BIO5'] <-'Max_Temperature_of_Warmest_Month' #1
names(env)[names(env) == 'BIO6'] <-'Min_Temperature_of_Coldest_Month' #1
names(env)[names(env) == 'BIO7'] <-'Temperature_Annual_Range' 
names(env)[names(env) == 'BIO8'] <-'Mean_Temperature_of_Wettest_Quarter' 
names(env)[names(env) == 'BIO9'] <-'Mean_Temperature_of_Driest_Quarter' #1
names(env)[names(env) == 'BIO10'] <-'Mean_Temperature_of_Warmest_Quarter' #1
names(env)[names(env) == 'BIO11'] <-'Mean_Temperature_of_Coldest_Quarter' #1
names(env)[names(env) == 'BIO12'] <-'Annual_Precipitation'
names(env)[names(env) == 'BIO13'] <-'Precipitation_of_Wettest_Month' #
names(env)[names(env) == 'BIO14'] <-'Precipitation_of_Driest_Month' #
names(env)[names(env) == 'BIO15'] <-'Precipitation_Seasonality'
names(env)[names(env) == 'BIO16'] <-'Precipitation_of_Wettest_Quarter' #
names(env)[names(env) == 'BIO17'] <-'Precipitation_of_Driest_Quarter' #
names(env)[names(env) == 'BIO18'] <-'Precipitation_of_Warmest_Quarter' #
names(env)[names(env) == 'BIO19'] <-'Precipitation_of_Coldest_Quarter' #

# Remove first row
env <- env[,-1]
# Get names of individuals
places <- read.csv("Population_number.csv",header=TRUE)
env$individual <- places$Individual
env$individual <- as.character(env$individual)
# sort
env <- env[order(env$individual),]
row.names(env) <- env$individual

#write.csv(gen.imp,"gen_2_12_18.csv")
# Read in data
gen.imp <- read.csv("gen_2_12_18.csv")

# Rename
# gen.imp[,1] contains the names and is called 'X'
str(gen.imp$X)
gen.imp[,1] <- sub("^([^_]*_[^_]*).*", "\\1", gen.imp[,1])
# Sort
gen.imp <- gen.imp[order(gen.imp[,1]),]
# Samples on row names
row.names(gen.imp) <- gen.imp[,1]
gen.imp <- gen.imp[,-1]

# Confirm that genotypes and environmental data are in the same order
gen.imp <- subset(gen.imp[which(row.names(gen.imp) %in% env$individual),])
# Reorder alphabetically
gen.imp <- gen.imp[order(row.names(gen.imp)),]
identical(rownames(gen.imp), env$individual) 

# Check for correlations in data
# |r| > 0.7 
pairs.panels(env[,4:23], scale=T)

# Try with: 'Annual_Mean_Temperature''Mean_Diurnal_Range''Isothermality' 'Temperature_Seasonality'
# 'Mean_Temperature_of_Wettest_Quarter' 'Annual_Precipitation''Precipitation_Seasonality'

# Test run with just these, but go back and check for correlations between env variables
pred <- subset(env, select=c(Long,Annual_Mean_Temperature,#Mean_Diurnal_Range,
                             #Isothermality,Temperature_Seasonality,
                             Mean_Temperature_of_Wettest_Quarter,
                             Annual_Precipitation,Precipitation_Seasonality))
# Check correlations again
pairs.panels(pred, scale=T)

# Run RDA
cabbage.rda <- rda(gen.imp ~ Long + Annual_Mean_Temperature +
                     Mean_Diurnal_Range + Isothermality + Temperature_Seasonality +
                     Mean_Temperature_of_Wettest_Quarter + Annual_Precipitation + Precipitation_Seasonality, data=pred, scale=T)
cabbage.rda

# CHECK WEBPAGE #

geno <- gen.imp

# Remove loci with frequency of minor allele < 0.05
#MAF <- 0.05
#frequencies <- colSums(geno)/(2*nrow(geno))
#maf <- which(frequencies > MAF & frequencies < (1-MAF))
#geno <- geno[,maf]

# Function that returns p-values and q-values for all the loci
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist 
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# Choose the best number of axes for the analyses
ggplot() +
  geom_line(aes(x=c(1:length(cabbage.rda_2$CCA$eig)), y=as.vector(cabbage.rda_2$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(cabbage.rda_2$CCA$eig)), y=as.vector(cabbage.rda_2$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) + ylab("Inertia") +
  theme_bw()

# Take the one just after the plataeu. In this case I would choose K=6
# Use radapt function with K=6
res_rdadapt <- rdadapt(cabbage.rda_2, 5)

# Get the names of the loci from the RDA and put with p- and q-values
rs <- names(cabbage.rda_2$colsum)
# Remove the X at the beginning of the snp names
rs <- gsub("X", "", paste(rs))
# Remove the bit after that isn't needed
rs <- sub("^([^_]*_[^_]*).*", "\\1", rs)
# Combine
snps_p_values_2_12_18 <- cbind(rs,res_rdadapt)

write.table(snps_p_values_2_12_18,"snps_p_values_2_12_18.txt",row.names = F)

# Visualise the results in a Manhattan plot that isn't against the reference (outliers: q.value < 0.1 are coloured in orange)
ggplot() +
  geom_point(aes(x=c(1:length(res_rdadapt[,1])), y=-log10(res_rdadapt[,1])),
             col = "gray83") + geom_point(aes(x=c(1:length(res_rdadapt[,1]))[which(res_rdadapt[,2] < 0.1)],
                                              y=-log10(res_rdadapt[which(res_rdadapt[,2] < 0.1),1])), col = "orange") +
  xlab("SNPs") + ylab("-log10(p.values)") + theme_bw()

# Outlier loci (q.values < 0.1)
outliers <- subset(snps_p_values_2_12_18[which(res_rdadapt[,2] < 0.1),])

# Projection of loci into RDA space
ggplot() +
  geom_point(aes(x=cabbage.rda_2$CCA$v[,1], y=cabbage.rda_2$CCA$v[,2]), col = "gray86") + 
  geom_point(aes(x=cabbage.rda_2$CCA$v[which(res_rdadapt[,2] < 0.1),1], 
                 y=cabbage.rda_2$CCA$v[which(res_rdadapt[,2] < 0.1),2]), col = "orange") + 
  geom_segment(aes(xend=cabbage.rda_2$CCA$biplot[,1]/10, yend=cabbage.rda_2$CCA$biplot[,2]/10, x=0, y=0),
               colour="black", size=0.5, linetype=1,
               arrow=arrow(length = unit(0.02, "npc"))) + 
  geom_text(aes(x=1.2*cabbage.rda_2$CCA$biplot[,1]/10, y=1.2*cabbage.rda_2$CCA$biplot[,2]/10,
                label = colnames(pred[,1:5]))) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_bw()


###################
# Manhattan plots #
###################

library(qqman)

file_path <- "/Users/umer/Documents/Results_Rad/manhattan/"
setwd(file_path)

# all_snps.txt file contains the following columns:
# chr -- the chromosome number
# rs -- the RAD_locus
# ps -- the position in the genome

# Positions of snps
pos <- read.table("all_snps.txt", header = T, sep="\t")
pos$chr_1 <- gsub("C", "", paste(pos$chr))
c<-as.numeric(pos$chr_1)
levels(as.factor(c))

# Snps with p- and q-values contains:
# rs -- the RAD_locus
# p-value
# q-value

snps_p_values <- read.table("snps_p_values.txt", header = T, sep = "")

library(dplyr)
df3 <- left_join(pos, snps_p_values, by = "rs")
df3[50700:50708,1:7]
df3$chr <- as.numeric(df3$chr)

# Outliers from above
outliers <- subset(outliers[which(outliers$ps %in% df3$ps),])
chrps <- outliers$rs
chrps <- chrps[chrps %in% df3$rs]

manhattan(df3, snp = "rs", chr = "chr", bp = "ps", 
          p = "q.values", logp=TRUE, ylab="Effect", ylim=c(0,5), highlight = chrps)
pdf("Manhattan.pdf", useDingbats=FALSE)


## Top outliers
top_snps <- subset(outliers[which(outliers$q.values > 0.098),])

gene_search_top_24 <- subset(df3[which(df3$rs %in% top_snps$rs),])
write.csv(gene_search_top_24,"gene_search_top_24.csv")

gene_names_top_24 <- read.csv("genes_top_24.csv")

single_gene_names_df_top_24 <- subset(gene_names_top_24, !duplicated(name))
single_gene_names_top_24 <- single_gene_names_df_top_24$name

single_gene_names_top_24 <- gsub("ID=", "", paste(single_gene_names_top_24))
single_gene_names_top_24 <- gsub("Parent=", "", paste(single_gene_names_top_24))
single_gene_names_top_24 <- gsub(";", "", paste(single_gene_names_top_24))

check_genes_top_24 <- cbind(single_gene_names_df_top_24,single_gene_names_top_24)

single_gene_names_top_24 <- unique(single_gene_names_top_24)

check_genes_top_24 <- subset(check_genes_top_24, !duplicated(single_gene_names_top_24))

write.csv(check_genes_top_24,"check_genes_top_24.csv")


### I want to match the outliers to the df3 dataframe so that I can obtain a start and end position
# Subset dataframe based on outliers and then write to csv file
gene_search <- subset(df3[which(df3$rs %in% chrps),])
write.csv(gene_search,"gene_search.csv")

file_path <- "/Users/umer/Documents/Results_Rad/bedtools_trials/"
setwd(file_path)
gene_names <- read.csv("genes.csv")

single_gene_names_df <- subset(gene_names, !duplicated(name))
single_gene_names <- single_gene_names_df$name

single_gene_names <- gsub("ID=", "", paste(single_gene_names))
single_gene_names <- gsub("Parent=", "", paste(single_gene_names))
single_gene_names <- gsub(";", "", paste(single_gene_names))

check_genes <- cbind(single_gene_names_df,single_gene_names)

single_gene_names <- unique(single_gene_names)

check_genes <- subset(check_genes, !duplicated(single_gene_names))

write.csv(check_genes,"check_genes.csv")

