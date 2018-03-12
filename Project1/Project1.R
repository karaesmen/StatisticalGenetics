wd = "~/Box Sync/Stat Anal of Genetic Data/Project1/"
setwd(wd)

library(genetics)
library(LDheatmap)

####################################################################
#### Load datasets & transform the data to the required format #####

# Gujarati Indian SNP data
gih <- read.table("HapMap_chrom22_20SNP_GIH.txt")

# look at data structure
dim(gih) # n sample = 101 
head(gih, 3)
tail(gih, 3) 

unique(gih$V8) # whole column 8 is replaced with TRUE istead if "T"
gih$V8 <- "T"  # fix the problem with R, reading T as TRUE on column 8
head(gih, 3)

# Mexican SNP data
mex <- read.table("HapMap_chrom22_20SNP_MEX.txt")

# look at data structure
dim(mex)  # n sample = 86
head(mex,3)


# Create an empty matrix to store the two alleles of a SNP into one column
gih.geno <- data.frame(matrix(nrow=nrow(gih), ncol=ncol(gih)/2))
mex.geno <- data.frame(matrix(nrow=nrow(mex), ncol=ncol(mex)/2))

# Replace N with empty string
combine.snp <- function(snp1, snp2){
  snp1 <- as.character(snp1)
  snp2 <- as.character(snp2)
  snp1[snp1=="N"] <- ""
  snp2[snp2=="N"] <- ""
  snp.geno <- genotype(snp1,snp2)
  return(snp.geno)
}

k <- 1
for(i in 1:ncol(gih.geno)){
  gih.geno[,i] <- combine.snp(gih[,k], gih[,k+1])
  k <- k+2
}

k <- 1
for(i in 1:ncol(mex.geno)){
  mex.geno[,i] <- combine.snp(mex[,k], mex[,k+1])
  k <- k+2
}

## name the columns the proper rsIDs of the 20 SNPs
snp.info <- read.xls("HapMap_chrom22_20SNP_info.xls")
head(snp.info,3)

colnames(gih.geno) <- as.character(snp.info$SNP.name)
colnames(mex.geno) <- as.character(snp.info$SNP.name)

head(gih.geno,3)
head(mex.geno,3)

##########################################################################################################

#####################################################
################ Questions #########################

#### Question 1 ####
# ?HWE.chisq

## GIH ##
gih.HWE.pval = NULL
for (i in 1:ncol(gih.geno)) gih.HWE.pval[i]=HWE.chisq(gih.geno[,i], B= 100000)$p.value

# SNPs deviating from HWE for GIH
cat("SNP(s) that deviate(s) from HWE in the GIH population:", colnames(gih.geno)[gih.HWE.pval < 0.05], sep="\n")

## MEX ##
mex.HWE.pval = NULL
for (i in 1:ncol(mex.geno)) mex.HWE.pval[i]=HWE.chisq(mex.geno[,i], B=100000)$p.value

# SNPs deviating from HWE for MEX
cat("SNP(s) that deviate(s) from HWE in the GIH population:", colnames(mex.geno)[mex.HWE.pval < 0.05], sep="\n")

#### Question 2 ####
?LD

## GIH ##
gih.LD = LD(gih.geno)
# names(gih.LD)
# head(gih.LD$`R^2`, 3)
# head(gih.LD$`P-value`, 3)

gih.LD.SNPs <- apply(gih.LD$`P-value`, 1,  function(x) names(which(x < 0.05)) )

# LD pairs for GIH
gih.LD.pairs = NULL
for(i in 1:length(gih.LD.SNPs)) gih.LD.pairs[i] = toString(gih.LD.SNPs[[i]])

## MEX ##
mex.LD = LD(mex.geno)
mex.LD.SNPs <- apply(mex.LD$`P-value`, 1,  function(x) names(which(x < 0.05)))

# LD pairs for MEX
mex.LD.pairs = NULL
for(i in 1:length(mex.LD.SNPs)) mex.LD.pairs[i] = toString(mex.LD.SNPs[[i]])


# any common SNPs that show LD?

# number of common SNPs
gih.LD.SNPs2 <- names(gih.LD.SNPs[lapply(gih.LD.SNPs, length) > 0])
mex.LD.SNPs2 <- names(mex.LD.SNPs[lapply(mex.LD.SNPs, length) > 0])
length(mex.LD.SNPs2[mex.LD.SNPs2 %in% gih.LD.SNPs2])

common.LD = list()
for(i in 1:length(gih.LD.SNPs)) common.LD[[names(gih.LD.SNPs)[i]]] = gih.LD.SNPs[[i]][gih.LD.SNPs[[i]] %in% mex.LD.SNPs[[i]]]

common.LD <- common.LD[lapply(common.LD, length) > 0]
common.LD.pairs = NULL
for(i in 1:length(common.LD)) common.LD.pairs[i] = toString(common.LD[[i]])
common.LD.df <- data.frame(Common_SNPS_that_show_LD= names(common.LD) , common.LD.pairs =common.LD.pairs , row.names=NULL)

# Table 1
xtable(common.LD.df, caption="SNPs with common LD pairs")

# Table 2
LD.df = data.frame(GIH = gih.LD.pairs, MEX = mex.LD.pairs, row.names = colnames(gih.geno) )
LD.df.tbl = xtable(LD.df, caption = "LD pairs for 20 SNPs in GIH and MEX populations")
print.xtable(LD.df.tbl, size ="small")


### Question 3 ####
data(CEUData)
CEUDist
head(snp.info)
head(CEUSNP)
head(gih.geno)

## r2
gih.r2 <- LDheatmap(gih.geno, snp.info$Position)
gih.D <- LDheatmap(gih.geno, snp.info$Position, LDmeasure = "D'")

rgb.palette <- colorRampPalette(rev(c("blue", "yellow", "red")), space = "rgb") #create color plot
LDheatmap(gih.r2, color=rgb.palette(100), SNP.name = snp.info$SNP.name, text=T)
LDheatmap(gih.D, color=rgb.palette(100), SNP.name = snp.info$SNP.name, LDmeasure = "D'")



LDheatmap(mex.geno, snp.info$Position, SNP.name = snp.info$SNP.name)
LDheatmap(mex.geno, snp.info$Position, SNP.name = snp.info$SNP.name, LDmeasure = "D'")



MyHeatmap<-LDheatmap(gih.geno, snp.info$Position) # plots the pairwise LD measures; default is r2; 
rgb.palette <- colorRampPalette(rev(c("blue", "yellow", "red")), space = "rgb") #create color plot
LDheatmap(MyHeatmap, color=rgb.palette(100))
names(MyHeatmap) # shows the objects saved in MyHeatmap
MyHeatmap$ LDmatrix # shows specifically the object LDmatrix, which has the pairwise LD measures.
#Following are some useful options for LDheatmap
LDheatmap(MyHeatmap, SNP.name=c("rs617160", "rs617337")) # labels the two listed SNPs; note the saved object MyHeatmap can be used, i.e, there is no need to create the map again.
LDheatmap.highlight(MyHeatmap, i=3, j=8, col="red", fill="grey") # highlight the region covered between the ith and jth SNP
LDheatmap(gih.geno, snp.info$Position, SNP.name = snp.info$SNP.name) # labels all the SNPs LDheatmap(CEUSNP, CEUDist, SNP.name=names(CEUSNP)[3:6])
# labels only the SNPs 3 to 6. LDheatmap(CEUSNP, CEUDist, LDmeasure="D'") # plots D’ instead of the default r2
