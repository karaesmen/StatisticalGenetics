setwd("~/Google Drive/AE/ProjectF/ezgi/")
library("SNPlocs.Hsapiens.dbSNP.20090506")
available.SNPs()

narac <- read.table("narac_chr6.txt", header = T, stringsAsFactors = F, na.strings = "NA")
snpid.narac <- colnames(narac)[3:ncol(narac)]

# make sure right package is used 2009 = hg18
snp.loc <- SNPlocs.Hsapiens.dbSNP.20090506::rsidsToGRanges(snpid.narac)

# Error in .snpid2rowidx(x, snpid) : SNP id(s) not found: 9258122
# exclude 9258122
snpid.narac = snpid.narac[-grep("9258122", snpid.narac)]

snp.loc <- SNPlocs.Hsapiens.dbSNP.20090506::rsidsToGRanges(snpid.narac)
head(snp.loc)
snp.loc <- data.frame(snp.loc)
head(snp.loc)
snp.loc <- cbind(rsID=paste0("rs", snp.loc$RefSNP_id), snp.loc[,c(1, 2)] )
colnames(snp.loc) <- c("rsID", "chr", "position" )

hapmap3.chr6 <- read.table("../hapmap3/genetic_map_chr6_combined_b36.txt", header=T)
head(hapmap3.chr6)
narac.hapmap3.map <- hapmap3.chr6[hapmap3.chr6$position %in% as.integer(snp.loc$start),]

narac.hapmap3.map <- merge(hapmap3.chr6, snp.loc, by.y = "position")
head(narac.hapmap3.map)

# only missing 3 SNPs
nrow(narac.hapmap3.map)
range(narac.hapmap3.map$position)

# create map for IMPUTE2
write.table(narac.hapmap3.map[,1:3], "narac_hapmap3_IMP2.map", sep="\t", quote = F, row.names = F)

### create .map & .ped files to generate .GEN file

# .map file
narac_hapmap_plink.map <- cbind(rep(6, nrow(narac.hapmap3.map)), 
                                narac.hapmap3.map[,c("rsID", "Genetic_Map.cM.", "position")] )

write.table(narac_hapmap_plink.map, "narac_hapmap3.map", sep="\t", 
            quote = F, row.names = F, col.names = F)

# .ped file

# get only matched SNPs from the original file
indx <- which(colnames(narac) %in% narac_hapmap_plink.map$rsID)
narac.ped <- narac[, c(1:2, indx)]
ncol(narac) ; ncol(narac.ped)

narac.cov <- read.table("cov.txt", header = F)
colnames(narac.cov) <- c("FID", "ID", "COV")

narac.ped <- merge(narac.cov[,1:2], narac.ped, by.y="ID")
ncol(narac.ped)
narac.ped <- cbind(narac.ped$FID, narac.ped$ID, 
                     rep(0, nrow(narac.ped)), 
                     rep(0, nrow(narac.ped)), 
                     rep(0, nrow(narac.ped)), 
                     narac.ped[,3:ncol(narac.ped)])


#install.packages('splitstackshape')
narac.ped <- data.frame(apply(narac.ped, 2, FUN=function(x){ x[is.na(x)] <- "0_0"; x}))

library(splitstackshape)
narac.ped <- data.frame(cSplit(narac.ped, colnames(narac.ped)[7:ncol(narac.ped)], sep="_"))

write.table(narac.ped, file="narac_hapmap3.ped", sep="\t", quote=F, col.names=F, row.names=F)

narac.cov[56,]

## read in the impute2_info


narac.impute.info <- read.table("narac.chr6.impute2_info", header = T, stringsAsFactors = F)
head(narac.impute.info)

narac.impute.info[order(narac.impute.info$position),]

narac.impute.info.excl <- subset(narac.impute.info, exp_freq_a1 < 0.01 & info < 0.3)



head(narac.impute.info.filt)

dim(narac.impute.info) ; dim(narac.impute.info.filt)

narac.impute.info.filt[order(narac.hapmap3.map)]

# #  read in imputed.map 
# 
# imp.map <- read.table("narac_imputed.map", header = F)
# imp.ped <- read.table("narac_imputed.ped", header = F)
#  
# dim(imp.ped) ; dim(imp.map)
# 

## 



#write.table(narac.impute.info.filt$rs_id, "include_rsids.txt", col.names=F, row.names=F, quote=F)


# 



#BiocInstaller::biocLite("VariantAnnotation")
library(VariantAnnotation)
narac.vcf <- readVcf("narac.chr6.vcf", genome="hg18")

narac.gt <- geno(narac.vcf)$GT
narac.gt <- t(narac.gt)
narac.gt[1:5,1:10]
