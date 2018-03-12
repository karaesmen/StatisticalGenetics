setwd("~/Box Sync/Stat Anal of Genetic Data/Project2/Description_Data/")

map.cM <- scan("chrom12_marker_info.txt")
map.cM <- c(0, map.cM)
map.M  <- map.cM/100
  
map.file <- data.frame(CHR = rep(12, length(map.cM)),
                       SNP = paste0("SNP", 1:length(map.cM)),
                       M   = map.M,
                       BP  = map.cM*10^6)
head(map.file)
write.table(map.file, "chrom12.map", sep="\t",col.names = F, row.names = F, quote = F)


### NARAC ####

## narac.ped
narac <- read.table("narac_chr6.txt", header = T, stringsAsFactors = F)
narac[is.na(narac)] <- "0_0"
ncol(narac)
narac[1:3, 1:5]

#anyNA(narac)


rsIDs = colnames(narac[ ,-c(1,2)])
temp2 = list()
for(i in rsIDs){
  temp = NULL
  temp <- narac[,i]
  temp = strsplit(temp, split="_")
  temp2[[i]] = do.call(rbind, temp)
}

genotypes <- do.call(cbind, temp2)
narac.ped = data.frame(FID= c(1:nrow(narac)),narac[,1:2], genotypes)
narac.ped[1:10, 1:10]
write.table(narac.ped, "narac.ped", col.names = F, row.names = F, quote = F)
  
## narac.map

narac.map <- data.frame(CHR = rep(6, length(rsIDs)),
           SNP = rsIDs,
           M   = rep(0, length(rsIDs)),
           BP  = 1:length(rsIDs)+10^6)
write.table(narac.map, "narac.map", col.names = F, row.names = F, quote = F)


### NARAC-Q3

tmax <- read.table("narac.model.best.mperm", header = T)
assoc <- read.table("narac.model", header = T)

dim(assoc)
dim(tmax)

nrow(subset(tmax, EMP2 < 0.05) )
nrow(subset(assoc, P < 0.05) )


## read in logis.assoc

log.assc <- read.table("./narac.assoc.logistic", header = T)
head(log.assc)

adj.assc <- read.table("./narac.assoc.logistic.adjusted", header = T)
head(adj.assc)

subset(adj.assc, BONF < 0.05)
snps<- adj.assc[adj.assc$BONF < 0.05, "SNP"]

subset(log.assc, TEST == "ADD" & SNP %in% snps)

nrow(subset(log.assc, TEST == "COV1" & P < 0.05) )

cov <- subset(log.assc, TEST == "COV1" & P < 0.05) 
summary(cov$P)


--no-fid --no-parents --no-sex
