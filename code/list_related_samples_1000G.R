ibd <- read.table("data/plink_pca_1000g/1000GxGDA.genome.gz", header = T, as.is = T)
head(ibd)

hist(ibd$PI_HAT, breaks = 100) # PI_HAT = Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
table(ibd$PI_HAT > 0.2)

exclusions = ibd[ibd$PI_HAT > 0.2, c('FID2','IID2')]
write.table( exclusions, file="data/plink_pca_1000g/related_samples.txt", col.names = F, row.names = F, quote = F)
