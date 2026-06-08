## Convert genlight to BayPass file ##
#Code adapted from Ahrens et al. (2021).https://doi.org/10.1111/1755-0998.13351

#load packages
library(gdata)

#Load genlight object
load("~/Desktop/Comparative-Landscape-Genomics/1_SNP_quality_control_filtering/1_1_BM/BM_gle.rdata")

#Convert genlight object to data frame
BM_df <- as.data.frame(BM_gle)
write.table(BM_df, "BM_012.txt", quote = FALSE)

#Import LONGLAT file, this should have all individuals in the same order as INPUT file
ind <- gle_order@ind.names
write.csv(ind, file = "BM_ind.csv")

pop <- gle_order@pop
write.csv(pop, file = "BM_pop.csv")

loc.names <- gle_order@loc.names
write.csv(loc.names, file = "BM_loc.names.csv")

longlat <- read.csv("/Volumes/Backup\ Plus/0_OneDrive/Thesis/3_Experiment/2_Chapter/2_2_Data/Genomic/BM/maf_2/run_2_2/gl2baypass/BM_lonlat.csv", header = TRUE)

#Create function for converting to BayPass file 
baypass_format <- function(INPUTFILE,LONGLAT,BAYPASSFILE){
  gen<-read.table(INPUTFILE, header = T, row.names=1)
  lonlat = LONGLAT
  allele1 = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (s)
      })
  })
  allelecount = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (n)
      })
  })
  allelefreq = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (s / n)
      })
  })
  allele2<-allelecount-allele1
  whole<-interleave(allele1,allele2)
  baypass<-t(whole)
  write.table(baypass,BAYPASSFILE, col.names = F, row.names = F, sep="\t")
}

# Run baypass conversion function
baypass_format(INPUTFILE = "BM_012.txt", LONGLAT = longlat, BAYPASSFILE = "BM_Baypass")
