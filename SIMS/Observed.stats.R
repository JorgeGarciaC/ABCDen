library(tidyverse)
source("~/scratch/rstudio-singularity/All_statistics.R")

args = commandArgs(trailingOnly=TRUE)

replicate <- 1
l <- as.integer(args[1]) # 41851654-42151654
chunk <- args[2] 
dir <- "/homes/users/jgarciac/scratch/ABC.COMPUTING"
haps <- read.table(paste0(dir, "/merged.chr4.", chunk, ".norm.subset.flip.haps"))[,c(-1,-2,-4,-5)]
sample <- read.table(paste0(dir, "/merged.chr4.", chunk, ".norm.subset.sample"))[c(-1,-2),1]
POS <- haps[,1]
gt <- haps[,-1]
colnames(gt) <- paste0(rep(sample,each=2), c("_1","_2"))
haps_mod <- read.table(paste0(dir, "/modern.chr4.", chunk, ".subset.flip.haps"))[,c(-1,-2,-3,-4,-5)]
colnames(haps_mod) <- paste0(rep(sample[c(-297,-296,-295,-294)],each=2), c("_1","_2"))
names_poplabs <- c("Sample","Pop") 
pat_pops <- paste0(dir, "/pop.labels.txt")
poplabels <- read_table(pat_pops, col_names = names_poplabs )
aux <- data.frame(Sample = paste0(poplabels$Sample, "_1"), Pop = poplabels$Pop)
aux2 <- data.frame(Sample = paste0(poplabels$Sample, "_2"), Pop = poplabels$Pop)
poplabels <- rbind(aux,aux2)

# Extract haplotypes for each pop
subset_haplotype <- function(gt,POP) {
  df <- gt %>% dplyr::select(all_of(unlist(poplabels[poplabels$Pop == POP,1]))) %>% 
    mutate(across(everything(), ~ as.numeric(.))) %>% as.matrix()
  return(df)
}

CHB <- subset_haplotype(gt, "CHB")
EUR <- subset_haplotype(gt,"GBR")
EAF <- subset_haplotype(gt, "LWK")
modCHB <- subset_haplotype(haps_mod, "CHB")
modEUR <- subset_haplotype(haps_mod,"GBR")
modEAF <- subset_haplotype(haps_mod, "LWK")
NEA <- gt %>% dplyr::select(all_of(c("AltaiNeandertal_1", "Vindija33.19_1", "Chagyrskaya-Phalanx_1","AltaiNeandertal_2", "Vindija33.19_2", "Chagyrskaya-Phalanx_2"))) %>% 
  mutate(across(everything(), ~ as.numeric(.))) %>% as.matrix()
DEN <- gt %>% dplyr::select(starts_with("Denisova")) %>% 
  mutate(across(everything(), ~ as.numeric(.))) %>% as.matrix()

popsize_EAF <- (dim(EAF)[2])/2
popsize_CHB <- (dim(CHB)[2])/2
popsize_EUR <- (dim(EUR)[2])/2
popsize_DEN <- (dim(DEN)[2])/2
popsize_NEA <- (dim(NEA)[2])/2

Fr_EAF <- rowSums(EAF,na.rm = T)/(dim(EAF)[2]) # Unadmixed P1
Fr_CHB <- rowSums(CHB, na.rm = T)/(dim(CHB)[2]) # admixed P2
Fr_EUR <- rowSums(EUR, na.rm = T)/(dim(EUR)[2]) # admixed P2
Fr_DEN <- rowSums(DEN, na.rm = T)/(dim(DEN)[2]) # Archaic introgressed P3
Fr_NEA <- rowSums(NEA, na.rm = T)/(dim(NEA)[2])
P0 <- 0 

# This is important, due to ancient DNA miss positions, simulations data has more 
# segregatting sites than the merged VCFs, but similar to the "raw" VCF from modern day 
Fr_EAFmod <- rowSums(modEAF,na.rm = T)/(dim(modEAF)[2]) # Unadmixed P1
Fr_CHBmod <- rowSums(modCHB, na.rm = T)/(dim(modCHB)[2]) # admixed P2
Fr_EURmod <- rowSums(modEUR, na.rm = T)/(dim(modEUR)[2]) # admixed P2

#### Calculate stats

DfD_CHB_DEN <- getD_stats(Fr_EAF,Fr_CHB,Fr_DEN)
DfD_EUR_DEN <- getD_stats(Fr_EAF,Fr_EUR,Fr_DEN)
U_CHB_DEN <- extract.U(Fr_EAF,Fr_CHB,Fr_DEN,0.1,0.2) # U.10.20.100
U_EUR_DEN <- extract.U(Fr_EAF,Fr_EUR,Fr_DEN,0.1,0.2) # U.10.20.100

DfD_CHB_NEA <- getD_stats(Fr_EAF,Fr_CHB,Fr_NEA)
DfD_EUR_NEA <- getD_stats(Fr_EAF,Fr_EUR,Fr_NEA)
U_CHB_NEA <- extract.U(Fr_EAF,Fr_CHB,Fr_NEA,0.1,0.2) # U.10.20.100
U_EUR_NEA <- extract.U(Fr_EAF,Fr_EUR,Fr_NEA,0.1,0.2) # U.10.20.100

# Segregating sites
S_EAF <- S_sites(Fr_EAFmod)
S_CHB <- S_sites(Fr_CHBmod)
S_EUR <- S_sites(Fr_EURmod)
S_DEN <- S_sites(Fr_DEN)
S_NEA <- S_sites(Fr_NEA)

PS_EAF <- sum(P_seg(Fr_EAFmod))/S_EAF
PS_CHB <- sum(P_seg(Fr_CHBmod))/S_CHB
PS_EUR <- sum(P_seg(Fr_EURmod))/S_EUR

dxy_CHB_DEN <- sum(Fr_DEN*(1-Fr_CHB)+(1-Fr_DEN)*Fr_CHB)/l 
dxy_EUR_DEN <- sum(Fr_DEN*(1-Fr_EUR)+(1-Fr_DEN)*Fr_EUR)/l 
dxy_EAF_DEN <- sum(Fr_DEN*(1-Fr_EAF)+(1-Fr_DEN)*Fr_EAF)/l
R_d_CHB_D <- dxy_CHB_DEN/dxy_EAF_DEN
R_d_EUR_D <- dxy_EUR_DEN/dxy_EAF_DEN


dxy_CHB_NEA <- sum(Fr_NEA*(1-Fr_CHB)+(1-Fr_NEA)*Fr_CHB)/l 
dxy_EUR_NEA <- sum(Fr_NEA*(1-Fr_EUR)+(1-Fr_NEA)*Fr_EUR)/l 
dxy_EAF_NEA <- sum(Fr_NEA*(1-Fr_EAF)+(1-Fr_NEA)*Fr_EAF)/l
R_d_CHB_N <- dxy_CHB_NEA/dxy_EAF_NEA
R_d_EUR_N <- dxy_EUR_NEA/dxy_EAF_NEA

tw_EAF <- thetaW(S_EAF,popsize_EAF)/l
tw_CHB <- thetaW(S_CHB,popsize_CHB)/l
tw_EUR <- thetaW(S_CHB,popsize_EUR)/l
tp_EAF <- thetapi(PS_EAF,popsize_EAF)/l
tp_CHB <- thetapi(PS_CHB,popsize_CHB)/l
tp_EUR <- thetapi(PS_EUR,popsize_EUR)/l

Hs_CHB <- get_H12(modCHB)
Hs_EUR <- get_H12(modEUR)
Hs_EAF <- get_H12(modEAF)

reticulate::use_python("/homes/users/jgarciac/.conda/envs/Python-3.12_msprime-1.3.1_tskit-0.5.6_pyslim-1.0.4_tspop-0.0.2/bin/python", required = TRUE)
#reticulate::import_from_path("/homes/users/jgarciac/scratch/jupyter-singularity/test/lib/python3.11/site-packages/allel/")
#reticulate::py_module_available("allel")
reticulate::source_python("/homes/users/jgarciac/scratch/rstudio-singularity/ZnS_observed.py", convert = T)
modVCF <- paste0(dir, "/modern.chr4.", chunk, ".subset.vcf.gz")
ancVCF <- paste0(dir, "/merged.chr4.", chunk, ".norm.subset.vcf.gz")
ZnS_CHB <- Calculate_ZnS_python(modVCF, pat_pops, "CHB")
ZnS_EUR <- Calculate_ZnS_python(modVCF, pat_pops, "GBR")
ZnS_EAF <- Calculate_ZnS_python(modVCF, pat_pops, "LWK")
ZnS_DEN <- Calculate_ZnS_python(ancVCF, pat_pops, "DEN")
ZnS_NEA <- Calculate_ZnS_python(ancVCF, pat_pops, "NEA")

# Because the frequence of the allele is only important on the CORE central
# part of the putatively introgressed segment
if (chunk == "CORE") {
	fr_allel_EUR <- Fr_EUR[POS == 42001654]
	fr_allel_EAF <- Fr_EAF[POS == 42001654]
	fr_allel_CHB <- Fr_CHB[POS == 42001654] 
} else { 
	# They will be remove in the following proccesing steps but for 
  # easing the next steps we just output the known allele frequencies
	fr_allel_EUR <- 0.7802198
	fr_allel_EAF <- 0.08080808
	fr_allel_CHB <- 0.9708738
}

vector.of.statistics <- c(DfD_CHB_DEN[1],DfD_CHB_DEN[2],DfD_EUR_DEN[1],DfD_EUR_DEN[2],
                          DfD_CHB_NEA[1],DfD_CHB_NEA[2],DfD_EUR_NEA[1],DfD_EUR_NEA[2],
                          R_d_CHB_D,R_d_EUR_D,R_d_CHB_N,R_d_EUR_N,
                          U_CHB_DEN,U_EUR_DEN,U_CHB_NEA,U_EUR_NEA,
                          S_CHB,S_EUR,S_EAF,S_DEN,S_NEA,
                          tw_CHB,tw_EUR,tw_EAF,
                          tp_CHB,tp_EUR,tp_EAF,
                          Hs_CHB[1],Hs_CHB[2],Hs_CHB[3],Hs_CHB[4],
                          Hs_EUR[1],Hs_EUR[2],Hs_EUR[3],Hs_EUR[4],
                          Hs_EAF[1],Hs_EAF[2],Hs_EAF[3],Hs_EAF[4],
                          fr_allel_CHB,fr_allel_EUR,fr_allel_EAF,
                          1,0,0,1,0
                          ,ZnS_CHB,ZnS_EUR,ZnS_EAF,ZnS_DEN,ZnS_NEA,chunk
)

write.table(t(vector.of.statistics), file = paste0(dir, "/obs.stats.chr4.", chunk, ".tsv"), quote = F, sep = "\t", append = F,row.names = F, col.names = F)
