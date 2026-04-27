library(readxl)
library(tidyverse)
library(ggrepel)
library(WRS2)

# Read the tables with coordinates and MAF as well as the values from ADMIXTURE
slc <- read_xlsx("/home/jorgeg/PROJECTS.JORGE/CORRELATIONS_SNPS/Gene_SLC30A9_References.xlsx", sheet = 1) %>% select(c(Population, MAF, MAC,NCHROBS))
ADMIXTURE_4k <- read_xlsx("/home/jorgeg/PROJECTS.JORGE/CORRELATIONS_SNPS/41586_2023_6770_MOESM4_ESM.xlsx", sheet = 3, skip = 1) 

slc <- slc %>% group_by(Population) %>% summarise(MAF = sum(MAC)/sum(NCHROBS))

# Process
EUR_anc <- ADMIXTURE_4k %>% select(Population,`Europe_British-GBR`) %>% separate(`Europe_British-GBR`, into = "EUR", "%")

all_anc <- merge(slc,EUR_anc)
all_anc$EUR <- as.numeric(all_anc$EUR)

EAS_anc <- ADMIXTURE_4k %>% select(Population, `EastAsia_Japanese-JPT`) %>% separate(`EastAsia_Japanese-JPT`, into = "EAS", "%")
df_eas <- merge(slc,EAS_anc)
df_eas$EAS <- as.numeric(df_eas$EAS)

# Treat EUR and EAS as a single Eurasian ancestry
df <- merge(all_anc,EAS_anc) 
df$Eurasia <- df$EUR+round(as.numeric(df$EAS),3)

# Perform robust correlation
pbcor <-wincor(as.numeric(df$MAF),as.numeric(df$Eurasia))

ggplot(df, aes(x=MAF, y=Eurasia)) +
  geom_point() + # Show dots
  geom_text_repel(
    aes(label=Population), 
    check_overlap = T, size = 3, max.overlaps = 15, 
  ) + geom_smooth(formula = y ~ x, method = "lm") + labs(x = "rs1047626 MAF", y = "% of Eurasian Ancestry") + theme_classic() 

ggsave("/home/jorgeg/PROJECTS.JORGE/CORRELATIONS_SNPS/correlation_allele_Eurasiaancestry.png", device = "png", units = "px", height = 3000, width = 3200)
