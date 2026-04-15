# HaplotypAssigner

library(ggplot2);
library(tidyr);
library(dplyr);
library(ggforce)

#rm(list=ls());
rm(list = setdiff(ls(), "data.t"));

# estimate for each position the maximum length of the haplotype with the same allele than in the reference
# haplotype. In practice, any SNV with diff 0 to the reference haplotype
weight.positions.by.similarity.in.region <- function(h, coor)
{
  # vector to store the output for each SNV. Length of the fragment that is identical to
  # the reference haplotype. Returns 1 if the allele is different than the reference.
  v <- rep(1,length(h));
  # identify the positions that are equal to the reference. These are the ones we have to evaluate
  # the ones that are different imply that the length of the shared region is by default 1 bp
  equal.to.ref <- which(h==0);
  # identify the positions that are different to the reference. These are the positions that will define
  # the borders of the shared haplotypes
  diff.to.ref <- which(h!=0);
  # Iterate over all the positions equals to reference
  for(pos in equal.to.ref)
  {
    # positions that are different from the reference and smaller than current position that is equal to reference
    left <- diff.to.ref[coor[diff.to.ref] < coor[pos]];
    closest.zero.left <- 0;
    if(length(left) > 0)
    {
      # the closest point that is not anymore equal to the reference from the left
      closest.zero.left <- min(abs(coor[left]-coor[pos]));      
    }
    else
    {
      closest.zero.left <- NA;
    }
    right <- diff.to.ref[coor[diff.to.ref] > coor[pos]];
    closest.zero.right <- 0;
    if(length(right)>0)
    {
      # the closest point that is not anymore equal to the reference from the right
      closest.zero.right <- min(abs(coor[right]-coor[pos]));      
    }
    else
    {
      closest.zero.right <- NA;
    }
    # total length of haplotype that is
    v[pos] <- closest.zero.right + closest.zero.left;
  }
  # return the vector of fragment length
  return(v);
}

# Function to assign to proxy populations
assign.to.proxy <- function(haplotype, parental)
{
  diff.parental <- abs(sweep(parental, 2, haplotype, FUN = "-"));
  # If at least one of the two alleles is different, then set it to different
  diff.parental[diff.parental==0.5] <- 1;
  coor <- as.numeric(sub("^X", "", colnames(parental)));
  # matrix. Rows are haplotypes, columns are positions
  k <- (apply(diff.parental,1,weight.positions.by.similarity.in.region,coor));
  best.m <- cbind(rep(0,nrow(k)), k[,which(colnames(k)!="san")]);
  sans <- k[,which(colnames(k)=="san")];
  #print(tail(k, n = 30))
  best.m[,1] <- apply(sans,1,max, na.rm = T); # Change na.rm to TRUE
  # This generated lots of NAs when there were NAs in the best.m matrix
  # It was problematic with Mandenka population
  colnames(best.m)[1] <- "san";
  #print(rowSums(best.m))
  m_bin <- best.m == apply(best.m, 1, max);
  m_bin <- 1 * m_bin;
  # must add to 1
  m_bin <- m_bin / rowSums(m_bin);  
  return(m_bin);
}

# Read the data
data.t <- read.table(file="~/PROJECTS.JORGE/DENISOVANS_HAPLOTYPE_ASSIGN/mergedMatrix.gnomad_HGDP_1000G_grch37_phased.txt",header=T, na.strings = "NA");
pops <- data.t[,1];

# position in grch37 rs1047626 [Homo sapiens]
#Chromosome: 4:42001654 (GRCh38)
position_GRCh37 = 42003671;
# only haplotypes
X <- as.matrix(data.t[,3:ncol(data.t)]);
rownames(X) <- pops;
colnames(X) <- names(data.t)[3:ncol(data.t)];
# position of the SNP of interest
snp.of.interest <- which(colnames(X)==paste("X",position_GRCh37,sep=""));
# Pops of interest
MEL <- "melanesian";
EUR <- "ceu";
AFR <- "san";
ASI <- "chb";
OCE <- "papuan";
AFRA <- "mandenka";
#ARC <- c("DenisovaPinky", "Vindija33.19"); # A check. Same results.
ARC <- c("Chagyrskaya-Phalanx", "AltaiNea", "DenisovaPinky", "Vindija33.19"); #"Mez1"
ANC <- c("Ust_Ishim","Loschbour");
# make genotypes for anc
X.anc <- X[which(rownames(X) %in% ANC),];
X.anc <- (X.anc[seq(from=1,to=nrow(X.anc),by=2),] + X.anc[seq(from=2,to=nrow(X.anc),by=2),])/2;
X.anc[X.anc==0.5] <- NA;
# retrieve Oceania
X.oceania <- X[which(rownames(X) %in% c(EUR,OCE,MEL, ASI)),];
# homozygote 2 for the derived allele
X.oceania <- X.oceania[X.oceania[,snp.of.interest]==1,]; # CHANGE SNPS WIHT OR WITHOUT
X.afra.0 <- X[rownames(X) == AFRA,];
X.afra.0 <- X.afra.0[X.afra.0[,snp.of.interest]==0,];
X.oceania <- rbind(X.oceania,X.afra.0);
#X.oceania <- rbind(X.oceania,X.anc);
# homozygote 0 for the ancestral allele (in fact, all are homozygote for that allele)
X.afr <- X[rownames(X)==AFR,];
X.afr <- X.afr[X.afr[,snp.of.interest]==0,];
# make genotypes for arc
X.arc <- X[which(rownames(X) %in% ARC),];
X.arc <- (X.arc[seq(from=1,to=nrow(X.arc),by=2),] + X.arc[seq(from=2,to=nrow(X.arc),by=2),])/2;
# variables with missing in all the columns of X.arc and X.oceania
not.missing.all <- colMeans(is.na(rbind(X.oceania,X.arc))) == 0;
# we will have to remove these markers from the full database
X.afr <- X.afr[,not.missing.all];
X.oceania <- X.oceania[,not.missing.all];
X.arc <- X.arc[,not.missing.all];
## merge the archaic with the african San
X.arc.afr <- rbind(X.afr, X.arc);
# Identify markers that are not variable in X.arc.afr
new.position <- which(colnames(X.oceania)==paste("X",position_GRCh37,sep=""));
# number of SNPs surrounding the position of interest
size <- 3500; # 3500

### Retrieve data
do_for_indiv <- function(i,haplotype, j, population) {
  length(haplotype);
  parental <- X.arc.afr[,(new.position-size):(new.position+size)];
  coor <- as.numeric(sub("^X", "", colnames(parental)));
  ppos <- which(coor==42003671);
  labs <- rownames(parental);
  # Only consider markers that are in, at least TWO haplotypes across the haplotype and the parental
  nw <- 1 + nrow(parental);
  #nw <- nrow(parental)
  freq <- colMeans(rbind(haplotype,parental));
  maf <- apply(cbind(freq,1-freq),1,min);
  haplotype <- haplotype[maf > 1/nw];
  parental <- parental[,maf>1/nw];
  coor <- as.numeric(sub("^X", "", colnames(parental)));
  ppos <- which(coor==42003671);
  
  m <- assign.to.proxy(haplotype, parental)
  rownames(m) <- coor;
  #m <- m[rowSums(is.na(m))==0,];
  coor <- as.numeric(sub("^X", "", rownames(m)));
  
  #m_bin <- m == apply(m, 1, max);
  #m_bin <- 1 * m_bin;
  
  m_prop <- m;
  
  df_long <- as.data.frame(m_prop)
  df_long$idx  <- seq_len(nrow(m_prop))
  df_long$coor <- coor
  
  df_long <- pivot_longer(
    df_long,
    cols = -c(idx, coor),
    names_to = "category",
    values_to = "value"
  )
  
  df_long$individual <- j
  
  df_long$pop <- population
  
  return(df_long)
}

# A function to apply the algorithm to each haplotype and store it
get_pop_df <- function(pop, pop_name) {
  list_of_dfs <- list()
  for (j in 1:sum(rownames(X.oceania)==pop) ) {
    i <- which(rownames(X.oceania)==pop)[j];
    haplotype <- X.oceania[i,(new.position-size):(new.position+size)];
    df_j <- do_for_indiv(i,haplotype,j,population = pop_name)
    list_of_dfs[[j]] <- df_j
  }
  
  df_pop <- dplyr::bind_rows(list_of_dfs)
  return(df_pop)
}

df_mel <- get_pop_df(MEL,"Melanesia")
df_CEU <- get_pop_df(EUR, "Central European in Utah")
df_CHB <- get_pop_df(ASI, "Han Chinese")
df_AFRA <- get_pop_df(AFRA, "Mandenka")
df_OCE <- get_pop_df(OCE, "Papuan")

#### PLOTTING FUNCTION
n_ind <- 70*3

to_plot <- function(df_pop,n_ind, filename,typp="PDF") {
  v_idx <- which.min(abs(unique(df_pop$coor) - 42003671))
  p <- ggplot(df_pop, aes(x = idx, y = value, fill = category)) +
    geom_col(width = 1) +
    facet_grid_paginate(individual ~ pop, 
                        scales = "free_x", 
                        space = "free_x",
                        ncol = n_ind,
                        nrow = length(unique(df_pop$pop)),
                        page = 1) +
    scale_fill_discrete(type = c("#008837","#a6dba0","#7b3294","darkgrey", "#c2a5cf")) +
    scale_x_continuous(
      breaks = unique(df_pop$idx)[seq(1, length(unique(df_pop$idx)), length.out = 5)],
      labels = unique(sort(df_pop$coor))[seq(1, length(unique(df_pop$idx)), length.out = 5)]
    ) +
    geom_vline(xintercept = v_idx, linetype = "dashed", linewidth = 0.5) +
    theme_minimal() +
    theme(
      axis.text.y  = element_blank() ,
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text.y = element_blank(),
      strip.background.y = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.size = unit(3, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.spacing.x = unit(1, "mm")
    ) 
  
  n_pages <- ceiling(length(unique(df_pop$individual))/n_ind)

  if (typp == "PDF") {
    
  
  pdf(filename, width = 8.27, height = 11.69, onefile = TRUE) #width = 14, height = 6)
  for(i in 1:n_pages){
    print(
      p + ggforce::facet_grid_paginate(
        individual ~ pop,
        scales = "free_y",
        space = "free_y",
        ncol = 1,
        nrow = n_ind,
        page = i
      ) )
    }
  }
  
      else if (typp == "PNG")
      {
      png(filename, width = 8.27, height = 11.69, units = "in", res = 120)
      print(p + ggforce::facet_grid_paginate(individual ~ pop, individual ~ pop,
                     scales = "free_y",
                     space = "free_y",
                     ncol = 1,nrow = n_ind) )
      
      }
  dev.off()
  
}

to_plot(df_pop = df_mel, n_ind, filename = "Mel_2.pdf", typp="PDF")
to_plot(df_pop = df_CHB, n_ind, filename = "CHB_2.pdf")
to_plot(df_pop = df_OCE, n_ind, filename = "Papuans_2.pdf")
to_plot(df_pop = df_AFRA, n_ind, filename = "Mandenka_2.pdf")
to_plot(df_pop = df_CEU, 300, filename = "CEU_2.pdf")

to_plot_png <- function(df_pop, filename) {
  
  v_idx <- which.min(abs(unique(df_mel$coor) - 42003671))
  p <- ggplot(df_pop, aes(x = idx, y = value, fill = category)) +
    geom_col(width = 1) +
    facet_grid(individual ~ pop, 
               scales = "free_y", 
               space = "free_y") +
    scale_fill_discrete(type = c("#008837","#a6dba0","#7b3294","darkgrey", "#386cb0")) +
    scale_x_continuous(
      breaks = unique(df_pop$idx)[seq(1, length(unique(df_pop$idx)), length.out = 5)],
      labels = unique(sort(df_pop$coor))[seq(1, length(unique(df_pop$idx)), length.out = 5)]
    ) +
    geom_vline(xintercept = v_idx, linetype = "dashed", linewidth = 0.5) +
    theme_minimal() +
    theme(
      axis.text.y  = element_blank() ,
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text.y = element_blank(),
      strip.background.y = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.size = unit(3, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.spacing.x = unit(1, "mm")
    ) 
  
  png(filename,width = 8.27, height = 11.69, units = "in", res = 120)
  print(p) 
  dev.off()
}




to_plot_png(df_CEU, filename = "CEU_PNG_2.png")
to_plot_png(df_pop = df_mel, filename = "Mel_PNG.png")
to_plot_png(df_pop = df_CHB, filename = "CHB_PNG.png")
to_plot_png(df_pop = df_OCE, filename = "Papuans_PNG.png")
to_plot_png(df_pop = df_AFRA, filename = "Mandenka_PNG.png")


#####
# Control 
#####

# just repeating the same but for the control haplotypes
# for which rs1047626 is ancestral
# We are only repeating the processing of the matrix
X.oceania <- X[which(rownames(X) %in% c(EUR,OCE,MEL, ASI)),];
X.oceania <- X.oceania[X.oceania[,snp.of.interest]==0,]; # rs1047626 ancestral
X.afra.0 <- X[rownames(X) == AFRA,];
X.afra.0 <- X.afra.0[X.afra.0[,snp.of.interest]==0,];
X.oceania <- rbind(X.oceania,X.afra.0);
X.afr <- X[rownames(X)==AFR,];
X.afr <- X.afr[X.afr[,snp.of.interest]==0,];
X.arc <- X[which(rownames(X) %in% ARC),];
X.arc <- (X.arc[seq(from=1,to=nrow(X.arc),by=2),] + X.arc[seq(from=2,to=nrow(X.arc),by=2),])/2;
# variables with missing in all the columns of X.arc and X.oceania
not.missing.all <- colMeans(is.na(rbind(X.oceania,X.arc))) == 0;
# we will have to remove these markers from the full database
X.afr <- X.afr[,not.missing.all];
X.oceania <- X.oceania[,not.missing.all];
X.arc <- X.arc[,not.missing.all];
## merge the archaic with the african San
X.arc.afr <- rbind(X.afr, X.arc);
# Identify markers that are not variable in X.arc.afr
new.position <- which(colnames(X.oceania)==paste("X",position_GRCh37,sep=""));

df_CEU_CONTROL <- get_pop_df(EUR, "Central European in Utah - CONTROL")
to_plot(df_pop = df_CEU_CONTROL, n_ind, filename = "CEU_CONTROL.pdf")
to_plot_png(df_pop = df_CEU_CONTROL, filename = "CEU_CONTROL_PNG.png")
# End
