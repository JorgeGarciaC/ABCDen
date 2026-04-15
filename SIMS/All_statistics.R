### ALL STATISTICS COMPUTED

options(dplyr.summarise.inform = FALSE)

getD_stats <- function(P1,P2,P3) {
  log_v <- (P1 == 0 & P2 == 0 & P3 ==0)
  P1 <- P1[!log_v]
  P2 <- P2[!log_v]
  P3 <- P3[!log_v]
  Pd <- pmax(P2,P3)
  P0 <- 0 
  
  ABBAs <- (1 - P1) * P2 * P3 * (1 - P0) # calculamos ABBAs
  BABAs <- P1 * (1 - P2) * P3 * (1 - P0) # Calculamos BABAs
  Max_ABBAs <- (1 - P1) * Pd * Pd * (1 - P0) 
  Max_BABAs <- P1 * (1 - Pd) * Pd * (1 - P0)
  
  if (sum(ABBAs + BABAs) == 0) {
    D <- 0
  } else { D <- sum(ABBAs - BABAs) / sum(ABBAs + BABAs) }
  
  if (sum(Max_ABBAs - Max_BABAs) == 0) {
    fD <- 0
  } else { fD <- sum(ABBAs - BABAs) / sum(Max_ABBAs - Max_BABAs)  }
  
  return(c(D,fD))
}

# U
extract.U <- function(bait,target,ref,w,x,y=1) {
  return(sum((bait < w) & (target > x) & (ref==1)) + # when ref == homozygote derived
           sum((1-bait < w) & (1-target > x) & ref==1-y)) # when ref == homozygote ancestral, we calculated for 1-P and the condition is the opposite
}

# Q
extract.Q <- function(P1,P2,P3,w,x=1,q) {
  # P are numeric vectors with the allele frequencies for each pop f# or each site
  return(quantile(P2[P1 < w & P3 ==x],q))
}

# N segergating sites

S_sites <- function (P) {
  S <- sum(P > 0 & P < 1)
  return(S)
}

P_seg <- function (P) {
  P_s <- P[P > 0 & P < 1]
  return(P_s)
}

# ThetaW Watterson Stimator
thetaW <- function (S_pop,popsize) {
  # *2 cause it's sequences number
  # The denominator is the harmonic number relative to n_seq
  thetaW <- S_pop/sum(1/1:((popsize*2)-1)) 
  return(thetaW)
}

# Thetapi
thetapi <- function(P,popsize){
  pi <- sum(2*P*(1-P)) # frequency of heterozygotes
  n <- popsize *2
  tpi <- pi*n/(n-1) # Heterozygosity
  return(tpi)
}

# theta H

#gene_diversity 
#options(dplyr.summarise.inform = FALSE)
#hcounts <- data.frame(t(target)) %>%
#  group_by_all() %>% # grouping by all variables we look at unique combinations
#  summarise(count = n(),) %>% 
#  ungroup() %>%
#  select(count)

#n <- popsize_target*2
#Hg (n/(n-1))*(1-(hcounts/sum(hcounts))^2)

#### thetaH
# An implementation of theta Nei following https://doi.org/10.1111/2041-210X.13643

get_thetaH <- function(haplotype_pop) {
  aux.df <- data.frame(t(haplotype_pop)) %>%
    group_by_all() %>% # grouping by all variables we look at unique combinations of hap
    summarise(count = n()) 
  
  dff <- function(rows,aux.df) {
    fr <- aux.df$count/sum(aux.df$count) # relative fr of each haplotype
    aux.df <- aux.df %>% dplyr::select(-count) 
    kij <- sum(!aux.df[rows[2],] == aux.df[rows[1],]) # number of sites different for each comb
    pi <- fr[rows[2]] 
    pj <- fr[rows[1]]
    return(c(kij*pi*pj))
  }
  
  if (length(aux.df$count) < 2) {
    return(0)
  } else {
    C <- combn(length(aux.df$count),2) # Obtain all pairwaise comp of hapl
    summa <- apply(C, 2, function(rows) {dff(rows,aux.df)}) # To each column fo the combinatorial object  we apply the function
    thetaH <- (sum(aux.df$count)/(sum(aux.df$count)-1)) * sum(summa) # eq 5  https://doi.org/10.1111/2041-210X.13643
    return(thetaH)
  }
  
}

# H1, H2, H12, H2/H1
get_H12 <- function(haplotype_pop) {
  # Take a matrix of haplotypes and gives a
  # vector with the statistics
  
  row_c <- data.frame(t(haplotype_pop)) %>%
    group_by_all() %>% # grouping by all variables we look at unique combinations
    summarise(count = n()) %>% 
    ungroup() %>%
    dplyr::select(count) %>% 
    arrange(desc(count))
  if (dim(row_c)[1] > 1) {
    fr_H <- row_c$count/sum(row_c$count) #fr of each unique haplotype
    H1 <- sum(fr_H**2) # great power to detect both 
    H2 <- H1-fr_H[1]**2 # Lower in hard sweeps
    H12 <- H1+2*fr_H[1]*fr_H[2] # Identify 
    H2H1 <- H2/H1 # Should increase for soft sweeps
    return(c(H1,H2,H12,H2H1))
  }
  else {
    # if there is only one haplotype, then H1 is 1, H2 is 0
    return(c(1,0,1,0))
  }
}

#### Linkage disequilibrium
# Based on: A Test of Neutrality Based on Interlocus Associations
# John K Kelly

ZnS_ld <- function(pop_hap){
  P_hap <- rowSums(pop_hap)/(dim(pop_hap)[2])
  aux <- pop_hap[P_hap > 0 & P_hap < 1,]
  S <- nrow(aux) # number of polymorphic sites
  
  # Number of possible combinations for the number of variants in the window
  C <- combn(nrow(aux),2) 
  
  # we use each row in the combination (C) as index for the matrix
  # we count the ocurrences for the whole population for each of the 4 states using table
  # We use apply in columns, cause in the combination object they are keep as columns,
  
  r2values <- apply(C, 2, function(pair) {
    hap1 <- aux[pair[1], ]
    hap2 <- aux[pair[2], ]
    
    table_states <- table(hap1,hap2)
    n <- sum(table_states)
    
    # The if states are verbose, but better being robust
    # Calculate the frequencies of the haplotype states
    pAB <- if ("1" %in% rownames(table_states) & "1" %in% colnames(table_states)) table_states["1","1"] / n else 0
    pAb <- if ("1" %in% rownames(table_states) & "0" %in% colnames(table_states)) table_states["1","0"] / n else 0
    paB <- if ("0" %in% rownames(table_states) &  "1" %in% colnames(table_states)) table_states["0","1"] / n else 0
    pab <- if ("0" %in% rownames(table_states) &  "0" %in% colnames(table_states)) table_states["0","0"] / n else 0
    
    # Calculate the observed frequencies
    pA = pAB+pAb
    pB = pAB+paB
    pa = paB+pab
    pb = pAb+pab
    
    # Reporting the normalized linkage disequilibrium D' Lewontin 1988
    if(pA*(1-pA)*pB*(1-pB) != 0) { 
      r2 <- ((pAB-pA*pB)**2)/(pA*(1-pA)*pB*(1-pB)) }
    else { r2 <-((pab-pa*pb)**2)/(pa*(1-pa)*pb*(1-pb))} 
    
    r2
  })
  
  Zns <- sum(r2values,na.rm = T)*2/(S*(S-1)) # na.rm = T, NA implies r2 is 0
  return(Zns)
}

# We need to calculate ZnS per windows
# As it looks at pairwise comparisons of snps, as the number of snps increases
# it scalates in an untractable way

ZnS_windowed <- function(pop_hap, n_window = 100) {
  P_hap <- rowSums(pop_hap)/(dim(pop_hap)[2])
  aux <- pop_hap[P_hap > 0 & P_hap < 1,]
  
  n_snps <- dim(aux)[1] 
  index_2_windows <- seq(1,n_snps,n_window) # We make indices to divide in windows of 100 snps
  
  st_time <- Sys.time()
  ZnS_per_window <- numeric(0)
  for (i in (1:(length(index_2_windows)-1))) {
    if (i == length(index_2_windows[i])-1) {
      # We calculate the last window till last snps to assure at least 100 snps in the window
      ZnS_per_window <- c(ZnS_per_window,ZnS_ld(aux[index_2_windows[i]:n_snps,]))  
    } else {
      # Calculate ZnS for the snps in position [i] to [i+1]
      ZnS_per_window <- c(ZnS_per_window,ZnS_ld(aux[index_2_windows[i]:index_2_windows[i+1],]))
    }
  }
  return(mean(ZnS_per_window))
}
