library(slendr,quietly = T, warn.conflicts = F)
library(ggplot2,quietly = T, warn.conflicts = F)
library(tidyverse,quietly = T, warn.conflicts = F)
library(scales,quietly = T, warn.conflicts = F)
library(reticulate,quietly = T, warn.conflicts = F)
library(logr,quietly = T, warn.conflicts = F)

rm(list=ls())

start_time <- Sys.time()
source("~/Rscripts/All_statistics.R")
args <- commandArgs(trailingOnly = T)

directory_name <- args[1] #"/home/jorge/PROJECTS.JORGE/SLENDR" #
name_ts <- args[2] #"twodn" 
replicate <- args[3] #1 
Q <- as.numeric(args[4]) # Reescaling factor)

save.statistics.file <- paste0("Q",Q,".statistics.tsv")
setwd(directory_name)
lf <- log_open(paste0(name_ts,"_",replicate, ".log"))

###### Extension template ######
extension_template <- r"(
initialize() {
    defineConstant("s_den", {{s_den}});
    defineConstant("s_nea", {{s_nea}});
    defineConstant("s_eur", {{s_eur}});
    defineConstant("s_chb", {{s_chb}});
    defineConstant("s_ooa", {{s_ooa}});
    defineConstant("onset_time", asInteger({{onset_time}}));
    defineConstant("target_pop", "{{target_pop}}");
    defineConstant("origin_pop", "{{origin_pop}}");
    defineConstant("output_file", "{{directory_to_save}}/{{replicate_n}}_traj_" + target_pop + "_" + origin_pop + ".tsv"); 
}

initialize() {
    initializeMutationType("m1", 0.5, "f", s_den); 
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, SEQUENCE_LENGTH);
    initializeMutationRate(0);
    initializeRecombinationRate(RECOMBINATION_RATE);
}

function (void) add_mutation(void) {
    // add a beneficial mutation of the given selection strength to a random genome of the
    // population with the name `origin_pop`
    target = sample(population(origin_pop).haplosomes, 1);
    mut = target.addNewDrawnMutation(m1, position = asInteger(SEQUENCE_LENGTH/2));
    defineGlobal("MUTATION", mut);

    write_log("adding beneficial mutation to population " + origin_pop);

    // write out a header of an output file
    writeFile(output_file, "time\tfreq_ANC\tfreq_NEA\tfreq_DEN\tfreq_VIN\tfreq_Eafr\tfreq_OOA\tfreq_CHB\tfreq_EUR");
}

// add a beneficial mutation at a given onset time
tick(onset_time) late() {
    add_mutation();
}

// make selection coefficient different for each population
POPULATIONS.getValue("tsplit_gen")[1]:POPULATIONS.getValue("tremove_gen")[1] mutationEffect(m1,p1){ return 1.0 + s_nea; } // nea pop
//POPULATIONS.getValue("tsplit_gen")[2]:POPULATIONS.getValue("tremove_gen")[2] mutationEffect(m1,p2){ return 1.0 + s_den; }
POPULATIONS.getValue("tsplit_gen")[3]:POPULATIONS.getValue("tremove_gen")[3] mutationEffect(m1,p3){ return 1.0 + s_nea; } // vindija pop
POPULATIONS.getValue("tsplit_gen")[4]:SIMULATION_END mutationEffect(m1,p4){ return 1.0; } 
POPULATIONS.getValue("tsplit_gen")[5]:POPULATIONS.getValue("tremove_gen")[5] mutationEffect(m1,p5){ return 1.0 + s_ooa; }
POPULATIONS.getValue("tsplit_gen")[6]:SIMULATION_END mutationEffect(m1,p6){ return 1.0 + s_chb; }
POPULATIONS.getValue("tsplit_gen")[7]:SIMULATION_END mutationEffect(m1,p7){ return 1.0 + s_eur; }

// check that the mutation is not lost in every tick after it appeared in the simulation
tick(onset_time):SIMULATION_END late() {
    // the mutation is not segregating and is not fixed either -- we must restart
    if (!MUTATION.isSegregating & !MUTATION.isFixed) {
        write_log("mutation lost -- finishing simulation");
        if (TS_PATH != "")
          defineConstant("OUTPUT_TS", TS_PATH + "/output.trees");
        else if (TS_PATH != "")
          defineConstant("OUTPUT_TS", TS_PATH);
        else // this should never happen because the block would be deregistered
          stop("No tree-sequence output path has been set\n");

          //save_ts(OUTPUT_TS);
        sim.simulationFinished();
    }

    // compute the frequency of the mutation of interest and save it (if the
    // mutation is missing at this time, save its frequency as NA)
    freq_anc = "NA";
    freq_nea = "NA";
    freq_den = "NA";
    freq_vin = "NA";
    freq_eafr = "NA";
    freq_ooa = "NA";
    freq_chb = "NA";
    freq_eur = "NA";
    if (population("ANC", check = T))
      freq_anc = population("ANC").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("NEA", check = T))
      freq_nea = population("NEA").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("DENI", check = T))
      freq_den = population("DENI").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("VIN", check = T))
      freq_vin = population("VIN").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("Eafr", check = T))
      freq_eafr = population("Eafr").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("OOA", check = T))
      freq_ooa = population("OOA").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("CHB", check = T))
      freq_chb = population("CHB").haplosomes.mutationFrequenciesInHaplosomes();
    if (population("EUR", check = T))
      freq_eur = population("EUR").haplosomes.mutationFrequenciesInHaplosomes();

    writeFile(output_file,
              model_time(community.tick) + "\t" +
              freq_anc+ "\t" +
              freq_nea + "\t" +
              freq_den + "\t" + 
              freq_vin + "\t" +
              freq_eafr + "\t" + 
              freq_ooa + "\t" +
              freq_chb + "\t" + 
              freq_eur, append = T);
})"

##### Initialize ####
init_env(quiet = T)
msp <<- reticulate::import("msprime", delay_load = TRUE) # to acces late functionlaties
origin_pop <- "DENI"
target_pop <- "OOA"
onset_time <- runif(1, min = 250e3, max = 400e3) # Time arising the mutation
extra_chunk <- 50e3 # extra chunk to avoid edge effects on the end of the chunk
l <- 70000 # sequence length
simulated_chunk <- 310000
# generate selection coefficients for each population
s_den <- runif(1,min = 0.0001, max = 0.05) * Q
s_nea <- runif(1,min = -0.05, max = 0.05) * Q
s_eur <- runif(1,min = 0.0001, max = 0.05) * Q
s_chb <- runif(1,min = 0.0001, max = 0.05) * Q
s_ooa <- runif(1,min = 0.0001, max = 0.05) * Q
rrate <- 1e-8*Q #
mrate <- 1e-8*Q  # 2.36e-8 [38],
output_n <- paste0(name_ts,"_",replicate,".ts")
MODEL <- "SLiM"

########################
# NEANDERTAL ORIGIN MODEL
########################
ANC <- population("ANC", time = 650e3, N = 30000 /Q, remove = 40e3) # Papuans out of africa (human + arcaic)
E_afr <- population("Eafr", parent = ANC, time = 90e3, N = 25e3 /Q) # Times from Serradell et al, Ne two stems
OOA <- population("OOA", parent = E_afr, time = 70e3, N = 8500 /Q, remove = 32e3) # Papuans out of africa
eur <- population("EUR", parent = OOA, time = 37e3, N = 7000 / Q) # Papuans out of africa
chb <- population("CHB", parent = OOA, time = 44e3, N = 9000 /Q) # Papuans out of africa
nea <- population("NEA", parent = ANC, time = 600e3, N = 13000 /Q, remove = 40e3) # Papuans (deni + nean)
den <- population("DENI", parent = nea, time = 400e3, N = 1000 /Q, remove = 25e3) #Choin et al 2021
vin <- population("VIN", parent = nea, time = 130e3, N = 1000 /Q, remove = 40e3) #

# Neanderthal population splitting at 600 ky ago from modern humans(becomes extinct by 40 ky ago)
gf <- list(gene_flow(from = eur, to = E_afr, start = 13e3, end = 12e3, rate = 0.12), # BACK TO AFRICA
           gene_flow(from = den, to = nea, start = 300e3, end = 60e3, rate = 0.005 ), # review peyregne
           gene_flow(from = nea, to = den, start = 300e3, end = 60e3, rate = 0.005 ),
           gene_flow(from = eur, to = chb, start = 36e3, end = 0, rate = 1-((1-3.14e-5)^(36e3/29))), # Change from migration rates to admixture fraction
           gene_flow(from = chb, to = eur, start = 36e3, end = 0, rate = 1-((1-3.14e-5)^(36e3/29))), # Stpopsim Pap
           gene_flow(from = nea, to = OOA, start = 47e3, end = 46e3, rate = 0.02)) 

extension <- substitute_values(
  extension_template,
  s_den = s_den, s_nea = s_nea, s_chb = s_chb, s_eur = s_eur, s_ooa = s_ooa,
  onset_time = onset_time,
  origin_pop = origin_pop, target_pop = target_pop,
  directory_to_save = directory_name, replicate_n = replicate)

stem_model <- compile_model(
  populations = list(
    ANC,nea,vin,den,OOA,E_afr, eur, chb),
  gene_flow = gf, generation_time = 29*Q, #We do this in order to not need to reescale all the years
  # Briefly, years/generation_time = Nºgenerations; years/(generation_time*Q) = Nºgenerations/Q; 
  path = paste0(tempfile(), "_introgression"),extension = extension
)

Altai_samples <- schedule_sampling(stem_model, times = 70000, list(nea, 1)) # 
present_samples <- schedule_sampling(stem_model, times = 0, 
                                     list(chb,100), list(E_afr,100), list(eur, 100))
den_samples <- schedule_sampling(stem_model, times = 60000, list(den,1))
char_samples <- schedule_sampling(stem_model, times = 80e3, list(vin,1))
vin_samples <- schedule_sampling(stem_model, times = 52000, list(vin,1))
old_deni <- schedule_sampling(stem_model, times = 200e3, list(den,1))

ts <- try( slim(stem_model, 
                sequence_length = simulated_chunk, recombination_rate = rrate,
                samples = rbind(Altai_samples, present_samples, den_samples,vin_samples,old_deni,char_samples),
                path = output_n) )

print(Sys.time()-start_time)
if (inherits(ts, 'try-error')) {
  if (ts[1] == "Error : Unfortunately SLiM terminated before a tree sequence was saved.\nSee the above for an indication of where things could have gone wrong.\n")
    log_print(paste0("Expected behaviour: Mutation lost, simulation number ",replicate," was cancelled"))
  write.table(t(c(rep(-Inf,53),replicate)), file = save.statistics.file, quote = F, sep = "\t", append = T,row.names = F, col.names = F)
  mutation <- "0"
} else {
  ts <- ts_read(paste0(output_n,"/slim.trees"),stem_model)
  #### Tree Sequence processing
  if (!ts_coalesced(ts)) {
    ts <- ts_recapitate(ts, recombination_rate = rrate, Ne = 30000 / Q  )
  } # Check is coalesced, is false recapitate
  
  ts_small <- ts_simplify(ts) # simplify to eliminate non informative trees.
  # must be the same as the one specified in the extension
  ts_m <- ts_mutate(ts_small, mutation_rate = mrate, mutation_model = msp$SLiMMutationModel(type=1L)) 
  ts_save(ts_m,paste0(name_ts,"_",replicate,"_", "mut",".ts") )
  ts_vcf(ts_m, paste0(name_ts, "_", replicate, "_","mut",".vcf.gz"), chrom = "4")
  
  log_print(paste0("Simulation number ", replicate," ended ---- Starting computing statistics"))
  gt <- ts_genotypes(ts_m) 
  gt <- gt[gt$pos >= 120000 & gt$pos <= 120000+70000,] # To calculate the CORE position
  
  CHB <- gt %>% dplyr::select(starts_with("CHB")) %>% as.matrix()
  EUR <- gt %>% dplyr::select(starts_with("EUR")) %>% as.matrix()
  EAF <- gt %>% dplyr::select(starts_with("Eafr")) %>% as.matrix()
  NEA <- gt %>% dplyr::select(starts_with("NEA")) %>% as.matrix()
  DEN <- gt %>% dplyr::select(starts_with("DEN")) %>% as.matrix()
  
  popsize_EAF <- (dim(EAF)[2])/2
  popsize_CHB <- (dim(CHB)[2])/2
  popsize_EUR <- (dim(EUR)[2])/2
  popsize_DEN <- (dim(DEN)[2])/2
  popsize_NEA <- (dim(NEA)[2])/2
  
  Fr_EAF <- rowSums(EAF)/(dim(EAF)[2]) # Unadmixed P1
  Fr_CHB <- rowSums(CHB)/(dim(CHB)[2]) # admixed P2
  Fr_EUR <- rowSums(EUR)/(dim(EUR)[2]) # admixed P2
  Fr_DEN <- rowSums(DEN)/(dim(DEN)[2]) # Archaic introgressed P3
  Fr_NEA <- rowSums(NEA)/(dim(NEA)[2])
  P0 <- 0 
  
  #### Calculate stats
  # Introgression stats
  DfD_CHB_DEN <- getD_stats(Fr_EAF,Fr_CHB,Fr_DEN)
  DfD_EUR_DEN <- getD_stats(Fr_EAF,Fr_EUR,Fr_DEN)
  U_CHB_DEN <- extract.U(Fr_EAF,Fr_CHB,Fr_DEN,0.1,0.2) # U.1.20.100
  U_EUR_DEN <- extract.U(Fr_EAF,Fr_EUR,Fr_DEN,0.1,0.2) # U.1.20.100
  DfD_CHB_NEA <- getD_stats(Fr_EAF,Fr_CHB,Fr_NEA)
  DfD_EUR_NEA <- getD_stats(Fr_EAF,Fr_EUR,Fr_NEA)
  U_CHB_NEA <- extract.U(Fr_EAF,Fr_CHB,Fr_NEA,0.1,0.2) # U.1.20.100
  U_EUR_NEA <- extract.U(Fr_EAF,Fr_EUR,Fr_NEA,0.1,0.2) # U.1.20.100
  
  # Segregating sites
  S_EAF <- S_sites(Fr_EAF)
  S_CHB <- S_sites(Fr_CHB)
  S_EUR <- S_sites(Fr_EUR)
  S_DEN <- S_sites(Fr_DEN)
  S_NEA <- S_sites(Fr_NEA)
  PS_EAF <- sum(P_seg(Fr_EAF))/S_EAF
  PS_CHB <- sum(P_seg(Fr_CHB))/S_CHB
  PS_EUR <- sum(P_seg(Fr_EUR))/S_EUR
  
  # Genetic divergence
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
  
  # Nucleotide diversity
  tw_EAF <- thetaW(S_EAF,popsize_EAF)/l
  tw_CHB <- thetaW(S_CHB,popsize_CHB)/l
  tw_EUR <- thetaW(S_CHB,popsize_EUR)/l
  tp_EAF <- thetapi(PS_EAF,popsize_EAF)/l
  tp_CHB <- thetapi(PS_CHB,popsize_CHB)/l
  tp_EUR <- thetapi(PS_EUR,popsize_EUR)/l
  # Selection statistics
  Hs_CHB <- get_H12(CHB)
  Hs_EUR <- get_H12(EUR)
  Hs_EAF <- get_H12(EAF)
  #Linkage disequilibrium
  vcf_from_ts <- paste0(directory_name,"/",name_ts,"_",replicate, "_mut.vcf.gz")
  reticulate::source_python("~/PROJECTS.JORGE/TREE_SEQUENCES/SCRIPTS/ZnS_calculation.py", convert = T)
  pos_init <- toString(120000)
  pos_finit <- toString(120000+70000)
  ZnS_values <- Calculate_ZnS_python(vcf_from_ts, init_pos = pos_init, final_pos = pos_finit)
  ZnS_CHB <- ZnS_values[[1]]
  ZnS_EUR <- ZnS_values[[2]]
  ZnS_EAF <- ZnS_values[[3]]
  ZnS_DEN <- ZnS_values[[4]]
  ZnS_NEA <- ZnS_values[[5]]
  # Final mutation frequency
  fr <- read.table(paste0(replicate,"_traj_",target_pop,"_",origin_pop,".tsv"),header = T,sep = "\t")
  final_fr <- fr[which.min(fr$time),]
  fr_al_OOA <- mean(fr[fr$time <= 50000 & fr$time >= 35000,"freq_OOA"])
  fr_al_nea <- (mean(fr[fr$time <= 120000 & fr$time >= 40000,"freq_NEA"])+mean(fr[fr$time <= 120000 & fr$time >= 40000,"freq_VIN"]))/2
  if (is.na(final_fr[,9])) fr_allel_EUR <- 0 else fr_allel_EUR <- final_fr[,9]
  if (is.na(final_fr[,6])) fr_allel_EAF <- 0 else fr_allel_EAF <- final_fr[,6]
  if (is.na(final_fr[,8])) fr_allel_CHB <- 0 else fr_allel_CHB <- final_fr[,8]
  adaptive_snp <- gt[gt$pos == simulated_chunk/2,]
  presence_deni_old <- sum(adaptive_snp[,2:3]/2)
  presence_charg <- sum(adaptive_snp[,4:5]/2)
  presence_altai <- sum(adaptive_snp[,6:7]/2)
  presence_deni <- sum(adaptive_snp[,8:9]/2)
  presence_vin <- sum(adaptive_snp[,10:11]/2)
  # Print to table
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
                            fr_al_OOA,fr_al_nea
                            ,presence_deni_old, presence_charg,presence_altai,presence_deni,presence_vin
                            ,ZnS_CHB,ZnS_EUR,ZnS_EAF,ZnS_DEN,ZnS_NEA,replicate
  )
  write.table(t(vector.of.statistics), file = save.statistics.file, quote = F, sep = "\t", append = T,row.names = F, col.names = F) 
  mutation <- "1"
}

#### Saving parameters
end_time <- Sys.time()
time_simulation <- difftime(end_time, start_time, units = "secs")
parameter.df <- data.frame(s_den/Q,s_nea/Q,s_chb/Q,s_ooa/Q,s_eur/Q,onset_time,replicate,mutation,time_simulation)
save.df <- paste0(name_ts,".Q",Q, ".parameters.tsv")

write.table(parameter.df, file = save.df, 
            quote = F, sep = "\t", append = T,row.names = F, col.names = F)
log_print(paste0("Simulation number ", replicate," finished without problems in ", time_simulation, " seconds. Finishing!"))
