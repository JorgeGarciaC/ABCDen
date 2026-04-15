# ABCDen

This repository contains code used in the the paper: "Denisovan introgression left differential selection regimes in Humans and Neanderthals on the *SLC30A9* gene"

## SIMS
This folder contains the scripts used to simulate genetic data and the scripts used to calculate the statistics in the observed data with the helper script *All_statistics.R*. The pipeline used to process the observed data can be found in the *Snakefile*. This models are built on slendr with slim custome script to model a dominant mutation arising in the specified lineage.  
*ZnS_Calculation.py* and *ZnS_observed.py* includes auxiliar scripts to compute ZnS, as it is way more efficient than in R (and due to the big computational time we needed to reduce time). Models were built in order to be run in parallel, for example, I used:

`seq 1 100000 \| parallel -j 20 Rscript Int_model_anc.R /home/jorge/PROJECTS.JORGE/SLENDR/ANC.NEG.TEST AncSimulation {} 5`


## ABC_RF
This folder contains the script used to fine-tune the Random Forest, the scripts used in the ABC analysis and auxiliar script with processing data functions. 
*RF_4models.NEG.R* includes the code used to compare neural networks versus Random Forest and to model the Random Forest and perform the feature extraction.
*ABC_4MODELS.NEG.R* includes the code of the ABC model comparison and ABC parameter inference using the extracted features from the random forest.

## HApsigN.R

Includes the developed haplotype painting algorithm. 


For more details, see Materials and Methods.
