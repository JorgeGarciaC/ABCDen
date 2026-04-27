# ABCDen

This repository contains code used in the the paper: "Denisovan introgression 
left differential selection regimes in Humans and Neanderthals on the *SLC30A9* gene"

## Dependencies

*SLiM.v5.0* can be download and installed from: https://messerlab.org/slim/
To install R package *slendr.1.1.0* we recommend to follow the available 
instructions: https://bodkan.net/slendr/articles/vignette-00-installation.html
slender will create a miniconda environment for the cross-platform utilities. We 
use reticulate to load a python script to calculate the population level linkage 
disequilibrium statistic ZnS. To allow so, *scikit-allele.1.3.13* must be 
installed in the *slendr* miniconda environment: 
"Python-3.12_msprime-1.3.4_tskit-0.6.4_pyslim-1.0.4_tspop-0.0.2". 
The scripts under the ABC folder requires the following R libraries: *ranger*, 
*caret*, *tidyverse*, *gridExtra* and *abc*.

## SIMS
This folder contains the scripts used to simulate genetic data and the scripts 
used to calculate the statistics in the observed data with the helper script 
*All_statistics.R*. The pipeline used to process the observed data can be found 
in the *Snakefile*. This models are built on *slendr* with slim custome script 
to model a dominant mutation arising in the specified lineage. *ZnS_Calculation.py* 
and *ZnS_observed.py* includes auxiliar scripts to compute ZnS, as it is way 
more efficient than in R (and due to the big computational time we needed to 
reduce time). Models were built in order to be run in parallel, for example, I used:

`seq 1 100000 | parallel -j 20 Rscript Int_model_anc.R /home/jorge/PROJECTS.JORGE/SLENDR/ANC.NEG.TEST AncSimulation {} 5`

This will the following output files: 
- SLiM simulations in tree format (".ts").
- SLiM simulations in VCF format (".vcf.gz"). 
- A ".tsv" file format with the allele frequency trajectory of the mutation under selection
- A ".tsv" file format with the parameters of each simulation.
- A ".tsv" file with the calculated statistics.


## ABC_RF
This folder contains the script used to fine-tune the Random Forest, the scripts
used in the ABC analysis and auxiliar script with processing data functions.
*RF_4models.NEG.R* includes the code used to compare neural networks versus 
Random Forest and to model the Random Forest and perform the feature extraction. 
Input files for this script are the training + evaluation parameters and 
statistics from simulations. The script writes the following output files: the 
hyperparameter grid performance in ".txt", accuracy of the cross-validation in 
".txt" format, and oxpots showing the best features for each variable. 
*ABC_4MODELS.NEG.R* includes the code of the ABC model comparison and ABC 
parameter inference using the extracted features from the random forest. Input 
files for this script are the test parameters and statistics from simulations. 
The script outputs the pca plots used in ABC validation.

## HApsigN.R

Includes the developed haplotype painting algorithm. As input, the script 
expects a haplotype matrix of 0,1 (reference, alternative) in which each column 
represents a polymorphism and each row an individual. 

## HaplotypeStructure.ipynb & Correlation_Eurasian.R

Scripts used to generate Figure 2b and Supplementary Figure S2

For more details, see Materials and Methods.
