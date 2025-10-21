# ABCDen

This repository contains code used in the the paper: "Denisovan introgression left differential selection regimes in Humans and Neanderthals on the *SLC30A9* gene"

Models were built in order to be run in parallel, for example, I used:

seq 1 100000 \| parallel -j 20 Rscript Int_model_anc.neg.local.R /home/jorge/PROJECTS.JORGE/SLENDR/ANC.NEG.TEST AncSimulation {} 5

Where 5 is the rescaling factor Q. This models are built on slendr with slim custome script to model a dominant mutation arising in the specified lineage.

Snakefile includes the workflow used to compute the observed statistics.

For more details, see Materials and Methods.
