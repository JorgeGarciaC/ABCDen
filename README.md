# ABCDen

This repository contains code used in the the paper: "Denisovan introgression left differential selection regimes in Humans and Neanderthals on the *SLC30A9* gene"

Models were built in order to be run in parallel, for example, I used:

seq 1 100000 \| parallel -j 20 Rscript Int_model_anc.R /home/jorge/PROJECTS.JORGE/SLENDR/ANC.NEG.TEST AncSimulation {} 5

where 5 is the rescaling factor Q. This models are built on slendr with slim custome script to model a dominant mutation arising in the specified lineage.

-   **Snakefile:** Includes the workflow used to compute the observed statistics.

-   **RF_4models.NEG.R:** Includes the code used to compare neural networks versus Random Forest and to model the Random Forest and perform the feature extraction

-   **ABC_4MODELS.NEG.R:** Includes the code of the ABC model comparison and ABC parameter inference using the extracted features from the random forest

-   **ZnS_Calculation.py:** Includes an auxiliar script to compute ZnS, as it is way more efficient than in R (and due to the big computational time we needed to reduce time)

For more details, see Materials and Methods.
