# Constraint-Based-Modeling
MATLAB code developed for Contraint-Based Modeling methods of Metabolic Networks. Methods are meant to be used with the MATLAB version of the [COBRA](https://opencobra.github.io/cobratoolbox/stable/) toolbox. The following files are contained in this repository.

## CORDA.m
Performs the Cost-Optimization Reaction Dependency Assessment algorithm. This algorithm tailors generalized metabolic reconstructions to context-specific metabolic models (e.g. specific cell-line or cell type) using gene or protein expression of metabolic enzymes. The algorithm is described in:

Schultz, A., & Qutub, A. A. (2016). Reconstruction of tissue-specific metabolic networks using CORDA. PLoS computational biology, 12(3), e1004808. [PMID: 26942765](https://www.ncbi.nlm.nih.gov/pubmed/26942765)

The following figure from the publication gives a brief overview of the algorithm:

![](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=10.1371/journal.pcbi.1004808.g001)

<sub>(A) Recon1 subnetwork involving water (h2o), oxygen (o2), hydrogen peroxide (h2o2) and superoxide anion (o2s) illustrating how standard oxygen (blue) and water (green) import pathways can be substituted by alternative, physiologically unlikely pathways (red and orange respectively). All metabolites and reactions are labeled as in Recon1. (B) Overview of the dependency assessment method. Each reaction in the reconstruction is associated with a specific cost through the addition of a pseudo-metabolite to the model. FBA is then performed while minimizing the cost production in order to identify high cost reactions which are favorable to the reaction being tested. (C) The CORDA tissue-specific algorithm. During each step, reaction groups being tested are outlined in blue, while reaction groups associated with a high cost are outlined in red.  </sub>

## CORDA2.m
Implements the second iteration of CORDA, which is faster and noise-independent. The algorithms is described in:

Schultz, A., Mehta, S., Hu, C. W., Hoff, F. W., Horton, T. M., Kornblau, S. M., & Qutub, A. A. (2017). Identifying cancer specific metabolic signatures using constraint-based models. In PACIFIC SYMPOSIUM ON BIOCOMPUTING 2017 (pp. 485-496). [PMID: 27897000](https://www.ncbi.nlm.nih.gov/pubmed/27897000)

## mfACHR.m
Implements the matrix-form Artificially Centered Hit and Run algorithm. When compared to the gpSampler implementation, mfACRH runs about two times faster and reduces the need for parallelization. The algorithm is described in:

Schultz, A., Mehta, S., Hu, C. W., Hoff, F. W., Horton, T. M., Kornblau, S. M., & Qutub, A. A. (2017). Identifying cancer specific metabolic signatures using constraint-based models. In PACIFIC SYMPOSIUM ON BIOCOMPUTING 2017 (pp. 485-496). [PMID: 27897000](https://www.ncbi.nlm.nih.gov/pubmed/27897000)

## corsoFBA.m
Implementation of the COst Reduced Sub-Optimal FBA algorithm. This algorithm predicts metabolic reaction fluxes in a sub-optimal space by minimizing reaction costs estimated based on protein levels and thermodynamic values. The algorithm is described in:

Schultz, A., & Qutub, A. A. (2015). Predicting internal cell fluxes at sub-optimal growth. BMC systems biology, 9(1), 18. [PMID: 25890056](https://www.ncbi.nlm.nih.gov/pubmed/25890056)

## sammif
Folder contains files to run the Semi-Automated Metabolic Map Illustrator from MATLAB. Add this folder to the MATLAB path to use this tool. To see the options for running the function type ```help sammi``` in the MATLAB command line. To test SAMMI within MATLAB run the command ```testSammi(n)``` where ```n``` ranges from zero to four. To view the code for these examples type ```edit testSammi```.

The SAMMI tool is available online at [www.sammitool.com](http://www.sammitool.com). Click [here](https://i9.ytimg.com/vi/YJ-0J4DysY4/mqdefault.jpg?sqp=COCZ5ucF&rs=AOn4CLCZ93CoPpQGO1OkSEcQqUmWnzJPSA&time=1559858445432)](https://youtu.be/YJ-0J4DysY4) for a short SAMMI tutorial:

