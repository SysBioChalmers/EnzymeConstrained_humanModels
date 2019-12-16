# Enzyme-constrained Human1 GEMs
This directory contains an automated pipeline for constructing cell-specific enzyme-constrained GEMs (ecGEMs) derived from Human-GEM (v1.3.0) based on transcriptomics and proteomics datasets.

### Required Software:
* A functional Matlab installation (MATLAB 7.3 and higher)
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN) (version 2.3.0)
* For generating some figures, a functional [R installation](https://www.r-project.org/) (version 3.6.1)

### Dependencies - Recommended Software:
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended)
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulation

## Regenerating the ecGEMs:
The ecGEMs are already present in the `models/` subdirectory, but the scripts and data necessary to regenerate the models are available here. The master script for generating the ecGEMs from the tINIT GEMs (provided in the `/models/humanGEM_cellLines/11models.mat` file) is `generate_human_ecModels_NCI60.m`, located in the `/ComplementaryScripts` directory. Run this `generate_human_ecModels_NCI60` script in MATLAB to regenerate the 11 ecGEMs.

## Flux variability analysis
Flux variability analysis (FVA) (corresponding to the results presented in Fig. 5B) can be run using the `comparativeFVA_humanModels.m` function in the `ComplementaryScripts/Simulation` subdirectory. Specify the name of the model (cell line) for which FVA is to be run; for example:

`results = comparativeFVA_humanModels('HOP92');`

## Prediction of growth and metabolite exchange rates:
To use the ecGEMs and non-ecGEMs to predict growth rates and metabolite exchange rates with increasing levels of constraints (as shown in Figs. 5C and 5D), run the `predict_cellLines_gRates.m` script in the `ComplementaryScripts` subdirectory.


## Generating plots for Figure 5
The R script used to generate plots shown in Fig. 5 of the main text is `plot_ecGEM_results.R`, located in the `ComplementaryScripts` subdirectory.
