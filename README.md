# EnzymeConstrained-human-GEM

Series of scripts for enhancing humanGEM based models with kinetic and proteomics constraints, and specialized simulation utilities. 

This repository contains an automated pipeline for constructing ecModels of human metabolism based on humanGEM for specific cell lines or tissues based on transcriptomics and proteomics datasets.

The generic genome-scale metabolic model of _Homo sapiens_ is available at:
https://github.com/SysBioChalmers/human-GEM

Keywords:

**GEM Category:** Cell lines/tissues; **Utilisation:** Proteomics data integration and simulation of human metabolism; **Field:** Enzyme-constrained metabolic network reconstruction; **Type of Model:** curated; **Model Source:** Hsa; **Omic Source:** Transcriptomics, Proteomics, Metabolomics; **Taxonomy:** _Homo sapiens_; **Metabolic System:** General Metabolism.

Last update: 2019-05-15

This repository is administered by [@IVANDOMENZAIN](https://github.com/IVANDOMENZAIN), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation

### Required Software:
* A functional Matlab installation (MATLAB 7.3 and higher).
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* Clone the model from [master](https://github.com/SysBioChalmers/) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers)

## Contributors
- Iván Domenzain [(@IVANDOMENZAIN)](https://github.com/IVANDOMENZAIN), Chalmers University of Technology, Gothenburg Sweden
- Raphaël Ferreira [(@raphDL)](https://github.com/raphDL), Chalmers University of Technology, Gothenburg Sweden
