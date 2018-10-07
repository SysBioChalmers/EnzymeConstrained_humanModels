%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
%
% Script that takes the already existing HMR2.0 GEM and extends it to an 
% Enzyme constrained GEM. The previous requeriments for the implementation
% of this script are the KEGG and Uniprot databases for H. sapiens,
% detailed on the GECKO's paper supplementary material and a list of
% equivalences between ENSEMBL human gene IDs and its correspondent short
% gene names.
%
% INPUT:
%   model       The latest version of HMR2.0 genome scale metabolic model
% OUTPUTS:
%   ecModel     Extended enzyme constrained model
%
% Ivan Domenzain.      Last edited: 2018-10-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ecModel,model_data,kcats] = enhanceHumanGEM(model,cellName)
current      = pwd;
org_name     = 'homo sapiens';
keggCode     = 'hsa';
GECKO_path   = '/Users/ivand/Documents/GitHub/GECKO';
%%%%%%%%%%%%%%%%%%%%%%%% model preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a new field in the model structure where the ENSEMBL gene IDs
% are converted to their short gene names in order to provide
% compatibility with the kcat matching algorithms
%model_modified = ravenCobraWrapper(model);
model_modified = modelModifications(model);
% Save models
cd (['../models/' cellName])
save('model.mat','model')
save('model_modified.mat','model_modified')
%%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve kcats & MWs for each rxn in the model:
cd ([GECKO_path '/geckomat/get_enzyme_data'])
model_data = getEnzymeCodes(model_modified);
%Tries to Match kinetic coefficients to every reaction with a non empty
%grRule
kcats = matchKcats(model_data,org_name);
%Save ecModel data
cd (current)
cd (['../models/' cellName '/Data'])
save('enzData.mat','model_data','kcats');
% Integrate enzymes in the model: ecModel is a draft model for proteomics
% data integration, if available.
cd (current)
model_data = removeFields(model_data);
cd ([GECKO_path '/geckomat/change_model'])
ecModel = readKcatData(model_data,kcats);
% Save output models:
cd (current)
cd (['../models/' cellName])
save('ecModel.mat','ecModel')
cd ([current '/GECKO'])
% Constrain total protein pool
Ptotal       = 0.5; %[g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;
[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);
cd (['../models/' cellName])
save('ecModel_batch.mat','ecModel_batch')
%%%%%%%%%%%%%%%%%%%%%%%% Matched Kcats analysis %%%%%%%%%%%%%%%%%%%%%%%
%Gets the model Kcats cumulative distributions and compares it to all
%the Kcat entries in BRENDA for Homo sapiens (just for natural
%substrates)
cd ([current '/KcatDistributions'])
kcat_distributions(ecModel,kcats,{'homo sapiens'})
cd (current)
%%%%%%%%%%%%%%%%%%%%%%%% Constrain enzyme pool %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Sensitivity analysis %%%%%%%%%%%%%%%%%%%%%%%%%

                                   
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%