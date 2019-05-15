function [ecModel,ecModel_batch] = enhanceGEM_cellLine(cellName)
%enhanceGEM_cellLine
% 
% Function that loads a cell-line or tissue specific human metabolism model
% and enhances it with enzyme-constraints with the use of the GECKO
% pipeline.
%
%   cellName        (String) name for the context-specific model. It should 
%                   be consistent with the name of the subfolder in which  
%                   the original model is stored. 
%
%   ecModel         Enzyme-constrained model structure without proteomics
%                   constraints.
%   ecModel_batch   Enzyme-constrained model structure with a constrained
%                   total protein pool.
%
% Usage: [ecModel,ecModel_batch] = enhanceGEM_cellLine(cellName)
%
% Ivan Domenzain.      Last edited: 2019-05-15

current  = pwd;
org_name = 'homo sapiens';
%Load original model 
model    = load(['../models/' cellName '/' cellName '.mat']);
eval(['model = model.' cellName])
%model = model.model;
%%%%%%%%%%%%%%%%%%%%%%%% model preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%
%model_modified = ravenCobraWrapper(model);
model_modified = modelModifications(model);
% Save models
cd (['../models/' cellName])
save('model.mat','model')
save('model_modified.mat','model_modified')
%%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve kcats & MWs for each rxn in the model from Uniprot database:
cd ([GECKO_path '/geckomat/get_enzyme_data'])
model_data = getEnzymeCodes(model_modified);
%Tries to Match kinetic coefficients to every reaction with a non empty
%grRule
kcats = matchKcats(model_data,org_name);
%Save ecModel data
cd (current)
cd (['../models/' cellName])
mkdir Data
cd Data
save('enzData.mat','model_data','kcats');
cd (current)
%GEt ecModel matlab structure
model_data = removeFields(model_data);
cd ([GECKO_path '/geckomat/change_model'])
ecModel = readKcatData(model_data,kcats);
% Save output models:
cd (current)
cd (['../models/' cellName])
save('ecModel.mat','ecModel')
cd ([GECKO_path '/geckomat/limit_proteins'])
% Constrain total protein pool
Ptotal       = 0.609; %HepG2 total protein content [g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;
[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);
cd (['../../../../models/' cellName])
save('ecModel_batch.mat','ecModel_batch')
cd (current)
end