function [ecModel,ecModel_batch] = enhanceGEM_cellLine(cellName,GECKO_path)
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
% Ivan Domenzain.      Last edited: 2019-11-20

current  = pwd;
org_name = 'homo sapiens';
%Load original model 
model = load(['../models/' cellName '/' cellName '.mat']);
eval(['model = model.' cellName])
%% model preprocessing 
model_modified = modelModifications(model);
% Save models
save(['../models/' cellName '/model.mat'],'model')
save(['../models/' cellName '/model_modified.mat'],'model_modified')
%% GECKO pipeline 
% Retrieve kcats & MWs for each rxn in the model from Uniprot database:
cd GECKO/geckomat/get_enzyme_data
model_data = getEnzymeCodes(model_modified);
%Tries to Match kinetic coefficients to every reaction with a non empty
%grRule
kcats = matchKcats(model_data,org_name);
%Save ecModel data
newDir = ['../../../../models/' cellName '/Data'];
mkdir(newDir)
save([newDir '/enzData.mat'],'model_data','kcats');
%GEt ecModel matlab structure
model_data = removeFields(model_data);
cd ../change_model
ecModel = readKcatData(model_data,kcats);
% Save output models:
save(['../../../../models/' cellName '/ecModel.mat'],'ecModel')
cd ../limit_proteins
% Constrain total protein pool
Ptotal       = 0.609; %HepG2 total protein content [g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;
[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);
save(['../../../../models/' cellName '/ecModel_batch.mat'],'ecModel_batch')
end