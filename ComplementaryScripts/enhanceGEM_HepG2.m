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
% Ivan Domenzain.      Last edited: 2018-01-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = enhanceGEM_HepG2(HepG2model)
    current      = pwd;
    org_name     = 'homo sapiens';
    keggCode     = 'hsa';
    GECKO_path   = 'Write your GECKO directory here';
    %%%%%%%%%%%%%%%%%%%%%%%% model preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creates a new field in the model structure where the ENSEMBL gene IDs
    % are converted to their short gene names in order to provide
    % compatibility with the kcat matching algorithms
    HepG2model_modified = modelModifications(HepG2model);
    % Save models  
    cd ../models
    save('HepG2model.mat','HepG2model')
    save('HepG2model_modified.mat','HepG2model_modified')
    %%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve kcats & MWs for each rxn in the model:
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_model_data = getEnzymeCodes(HepG2model_modified);
    %Tries to Match kinetic coefficients to every reaction with a non empty
    %grRule
    HepG2_kcats = matchKcats(HepG2_model_data,org_name);
    %Save ecModel data
    cd (current)
    cd ../models/HepG2/Data
    save('HepG2_enzData.mat','HepG2_model_data','HepG2_kcats');
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    HepG2_ecModel = readKcatData(HepG2_model_data,HepG2_kcats);
    % Save output models:
    cd (current)
    cd ../models/HepG2/ecHepG2
    save('HepG2_ecModel.mat','HepG2_ecModel')
    %%%%%%%%%%%%%%%%%%%%%%%% Matched Kcats analysis %%%%%%%%%%%%%%%%%%%%%%%
    %Gets the model Kcats cumulative distributions and compares it to all
    %the Kcat entries in BRENDA for Homo sapiens (just for natural
    %substrates)
    cd ([current '/KcatDistributions']) 
    kcat_distributions(HepG2_ecModel,HepG2_kcats,{'homo sapiens'})
    cd (current)
    %%%%%%%%%%%%%%%%%%%%%%%% Constrain enzyme pool %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Sensitivity analysis %%%%%%%%%%%%%%%%%%%%%%%%%
                  
                                   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%