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
function [ecModel,model_data,kcats] = enhanceGEM_ECHMR(HMR_modified)
    current      = pwd;
    org_name     = 'homo sapiens';
    keggCode     = 'hsa';
    GECKO_path   = 'Write your GECKO directory here';
    %%%%%%%%%%%%%%%%%%%%%%%% model preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creates a new field in the model structure where the ENSEMBL gene IDs
    % are converted to their short gene names in order to provide
    % compatibility with the kcat matching algorithms
    HMR_modified  = substitute_grRules(HMR_modified); 
    HMR_modified  = substituteEnsemblGeneIDs(HMR_modified); 
    cd (current)
    % Standardizes the metabolites names
    HMR_modified  = modifyMetNames(HMR_modified);
    % Substitute biomass reaction
    HMR_modified  = substituteBiomassRxns(HMR_modified,false);
    % Save models  
    cd ../models
    save('HMR_modified.mat','HMR_modified')
    %%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve kcats & MWs for each rxn in the model:
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HMR_model_data = getEnzymeCodes(HMR_modified);
    cd (current)
    cd ../models/Data
    %Save ecModel data
    save('HMR_enzData.mat','HMR_modified','HMR_model_data');
    %Tries to Match kinetic coefficients to every reaction with a non empty
    % grRule
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HMR_kcats  =  matchKcats(HMR_model_data,org_name);
    cd (current)
    cd ../models/Data
    save('HMR_kcats.mat','HMR_kcats');    
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    HMR_ecModel = readKcatData(HMR_model_data,HMR_kcats);
    % Save output models:
    cd (current)
    cd ../models/EC_HMR
    save('HMR_ecModel.mat','HMR_ecModel')
    %%%%%%%%%%%%%%%%%%%%%%%% Matched Kcats analysis %%%%%%%%%%%%%%%%%%%%%%%
    %Gets the model Kcats cumulative distributions and compares it to all
    %the Kcat entries in BRENDA for Homo sapiens (just for natural
    %substrates)
    cd ([current '/KcatDistributions']) 
    kcat_distributions(HMR_ecModel,HMR_kcats,{'homo sapiens'})
    cd (current)
    %%%%%%%%%%%%%%%%%%%%%%%% Constrain enzyme pool %%%%%%%%%%%%%%%%%%%%%%%%
    % Create and constrain the total enzymes pool
    cd /Data_integration
    [Ptot,f]            = modelProteinContent(model);
    sigma               = 0.5;    %Average enzyme saturation factor
    cd ([GECKO_path '/Matlab_Module/limit_proteins'])
    HMR_ecModel_batch = constrainEnzymes(HMR_ecModel,0.Ptot,0.5);%,sigma,f);
    %%%%%%%%%%%%%%%%%%%%%%%% Sensitivity analysis %%%%%%%%%%%%%%%%%%%%%%%%%
    cd ([GECKO_path '/Matlab_Module/Kcat_sensitivity_analysis'])
    prevUnicodes      = [];
    exp_gRate         = 1; %[g/gDwh]
	HMR_ecModel       = modifyKcats(HMR_ecModel,HMR_ecModel_batch,exp_gRate);
    % Constrain the total enzyme pool once that the kinetic coefficients
    % were modified
    cd ([GECKO_path '/Matlab_Module/limit_proteins'])
    HMR_ecModel_batch = constrainEnzymes(HMR_ecModel,0.Ptot,0.5);%,sigma,f);
    cd ../models/EC_HMR
    save('HMR_ecModel_batch.mat','HMR_ecModel_batch')                            
                                   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%