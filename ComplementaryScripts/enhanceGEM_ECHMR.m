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
function [ecModel,model_data,kcats] = enhanceGEM_ECHMR(hsa_model)
    current      = pwd;
    org_name     = 'homo sapiens';
    keggCode     = 'hsa';
    GECKO_path   = 'Write your GECKO directory here';
    %%%%%%%%%%%%%%%%%%%%%%%% model preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creates a new field in the model structure where the ENSEMBL gene IDs
    % are converted to their short gene names in order to provide
    % compatibility with the kcat matching algorithms
    hsa_model  = substitute_grRules(hsa_model); 
    hsa_model  = substituteEnsemblGeneIDs(hsa_model); 
    cd (current)
    % Standardizes the metabolites names
    hsa_model  = modifyMetNames(hsa_model);
    % Add HepG2 biomass reaction
    HepG2model          = addHepG2BiomassRxns(hsa_model,false);
    HepG2model.rxnNames = HepG2model.rxnNames;
    % Save models  
    cd ../HMR2.0
    save('HepG2.mat','hsa_model','HepG2model')
    %%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve kcats & MWs for each rxn in the model:
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_model_data = getEnzymeCodes(HepG2model);
    cd (current)
    cd ../HMR2.0
    %Save ecModel data
    save('HepG2_enzData.mat','HepG2model','HepG2_model_data');
    %Tries to Match kinetic coefficients to every reaction with a non empty
    % grRule
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_kcats  =  matchKcats(HepG2_model_data,org_name);
    cd (current)
    cd ../HMR2.0
    save('HepG2_kcats.mat','HepG2_kcats');    
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    HepG2_ecModel = readKcatData(HepG2_model_data,HepG2_kcats);
    % Save output models:
    cd (current)
    cd ../HMR2.0/EC_HMR
    save('EC_HepG2.mat','HepG2_ecModel','HepG2_model_data','HepG2_kcats')
    %%%%%%%%%%%%%%%%%%%%%%%% Matched Kcats analysis %%%%%%%%%%%%%%%%%%%%%%%
    %Gets the model Kcats cumulative distributions and compares it to all
    %the Kcat entries in BRENDA for Homo sapiens (just for natural
    %substrates)
    cd ([current '/KcatDistributions']) 
    kcat_distributions(model,HepG2_kcats,{'homo sapiens'})
    cd (current)
    %%%%%%%%%%%%%%%%%%%%%%%% Constrain enzyme pool %%%%%%%%%%%%%%%%%%%%%%%%
    % Create and constrain the total enzymes pool
    cd /Data_integration
    [Ptot,f]            = modelProteinContent(model);
    sigma               = 0.5;    %Average enzyme saturation factor
    cd ([GECKO_path '/Matlab_Module/limit_proteins'])
    HepG2_ecModel_batch = constrainEnzymes(HepG2_ecModel,Ptot,sigma,f);
    %%%%%%%%%%%%%%%%%%%%%%%% Sensitivity analysis %%%%%%%%%%%%%%%%%%%%%%%%%
    cd ([GECKO_path '/Matlab_Module/Kcat_sensitivity_analysis'])
    prevUnicodes        = [];
    exp_gRate           = 1; %[g/gDwh]
	HepG2_ecModel       = modifyKcats(HepG2_ecModel,ecModel_batch,exp_gRate);
    % Constrain the total enzyme pool once that the kinetic coefficients
    % were modified
    HepG2_ecModel_batch = constrainEnzymes(HepG2_ecModel,Ptot,sigma,f);

    cd ../HMR2.0/EC_HMR
    save('ecHepG2_batch.mat','HepG2_ecModel','HepG2_model_data',...
                                    'HepG2_kcats','HepG2_ecModel_batch')                            
                                   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%