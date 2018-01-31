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
    GECKO_path   = '/Users/ivand/Desktop/GECKO2';
    Protdatabase = 'hsa_ProtDatabase.mat';
    
    toolbox      = 'COBRA';
    % COBRA toolbox required
        initCobraToolbox;
    cd ../HMR2.0
    
    % Update protein databases (KEGG and uniprot)
    cd ../ComplementaryScripts
    DB_path =  '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM/Databases';
    updateDatabases(GECKO_path,DB_path,keggCode);
    cd (current)
    
    % Creates a new field in the model structure where the ENSEMBL gene IDs
    % are converted to their short gene names in order to provide
    % compatibility with the kcat matching algorithms
    hsa_model  = substitute_grRules(hsa_model); 
    hsa_model  = substituteEnsemblGeneIDs(hsa_model); 
    cd (current)
    % Standardizes the metabolites names
    hsa_model  = modifyMetNames(hsa_model);

    % Add HepG2 biomass reaction
    HepG2model = addHepG2BiomassRxns(hsa_model,false);
    cd ../HMR2.0
    save('HepG2.mat','hsa_model','HepG2model');

    % Retrieve kcats & MWs for each rxn in model:
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_model_data = getEnzymeCodes(HepG2model);
    cd (current)
    cd ../HMR2.0
    save('HepG2_enzData.mat','HepG2model','HepG2_model_data');
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_kcats      =  matchKcats(HepG2_model_data,org_name);
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
    [Ptot,f] = modelProteinContent(model)

    % Create and constrain the total enzymes pool
    cd ([GECKO_path '/Matlab_Module/limit_proteins'])
    sigma               = 0.5;      %Average enzyme saturation factor
    HepG2_ecModel_batch = constrainEnzymes(HepG2_ecModel,Ptot,sigma);
    cd (current)
    cd ../HMR2.0/EC_HMR
    save('ecHepG2_batch.mat','HepG2_ecModel','HepG2_model_data',...
                                    'HepG2_kcats','HepG2_ecModel_batch')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%