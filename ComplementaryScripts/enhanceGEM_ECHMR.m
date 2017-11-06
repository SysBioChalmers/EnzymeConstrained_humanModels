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
% Ivan Domenzain.      Last edited: 2017-10-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecModel,model_data,kcats] = enhanceGEM_ECHMR

    org_name     = 'homo sapiens';
    keggCode     = 'hsa';
    GECKO_path   = '/Users/ivand/Desktop/GECKO/GECKO-master';
    KEGG_path    = '/Volumes/ftp.bioinformatics.jp/kegg';
    EC_HMR_path  = '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM';
    Protdatabase = 'hsa_ProtDatabase.mat';
    toolbox      = 'COBRA';
    format short e    
    % COBRA toolbox required
    initCobraToolbox;
    cd ../HMR2.0
    hsa_model = importModel('HMRdatabase2_00_Cobra.xml',true,true);
    % Update protein databases (KEGG and uniprot)
    cd ../ComplementaryScripts
    DB_path =  '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM/Databases';
    updateDatabases(GECKO_path,DB_path,keggCode);
    
    % Creates a new field in the model structure where the ENSEMBL gene IDs
    % are converted to their short gene names in order to provide
    % compatibility with the kcat matching algorithms
    cd ([EC_HMR_path '/ComplementaryScripts'])
    hsa_model = substituteEnsemblGeneIDs(hsa_model);     
    % Standardizes the metabolites names
    hsa_model = modifyMetNames(hsa_model);
    cd ([EC_HMR_path '/HMR2.0'])
    save('human.mat','hsa_model');
    
    % Retrieve kcats & MWs for each rxn in model:
    cd ([EC_HMR_path '/ComplementaryScripts'])
    hsa_model_data = getEnzymeCodes_ECHMR(hsa_model,Protdatabase,GECKO_path);
    cd ([EC_HMR_path '/HMR2.0'])
    save('human_enzData.mat','hsa_model','hsa_model_data');
    cd ../ComplementaryScripts
    hsa_kcats      = matchKcats(hsa_model_data,org_name,GECKO_path,KEGG_path);
    cd ([EC_HMR_path '/HMR2.0'])
    save('human_kcats.mat','hsa_kcats');
     
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    hsa_ecModel = readKcatData(hsa_model_data,hsa_kcats);
    % Save output models:
    cd ([EC_HMR_path '/HMR2.0/EC_HMR'])
    save('EC_HMR.mat','hsa_ecModel','hsa_model_data','hsa_kcats')
    cd '../../xml'
    exportModel(hsa_ecModel,'hsa_ecModel.xml',true)
    pathTXT = [EC_HMR_path '/txt'];
    exportToTabDelimited(hsa_ecModel,pathTXT);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%