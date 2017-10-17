%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
% 
%  Benjam?n J. S?nchez. Last edited: 2017-04-12
%  Ivan Domenzain.      Last edited: 2017-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%function [ecModel,model_data,kcats] = enhanceGEM_ECHMR(model,toolbox,name)
    GECKO_path   = '/Users/ivand/Desktop/GECKO-IVAN';
    EC_HMR_path  = '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM';
    org_name     = 'homo sapiens';
    Protdatabase = 'hsa_ProtDatabase.mat';
    toolbox = 'COBRA';
    format short e
    if strcmp(toolbox,'COBRA')
       initCobraToolbox;
    end
    % 
    cd ../HMR2.0
    hsa_model = importModel('HMRdatabase2_00_Cobra.xml',true,true)
    save('human.mat','hsa_model');

    %Add some RAVEN fields for easier visualization later on:
    %hsa_model      = standardizeModel(sce_model,toolbox);

    cd ../ComplementaryScripts
    %updateDatabases(GECKO_path)
    cd ([EC_HMR_path '/ComplementaryScripts'])
    [hsa_model] = substituteEnsemblGeneIDs(hsa_model);
    %Retrieve kcats & MWs for each rxn in model:
    cd ([EC_HMR_path '/ComplementaryScripts'])
    hsa_model_data = getEnzymeCodes_ECHMR(hsa_model,Protdatabase,GECKO_path);
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    hsa_kcats      = matchKcats(hsa_model_data);
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    hsa_ecModel = readKcatData(hsa_model_data,hsa_kcats);
   
    %Save output models:
    cd ([EC_HMR_path '/HMR2.0/EC_HMR'])
    save(['EC_HMR.mat'],'hsa_ecModel','hsa_model_data','hsa_kcats')
    %saveECmodelSBML(hsa_ecModel,'hsa_ecModel');
    cd ([GECKO_path '/Matlab_Module'])

%end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%