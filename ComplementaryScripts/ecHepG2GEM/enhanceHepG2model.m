function HepG2_ecModel = enhanceHepG2model
% enhanceHepG2model
%   Converts the HepG2 specific GEM to an Enzyme constrained version
%
%
%   Ivan Domenzain, 2018-04-05
%
    current    = pwd;
    org_name   = 'homo sapiens';
    GECKO_path = '/Users/ivand/Documents/GitHub/GECKO';
    cd ../models/HepG2
    %Import HepG2 model xml file and save it as a matlab structure
    HepG2model = importModel('Hep-G2.xml');
    save('HepG2model.mat','HepG2model')
    cd ..
    [grRules,rxnGeneMat]  = standardizeGrRules(HepG2model);
    HepG2model.grRules    = grRules;
    HepG2model.rxnGeneMat = rxnGeneMat;
    HepG2model_modified   = substitute_grRules(HepG2model); 
    HepG2model_modified   = substituteEnsemblGeneIDs(HepG2model_modified); 
    % Standardizes the metabolites names
    cd (current)
    HepG2model_modified   = modifyMetNames(HepG2model_modified);
    cd (current)
    cd ../models/HepG2
    save('HepG2model_modified.mat','HepG2model_modified')
    %%%%%%%%%%%%%%%%%%%%%%%% GECKO modifications %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve kcats & MWs for each rxn in the model:
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    HepG2_model_data = getEnzymeCodes(HepG2model_modified);
    %Tries to Match kinetic coefficients to every reaction with a non empty
    % grRule
    HepG2_kcats      =  matchKcats(HepG2_model_data,org_name);
    cd (current)
    cd ../models/HepG2/Data
    save('HepG2_enzData.mat','HepG2_model_data','HepG2_kcats');
    % Integrate enzymes in the model:
    cd ([GECKO_path '/Matlab_Module/change_model'])
    HepG2_ecModel     = readKcatData(HepG2_model_data,HepG2_kcats);
    [HepG2_ecModel,~] = manualModificationsGeneral(HepG2_ecModel);
    % Save output models:
    cd (current)
	cd ../models/HepG2/ecHepG2
	save('HepG2_ecModel.mat','HepG2_ecModel')
end