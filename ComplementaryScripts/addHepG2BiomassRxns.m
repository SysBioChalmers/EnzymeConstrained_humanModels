function HepG2model = addHepG2BiomassRxns(model,ecFlag)
%ADDHEPG2BIOMASSRXNS  Add HepG2 biomass reactions to (ec)HMR2 model.
%
% HEPG2MODEL = ADDHEPG2BIOMASSRXNS(MODEL) adds reactions associated with
% biomass production to the HMR2 genome scale metabolic model. If the model
% is enzyme constrained (ecHMR), then it inserts the new reactions and
% metabolites into the model at indices just after the original HMR2
% content, and before the added EC content.


%% Remove boundary metabolites, if present
% boundary metabolites (those in compartment "x" for RAVEN-type models)
% should be removed from the stoichiometry matrix before running any
% simulations, as they prevent mass balancing. Currently, this just sets
% their stoichiometric coefficients to zero, but probably should be revised
% to remove them from model entirely.

% just in case, confirm that this is not a Cobra-style model
if all(ismember(model.comps,{'c','e','g','l','m','n','r','x','b'}))
    warning('Model appears to be Cobra-like: boundary metabolites will not be removed.');
else
    % set S-matrix to zero for all mets in "x" boundary compartment
    [~,boundCompInd] = ismember('x',model.comps);
    boundMets = (model.metComps == boundCompInd);
    model.S(boundMets,:) = 0;
end


%% Add new metabolites for biomass reactions
% metabolite ID numbering started arbitrarily from m90000
addMetData = {
'm90000c',	'proteinPool',              'c'
'm90001c',	'growthMaintenance',        'c'
'm90002c',	'biomass',                  'c'
'm90002c',	'biomass',                  's'
'm90003c',	'human_RNAPool',            'c'
'm90004c',	'human_DNAPool',            'c'
'm90005c',  'phosphatidate',            'c'
'm90006c',  'CDP-diacylglycerol',       'c'
'm90007c',  'phosphatidylserine',       'c'
'm90008c',  'phosphatidylcholine',      'c'
'm90009c',  'phosphatidylethanolamine',	'c'
'm90010c',  'phosphatidylPool',         'c'
'm90011c',  'fattyAcidPool',            'c'
'm90012c',  'lipidPool',                'c'};

% insert new metabolites into model structure where the original (non-EC) metabolites ended

% search for max index of model.met of the format "m###..."
insMetInd = find(~cellfun(@isempty,regexp(model.mets,'^m\d+')),1,'last');  

% insert info into model
model.mets = [model.mets(1:insMetInd); addMetData(:,1); model.mets(insMetInd+1:end)];
model.b = [model.b(1:insMetInd); zeros(size(addMetData,1),1); model.b(insMetInd+1:end)];
model.metNames = [model.metNames(1:insMetInd); addMetData(:,2); model.metNames(insMetInd+1:end)];
model.S = [model.S(1:insMetInd,:); zeros(size(addMetData,1),size(model.S,2)); model.S(insMetInd+1:end,:)];
model.metFormulas = [model.metFormulas(1:insMetInd); repmat({''},size(addMetData,1),1); model.metFormulas(insMetInd+1:end)];

% determine metabolite compartment numbers and add to model
[~,addMetComps] = ismember(addMetData(:,3),model.comps);
model.metComps = [model.metComps(1:insMetInd); addMetComps; model.metComps(insMetInd+1:end)];



%% Addition of biomass component reactions
addRxnData = {
'human_ATPMaintenance',     'ATP[c] + H2O[c] => ADP[c] + Pi[c]'
'human_proteinPool',        '0.0937 alanine[c] + 0.0507 arginine[c] + 0.04 asparagine[c] + 0.04 aspartate[c] + 0.0142 cysteine[c] + 0.1118 glutamine[c] + 0.1831 glycine[c] + 0.0198 histidine[c] + 0.0309 isoleucine[c] + 0.0664 leucine[c] + 0.0571 lysine[c] + 0.0156 methionine[c] + 0.029 phenylalanine[c] + 0.0853 proline[c] + 0.0491 serine[c] + 0.0402 threonine[c] + 0.019 tyrosine[c] + 0.0471 valine[c] + 0.0072 tryptophan[c] => proteinPool[c]'
'HumanGrowth',              'glycogen[c] + 0.3 lipidPool[c] + 6 proteinPool[c] + 80 growthMaintenance[c] + 0.10 human_DNAPool[c] + 0.10 human_RNAPool[c] => biomass[c]'
'human_GrowthMaintenance',  'ATP[c] + H2O[c] => ADP[c] + Pi[c] + growthMaintenance[c]'
'humanGrowthOut',           'biomass[s] => '
'humanGrowthTransport',     'biomass[c] => biomass[s]'
'human_DNAPool',            '0.3 dAMP[c] + 0.2 dCMP[c] + 0.2 dGMP[c] + 0.3 dTMP[c] => human_DNAPool[c]'
'human_RNAPool',            '0.18 AMP[c] + 0.30 CMP[c] + 0.34 GMP[c] + 0.18 UMP[c] => human_RNAPool[c]'
'ApproxPhosphatidate',      'ATP[c] + 2 palmitate[c] + sn-glycerol-3-phosphate[c] => AMP[c] + PPi[c] + phosphatidate[c]'
'ApproxCDPdiacylglycerol',  'CTP[c] + phosphatidate[c] => PPi[c] + CDP-diacylglycerol[c]'
'ApproxPSerine',            'serine[c] + CDP-diacylglycerol[c] => CMP[c] + phosphatidylserine[c]'
'ApproxPCholine',           'choline[c] + CDP-diacylglycerol[c] => CMP[c] + phosphatidylcholine[c]'
'ApproxPEthanolAmine',      'ethanolamine[c] + CDP-diacylglycerol[c] => CMP[c] + phosphatidylethanolamine[c]'
'PhosphatidylPool',         '0.189 phosphatidylserine[c] + 0.223 phosphatidylcholine[c] + 0.588 phosphatidylethanolamine[c] => phosphatidylPool[c]'
'FattyAcidPool',            'H2O[c] + 0.003 margaric acid[c] + 0.126 myristic acid[c] + 1.056 oleate[c] + 1.308 palmitate[c] + 0.222 palmitolate[c] + 0.012 pentadecylic acid[c] + sn-glycerol-3-phosphate[c] + 0.27 stearate[c] => Pi[c] + fattyAcidPool[c]'
'lipidPool',                '0.331 cholesterol[c] + 0.384 phosphatidylPool[c] + 0 fattyAcidPool[c] => lipidPool[c]'};

nRxns = size(addRxnData,1);  % number of new reactions to be added

% organize new reaction data into a structure
addRxnDataStruct.rxns = addRxnData(:,1);
addRxnDataStruct.rxnNames = addRxnData(:,1);
addRxnDataStruct.equations = addRxnData(:,2);
addRxnDataStruct.lb = zeros(nRxns,1);
addRxnDataStruct.ub = Inf(nRxns,1);
addRxnDataStruct.subSystems = repmat({'Artificial'},nRxns,1);

% generate new temporary model with reaction info added
temp_model = addRxns(model,addRxnDataStruct,3,[],false);

% Remove the previous biomass rxn
massPos = [];
for rxnText = {'biomass','cofactors_vitamins','vitaminA','vitaminD','vitaminE'}
    massPos = [massPos; find(~cellfun(@isempty,strfind(temp_model.rxns,rxnText)))];
end

if ~isempty(massPos)
    temp_model = removeReactions(temp_model,massPos);
end

% The addRxns function adds new reaction data to the end of each model
% field. However, we want to insert the reaction data right before the new
% reactions that were added by the GECKO algorithm.
if ecFlag
    % search in model.rxns for first appearance of "..._REV" reaction name
    insRxnInd = find(~cellfun(@isempty,regexp(model.rxns,'_REV$')),1,'first')-1;

    % transfer new reaction fields from TEMP_MODEL to proper location in MODEL
    %rxnFields = {'rxns';'rxnNames';'grRules';'subSystems';'rules';'S';'lb';'ub';'rev';'c';'rxnGeneMat'};
    rxnFields = {'rxns';'rxnNames';'grRules';'subSystems';'S';'lb';'ub';'rev';'c';'rxnGeneMat'};

    for i = 1:length(rxnFields)
        if strcmp(rxnFields{i},'S')
            % S matrix has rxns as columns instead of rows, so it needs to be
            % processed differently than the other model fields
            model.S = [model.S(:,1:insRxnInd), temp_model.S(:,end-(nRxns-1):end), model.S(:,insRxnInd+1:end)];
        else
            model.(rxnFields{i}) = [model.(rxnFields{i})(1:insRxnInd,:); ...
                temp_model.(rxnFields{i})(end-(nRxns-1):end,:); model.(rxnFields{i})(insRxnInd+1:end,:)];
        end
    end
else
    model = temp_model;
end


%% Avlant's stoichiometry changes/notes regarding HepG2 biomass

% % I run the following lines to update the biomass equation for my
% % current project, the most important part is that it sets the lipidpool,
% % to 0, everything else should be qualitative.
% model = configureSMatrix(model, 80, 'HumanGrowth', 'human_growthMaintenance[c]');
% model = configureSMatrix(model, 6, 'HumanGrowth', 'human_protein_pool[c]');
% model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_RNAPool[c]');
% model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_DNAPool[c]');
% model = configureSMatrix(model, 0, 'lipidPool', 'fattyAcidPool[c]');
% model = configureSMatrix(model, 0.3, 'HumanGrowth', 'lipidPool[c]');
% 
% % Maybe consider setting glycogen to 0. I know Gatto once predicted that
% % knocking out glycogen production would be lethal to cells. But its not
% % obvious that they need it.
% model = configureSMatrix(model, 1, 'HumanGrowth', 'glycogen[c]');

%% Finalize changes

% set objective as 'HumanGrowth' reaction
model.c(:) = 0;
model.c(ismember(model.rxns,'humanGrowthOut')) = 1;

% assign output
HepG2model = model;






