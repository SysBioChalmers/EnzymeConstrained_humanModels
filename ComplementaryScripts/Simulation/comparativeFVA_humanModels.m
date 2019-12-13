function results = comparativeFVA_humanModels(cellLine)
%comparativeFVA_humanModels
%
% Function that runs a comparative flux variabiability analysis between a
% GEM and its enzyme constrained counterpart for a given cell-line specific
% model.
%
% cellLine  (string) cell-line name (It should be consistent with the name
%           of the subfolder in which the model is stored.
%
% results   (table) Contains rxnIDs, rxn formulas, Flux variability 
%           distributions for both GEM and ecGEM in [mmol/gDw h] and
%           metabolic subSystems for all rxns for which a variability range
%           in bothmodels was calculated
%
% Usage: FVA_Dists = comparativeFVA_humanModels(cellLine)
%
%   Ivan Domenzain, 2019-12-05

current    = pwd;
%Clone GECKO and pull comparativeFVa branch
system('git clone https://github.com/SysBioChalmers/GECKO.git');
cd GECKO
GECKO_path = pwd;
system('git checkout 388d9ad');
%Load GEM and ecGEM
load(['../../../models/' cellLine '/' cellLine '.mat'])
load(['../../../models/' cellLine '/ecModel_batch.mat'])
eval(['model = ' cellLine ';'])
%Set medium constraints
cd (current)
model   = removeMetFields(model);
model   = setHamsMedium(model,false);
ecModel = setHamsMedium(ecModel_batch,true);
evalin( 'base', 'clear(''model_modified'')' )
evalin( 'base', 'clear(''ecModel_batch'')' )
%Use GECKO built-in function for FVA comparative analysis
cd ([GECKO_path '/geckomat/utilities/FVA'])
CsourceUptk = 'HMR_9034';
[FVA_Dists,indexes,~,~] = comparativeFVA(model,ecModel,CsourceUptk,false,1E-8);
cd (current)
%Write results in a table and save it as .txt file
mkdir('../../Results/FVA')
variables = {'rxns' 'formulas' 'model_ranges' 'ecModel_ranges' 'subSystems'};
formulas  = constructEquations(model,indexes);
results   = table(model.rxns(indexes),formulas,FVA_Dists{1},FVA_Dists{2},model.subSystems(indexes),'VariableNames',variables);
writetable(results,['../../Results/FVA/FVA_comp_' cellLine],'Delimiter','\t','QuoteStrings',false)
end
%--------------------------------------------------------------------------
function model = removeMetFields(model)
if isfield(model,'inchis')
    model = rmfield(model,'inchis');
end
if isfield(model,'metMiriams')
    model = rmfield(model,'metMiriams');
end
if isfield(model,'metCharges')
    model = rmfield(model,'metCharges');
end
if isfield(model,'metFrom')
    model = rmfield(model,'metFrom');
end
end