%function FVA_Dists = comparativeFVA_humanModels(cellLine)
%comparativeFVA_humanModels
%
% Function that runs a comparative flux variabiability analysis between a
% GEM and its enzyme constrained counterpart for a given cell-line specific
% model.
%
% cellLine  (string) cell-line name (It should be consistent with the name
%           of the subfolder in which the model is stored.
%
% FVA_Dists (2x1 cell) Flux variability distributions for both GEM and
%           ecGEM in [mmol/gDw h]
%
% Usage: FVA_Dists = comparativeFVA_humanModels(cellLine)
%
%   Ivan Domenzain, 2019-05-16

current    = pwd;
%Clone GECKO and pull comparativeFVa branch
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
GECKO_path = pwd;
git ('checkout fix/comparativeFVA')
git pull
%Load GEM and ecGEM
load(['../../../models/humanGEM_cellLines/' cellLine '/model_modified.mat'])
load(['../../../models/humanGEM_cellLines/' cellLine '/ecModel_batch.mat'])
%Set medium constraints
glucBound   = 1;
oxBound     = 1000;
stdBound    = 1;
cd (current)
model       = removeMetFields(model_modified);
[model,~]   = setEMEMmedium(model,glucBound,oxBound,stdBound,false);
[ecModel,~] = setEMEMmedium(ecModel_batch,glucBound,oxBound,stdBound);
evalin( 'base', 'clear(''model_modified'')' )
evalin( 'base', 'clear(''ecModel_batch'')' )
%Use GECKO built-in function for FVA comparative analysis
cd ([GECKO_path '/geckomat/utilities/FVA'])
CsourceUptk       = 'HMR_9034';
[FVA_Dists,~,~,~] = comparativeFVA(model,ecModel,CsourceUptk,false);
cd (current)
cd (['../../models/humanGEM_cellLines/' cellLine])
save('FVA_results.mat','FVA_Dists','indexes')
%end
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