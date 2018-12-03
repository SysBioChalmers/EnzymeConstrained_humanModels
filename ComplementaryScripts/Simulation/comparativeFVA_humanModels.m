function [FVA_Dists,indexes,stats] = comparativeFVA_humanModels(cellLine)
current = pwd;
GECKO_path = '/Users/ivand/Documents/GitHub/GECKO';
load(['../../models/' cellLine '/model_modified.mat'])
load(['../../models/' cellLine '/ecModel_batch.mat'])
glucBound = 1;
oxBound   = 20;
stdBound  = 20;
[model,~] = setEMEMmedium(model_modified,glucBound,oxBound,stdBound,false);
[ecModel,]= setEMEMmedium(ecModel_batch,glucBound,oxBound,stdBound);
evalin( 'base', 'clear(''model_modified'')' )
evalin( 'base', 'clear(''ecModel_batch'')' )
%Use GECKO built-in function for FVA comparative analysis
cd ([GECKO_path '/geckomat/utilities/FVA'])
[FVA_Dists,indexes,stats] = comparativeFVA(model,ecModel,'HMR_9034',false,0);
cd (current)
cd (['../../models/' cellLine '/Results'])
save('FVA_results.mat','FVA_Dists','indexes')
end