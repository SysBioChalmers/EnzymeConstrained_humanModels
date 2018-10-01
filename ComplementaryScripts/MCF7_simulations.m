%Load MCF7 ecModel
cd ../models/MCF-7/ecMCF7
load('MCF7_ecModel_batch.mat')
model = MCF7_ecModel_batch;
cd ../../../ComplementaryScripts
%Set EMEM medium
[model,essential,sensitivities] = setEMEMmedium(model,0.0328);
%Set suboptimal growth rate
gRate = 0.0273;
metList = {'acetyl-CoA','malonyl-CoA','CoA'};
[resultsProduction_OE,~] = metEng_TargetsFinder(model,'HMR_9034',metList,'production','cytosol','deletion',gRate);
metList = [metList,'growth'];
metList = strrep(metList,'-','_');
rows    = [{'Yields'};model.genes];
%Write table file in the results folder
resultsProduction_OE = num2cell(resultsProduction_OE);
T = cell2table(resultsProduction_OE,'VariableNames',metList,'RowNames',rows);
cd ../Results
writetable(T,'resultsProduction_OE.txt')
%Repeat again but using the consumption strategy
cd ../ComplementaryScripts
metList = {'pyruvate','acetate','citrate','acetyl-CoA','malonyl-CoA','CoA'};
[resultsConsumption_OE,~] = metEng_TargetsFinder(model,'HMR_9034',metList,'consumption','cytosol','deletion',gRate);
metList = [metList,'growth'];
metList = strrep(metList,'-','_');
rows    = [{'Yields'};model.genes];
resultsConsumption_OE = num2cell(resultsConsumption_OE);
T = cell2table(resultsConsumption_OE,'VariableName',metList);
cd ../Results
writetable(T,'resultsConsumption_OE.txt')