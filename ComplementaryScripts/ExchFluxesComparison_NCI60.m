function [meanGRerror,minGRerror, maxGRerror] = ExchFluxesComparison_NCI60(constLevel,ecFlag,fixedBounds)
%ExchFluxesComparison_NCI60
% 
% Function that performs pFBA simulations on the 11 cell-line specific
% models (from the NCI60 cell-lines) and computes errors in exchange fluxes
% and growth rate predictions compared to experimental data.
%
%   constLevel      (integer) 0 for EMEM constraints, 1 for level0+addition
%                   of GUR constraint, 2 for Level1+threonin, 3 for Level2+lysine, 
%                   4 for Level 3+glutamine and 5 for Level 4+serine
%                   constraint.
%   ecFlag          True if ecModels are desired for simulation or false
%                   for normal models.
%   fixedBounds     True if bounds are desired to be fixed (lb and ub).
% 
%   meanGRerror     Mean growth rate prediction relative error across the
%                   11 cell-line models.
%   minGRerror      Min growth rate prediction error
%   maxGRerror      Max growth rate prediction error
%
% Usage: [meanGRerror,minGRerror, maxGRerror] = ExchFluxesComparison_NCI60(constLevel,ecFlag,fixedBounds)
%
% Ivan Domenzain.      Last edited: 2019-05-15

clc
close all
current      = pwd;
cd ../models/humanGEM_cellLines
modelsFolder = pwd; 
%Assign specific color code for the different cell-lines
[folders, colorS] = assignNamesAndColors;
%Open exchange fluxes data file
fID     = fopen('../../DataFiles/exchangeFluxes_NCI60.txt');
expData = textscan(fID,'%s %s %f %f %f %f %f %f %f %f %f %f %f',...
                   'Delimiter','\t','HeaderLines',1);
%%Save experimental data 
expIDs = expData{1};expMets = expData{2};expData = expData(3:end);
%Initialize variables
legendStr       = {};
predictedFluxes = zeros(length(expIDs)-1,length(folders));
errors_ExFlux   = zeros(length(folders),3);
errors_GRates   = zeros(length(folders),1);
%For each cell line
for i=1:length(folders)
    %Measured fluxes for the i-th cell-line
    measuredFluxes = expData{i};
    cellLineStr    = strrep(folders{i},'_','-');
    cd (modelsFolder)
    cd (folders{i})
    mkdir Results
    disp(folders{i})
    %Load ecModel and add HepG2 biomass reaction
    load ecModel_batch.mat
    load model.mat
    model.b = zeros(length(model.mets),1);
    cd ../../../ComplementaryScripts
    ecModel_batch  = substituteBiomassRxns(ecModel_batch);
    model_modified = substituteBiomassRxns(model);
    %Allow any level of average enzyme saturaion (0 - 100%)
    P_index = find(strcmpi(ecModel_batch.rxnNames,'prot_pool_exchange'));
    ecModel_batch.ub(P_index) = 0.609*0.5;
    %Save models
    cd ([modelsFolder, '/', folders{i}])
    save ('ecModel_batch.mat', 'ecModel_batch')
    save ('model_modified.mat', 'model_modified')
    %Run FBA 
    ecMsol = solveLP(ecModel_batch);
    GEMsol = solveLP(model_modified);
    %If both models are feasible
    if ~isempty(ecMsol.x) & ~isempty(GEMsol.x)
        cd ../../../ComplementaryScripts/Simulation
        %Set EMEM constraints
        [ecModel_batch,~]  = setEMEMmedium(ecModel_batch,1000,1000,1000,true);
        [model_modified,~] = setEMEMmedium(model_modified,1000,1000,1000,false);
        %Set fixed constraints (nutrients uptakes) 
        ecModel_batch      = setDataConstraints(ecModel_batch,expData{i},expIDs,true,constLevel,fixedBounds);
        model_modified     = setDataConstraints(model_modified,expData{i},expIDs,false,constLevel,fixedBounds);
        %Get exchange fluxes and metabolites IDs in the original model
        model_modified       = rmfield(model_modified,'unconstrained');
        [EX_IDs,exc_Indexes] = getExchangeRxns(model_modified);
        exchangeMets = {};
        for j=1:length(expIDs)-1
            metJ         = find(strcmpi(model_modified.rxns,expIDs(j)));
            exchangeMets = [exchangeMets;...
                model_modified.metNames(find(model_modified.S(:,metJ),1))];
        end
        %Get parsimonious solutions for both models
        ecMsol = minProtSimulation(ecModel_batch);
        GEMsol = solveLP(model_modified,1);
        %If both models are feasible then save results
        if ~isempty(ecMsol.x) & ~isempty(GEMsol.x)
            %Extract values for exchange fluxes from GEMsol
            exc_Model = GEMsol.x(exc_Indexes);
            cd ..
            exc_ecModel =[];
            %Map each exchange rxn to ecModel_batch
            for j=1:length(exc_Indexes)
                indexes = [];
                indexes = rxnMapping(EX_IDs{j},ecModel_batch,true);
                if length(indexes)>1
                    exc_ecModel = [exc_ecModel; ecMsol.x(indexes(1))-ecMsol.x(indexes(2))];
                else
                    exc_ecModel = [exc_ecModel; ecMsol.x(indexes)];
                end
            end
            %Keep just the exchange fluxes that are also part of the
            %dataset and compare both GEM and ecGEM predictions
            [~,toKeep,order]  = intersect(EX_IDs,expIDs);
            str               = strrep(folders{i},'_','-');
            if length(toKeep) == length(expIDs(1:end-1))
                if ecFlag
                    predictions = exc_ecModel(toKeep);
                    fileName1 =['ecGEM_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
                    fileName2 = ['ecGEM_const_' num2str(constLevel) '_errorMetrics.txt'];
                    fileName3 = ['ecGEM_const_' num2str(constLevel) '_error_GRate.txt'];
                else
                    predictions = exc_Model(toKeep);
                    fileName1 = ['GEM_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
                    fileName2 = ['GEM_const_' num2str(constLevel) '_errorMetrics.txt'];
                    fileName3 = ['GEM_const_' num2str(constLevel) '_error_GRate.txt'];
                end
                predictedFluxes(:,i) = predictions;
                experimental         = measuredFluxes(order);
                exchangeMets         = exchangeMets(order);
                exchangeIDs          = expIDs(order);
                [Xvalues,Yvalues,direction,errors_ExFlux(i,:)] = getPlotValues(experimental,predictions);
                %Calculate mean absolute error for Growth rate predictions
                errors_GRates(i)     = computeErrorMetric(predictions(end),experimental(end),'MRE');
                disp(['Relative error for GRate prediction: ' num2str(errors_GRates(i))])
                plotDataPoints(Xvalues,Yvalues,colorS(i,:),direction)
                legendStr = [legendStr; [cellLineStr '/ RSME = ' num2str(errors_ExFlux(i,3))]];
                hold on
                
            else
                disp('Inconsistent mapping')
            end
        else
            disp(['Not feasible 2: ' folders{i}])
        end
    else
        disp(['Not feasible 1: ' folders{i}])
    end
    %write file with experimental and predicted exchange fluxes
    variables = {'Rxn_ID' 'Metabolite' 'experimental' 'predictions'};
    T = table(exchangeIDs,exchangeMets,experimental,predictions,'VariableNames',variables);
    writetable(T,[modelsFolder, '/', folders{i},'/Results/',fileName1],'QuoteStrings',false,'Delimiter','\t')
end
%Plot exchange fluxes prediction comparison
x1 = linspace(-9,1,100);
plot(x1,x1)
legend(legendStr)
hold on
%Write output files
cd (modelsFolder)
cd ../../Results
mkdir 11_cellLines_NCI60
cd 11_cellLines_NCI60
%write file with different error metrics for all cell lines
variables = {'cell_Line' 'pearson' 'MAE' 'RMSE'};
T         = table(folders',errors_ExFlux(:,1),errors_ExFlux(:,2),errors_ExFlux(:,3),'VariableNames',variables);
writetable(T,fileName2,'QuoteStrings',false,'Delimiter','\t')
%write file with GRate predictions MAE
variables = {'cell_Line' 'MAE'};
T         = table(folders',errors_GRates,'VariableNames',variables);
writetable(T,fileName3,'QuoteStrings',false,'Delimiter','\t')
meanGRerror = mean(errors_GRates);
minGRerror  = min(errors_GRates);
maxGRerror  = max(errors_GRates);
cd (current)
end
%--------------------------------------------------------------------------
function [cellNames, colors] = assignNamesAndColors
cellNames = [{'HS_578T'} {'RPMI_8226'} {'HT29'} {'MALME_3M'} {'SR'}...
           {'UO_31'} {'MDAMB_231'} {'HOP62'} {'NCI_H226'} {'HOP92'}...
           {'O_786'}];

colors    = [0    0    128             %dark blue
             0    0    255             %blue
             0    170  255             %light blue
             0    128  0               %forest green
             128  255  0               %lime green
             255  105  180             %pink
             255  200  050             %yellow
             191  0    255             %purple
             255  128  0               %orange
             0    0    0               %black
             255  0    0]./255;        %red
end     
%--------------------------------------------------------------------------
function [Xvalues,Yvalues,direction,errors] = getPlotValues(predictions,experimental)
Xvalues      = (log10(abs(experimental)+1E-9));
Yvalues      = (log10(abs(predictions)+1E-9));
direction    = sign(experimental);
MRE          = computeErrorMetric(experimental,predictions,'MRE');
RMSE         = computeErrorMetric(experimental,predictions,'RMSE');
pearson      = computeErrorMetric(experimental,predictions,'pearson');
errors       = [pearson,MRE,RMSE];
end
%--------------------------------------------------------------------------
function plotDataPoints(Xvalues,Yvalues,color,direction)
sizes   = [];
colorS  = zeros(length(Xvalues),3);
for i=1:length(Xvalues)
    if direction(i)<=0
        sizes = [sizes; 75];
    else
        sizes = [sizes; 150];
    end
    colorS(i,:)  = color;
end
scatter(Xvalues,Yvalues,sizes,colorS,'d','filled')
end
%--------------------------------------------------------------------------
function avg_Eps = computeErrorMetric(data,predictions,method)
data(data==0) = 1E-9;
switch method
    case 'pearson'
        avg_Eps  = corrcoef(predictions,data);
        avg_Eps  = avg_Eps(1,2);
    case 'MRE'
        relErr   = (data-predictions)./data;
        relErr   = (abs(relErr));
        avg_Eps  = mean(relErr);
    case 'RMSE'
        total = 0;
        T     = length(data);
        for i=1:T
            total = total + (data(i) - predictions(i))^2;
        end
        avg_Eps = sqrt(total/(T));
end
end
%--------------------------------------------------------------------------
function  solution = minProtSimulation(model)
sol      = solveLP(model,1);
objIndex = find(model.c);
objValue = sol.x(objIndex);
%Fix objective value
model.ub(objIndex) = 1.001*objValue;
model.lb(objIndex) = 0.999*objValue;
%Find protein pool exchange pseudoreaction and set it as an objetive to
%minimize
model.c           = zeros(length(model.c),1);
protIndx          = find(strcmpi(model.rxns,'prot_pool_exchange'));
model.c(protIndx) = -1;
solution          = solveLP(model,1);
end
%--------------------------------------------------------------------------
function model = setDataConstraints(model,fluxes,expIDs,ecFlag,const_level,fixed)
%Biomass UB
GrowthRate = fluxes(end-1);
%value   = GrowthRate;
GrowthIndx = find(strcmpi(model.rxns,expIDs(end-1)));
%model   = setBounds(model,GrowthIndx,value,ecFlag,true);
if ecFlag
    Prot_biomass = fluxes(end);         %Total protein content in biomass [g prot/g DW]
    Prot_pool    = fluxes(end)*0.5*0.75; %Amount of enzymes available for biochemical reactions
    model        = rescaleBiomassProtein(model,Prot_pool,ecFlag);
end
%Glucose exchange bound
if const_level >0
    %value   = GrowthRate;
    %model   = setBounds(model,GrowthIndx,value,ecFlag,true);
    index   = strcmpi(expIDs,'HMR_9034');
    rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
    value   = fluxes(index);
    model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
    %Threonine exchange bound
    if const_level>1
        index   = strcmpi(expIDs,'HMR_9044');
        rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
        value   = fluxes(index);
        model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
        %Lysine exchange bound
        if const_level>2 
            index   = strcmpi(expIDs,'HMR_9041');
            rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
            value   = fluxes(index);
            model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
            %Glutamine exchange bound
            if const_level>3
                index   = strcmpi(expIDs,'HMR_9063');
                rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
                value   = fluxes(index);
                model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
                %Serine exchange bound
                if const_level>4
                    index   = strcmpi(expIDs,'HMR_9069');
                    rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
                    value   = fluxes(index);
                    model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
                end
            end
        end
    end
end
end
%--------------------------------------------------------------------------
function model = setBounds(model,index,value,ecFlag,fixed)
if fixed
    lowValue = 0.9;
else
    lowValue = 0;
end
direction = sign(value);
if direction <0
    if ecFlag
        %Block the opposite direction (production)
        model.ub(index) = 0;
        %Find uptake reaction
        rxn   = [model.rxns{index} '_REV'];
        index = find(strcmpi(model.rxns,rxn));
        model.ub(index) = abs(value);
        model.lb(index) = lowValue*abs(value);
    else
        model.lb(index) = value;
        model.ub(index) = lowValue*value;
    end
else
    model.ub(index) = value;
    model.lb(index) = lowValue*value;
end
end
%--------------------------------------------------------------------------
function model = rescaleBiomassProtein(model,Ptot,constrainPool)
if nargin <3 
    gIndex    = find(strcmpi(model.rxns,'HumanGrowth'));
    protIndex = find(strcmpi(model.metNames,'proteinPool'));
    coeff     = model.S(protIndex,gIndex); %Extract coeff corresponding to 0.609 g Prot / g Biomass
    newCoeff  = coeff*Ptot/0.609;
    %Rescale stoichiometric coefficient of proteinPool in biomass equation
    model.S(protIndex,gIndex) = newCoeff;
else
    if constrainPool
        protPool           = find(strcmpi(model.rxns,'prot_pool_exchange'));
        model.ub(protPool) = Ptot;
    end
end
end
