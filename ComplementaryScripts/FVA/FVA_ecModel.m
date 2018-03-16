%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rangeDist,rangeDist_EC] = FVA_ecirrevModel(irrevModel,ecirrevModel,BlockFlag)
%  
% This function goes through each of the rxns in a metabolic irrevModel and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
% Raphael Ferreira     Last edited: 2018-03-14
% Ivan Domenzain.      Last edited: 2018-03-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [rangeDist,rangeDist_EC] = FVA_ecirrevModel(irrevModel,ecirrevModel,BlockFlag,path)
    path = '/Users/ivand/Documents/GitHub/GECKO';
    current          = pwd;
    range_model      = [];
    optimizedIndxs   = [];
    optimizedECIndxs = [];
    range_ecModel    = [];   
    %Get the index for all the non-objective rxns in the original irrevModel
    rxnsIndxs    = find(model.c~=1); 
    %Convert to irreversible irrevModel
    cd ([path '/Matlab_Module/change_model'])
    irrevModel = convertToIrreversibleModel(model);
    cd (current)
    %Block glucose and oxygen production
	%if nargin>2
		Blocked_mets = {'glucose', 'oxygen'};
		irrevModel   = block_production(irrevModel,Blocked_mets,true);
        ecModel      = block_production(ecModel,Blocked_mets,true);
    %end    
     %Gets the optimal value for ecirrevModel and fixes the objective value to 
     %this for both irrevModels
     [OptimalValue,basalFluxDist,ecModel] = fixObjective(ecModel);
     [OptimalValue,~,irrevModel]          = fixObjective(irrevModel,OptimalValue);
     % Get the variability range for each of non-objective reactions in the
     % original irrevModel
     for i=1:length(rxnsIndxs)
         indx        = rxnsIndxs(i);  
         rxnID       = model.rxns(indx);
         mappedIndxs = rxnMapping(rxnID,model,false);
         FixedValues = [];
         range       = MAXmin_Optimizer(model,mappedIndxs,FixedValues);
         %If max and min were feasible then the optimization proceeds with
         %the ecModel
         if ~isempty(range)
            mappedIndxs = rxnMapping(rxnID,ecModel,true);
            FixedValues = basalFluxDist(mappedIndxs);
            rangeEC     = MAXmin_Optimizer(ecModel,mappedIndxs,FixedValues);
            if ~isempty(rangeEC)
                range_model   = [range_model; range];
                range_ecModel = [range_ecModel; rangeEC];
            end
         end
         disp(['ready with #' num2str(i)])
    end
     %Plot FV cumulative distributions
     %distributions          = {range_model(find(range_model)), range_ecModel(find(range_ecModel))};
     distributions          = {range_model, range_ecModel};
     legends                = {'Yeast model', 'ecYeast_ v1.5'};
     title                  = 'Flux variability cumulative distribution';
     [y_param, stats_param] = plotCumDist(distributions,legends,...
                               'Flux variability cumulative distribution');
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, irrevModel] = fixObjective(irrevModel,priorValue)
    % Optimize and fixes objective value for GEM
     objIndx  = find(irrevModel.c~=0);
     if nargin ==2
        irrevModel.lb(objIndx) = 0.99*priorValue;
        irrevModel.ub(objIndx) = priorValue;
     else
         sol               = solveLP(irrevModel);
         irrevModel.lb(objIndx) = 0.99*sol.x(objIndx);
         irrevModel.ub(objIndx) = sol.x(objIndx);
     end
     sol = solveLP(irrevModel,1);
     if ~isempty(sol.f)
         OptimalValue = sol.x(objIndx);
         optFluxDist  = sol.x;
     end
     disp(['The optimal value is ' num2str(OptimalValue)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_param, stats_param] = plotCumDist(E_parameter,legends,titlestr)
   figure
   for i=1:length(E_parameter)
        str{i} = horzcat(legends{i},' (',num2str(length(E_parameter{i})),...
                              ' / ',num2str(median(E_parameter{i})), ')');
        [y_param(i), stats_param(i)] = cdfplot(E_parameter{i});
        title(titlestr)
        ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
        xlabel('K_{cat} [s^{-1}]','FontSize',30,'FontWeight','bold');
        set(gca, 'XScale', 'log')
        hold on
   end
   legend(y_param,str);
end

