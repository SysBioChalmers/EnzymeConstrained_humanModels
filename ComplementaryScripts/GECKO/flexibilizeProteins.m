function [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,glucUptakeExp,glucUptkIndx,glnIndex,lacIndex,gluIndx)
% flexibilizeProteins
%   Function that takes an ecModel with proteomic constraints and, if it is
%   overconstrained with respect to the provided experimental growth rate,
%   iterates finding the top growth-limiting enzyme or enzyme complex. The
%   exchange rate upper bound for each of the identified enzyme or subunits 
%   is then set to infinity and after all iterations they are set equal to 
%   the usage value provided by an overall enzyme usage minimization 
%   simulation (subject to the provided growth rate and nutrient uptake 
%   constraints).
%   
%   model          ecModel with proteomic constraints (individual enzyme
%                  levels)
%   gRate          Growth rate [g/gDw h]observed on the experiment in which
%                  the proteomic measurements were taken
%   glucUptakeExp  (Opt) Experimentally measured glucose uptake rate 
%                  [mmol/gDw h]
%
%   model          ecModel with calibrated enzyme usage upper bounds
%   enzUsages      Calculated enzyme usages after final calibration 
%                  (enzyme_i demand/enzyme_i upper bound)
%   modifications  Table with all the modified values 
%                  (Protein ID/old value/Flexibilized value)
%
%   Usage: [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,glucUptakeExp)
%
%   Ivan Domenzain, 2018-06-11
%
current       = pwd;
flexFactor    = 100;
flexProts     = {};
enzUsages     = [];
modifications = {};

%constrain glucose uptake if an experimental measurement is provided
if nargin>2
    model.ub(glucUptkIndx) = 1.001*glucUptakeExp;
    model.ub(glnIndex) = 0;
    model.ub(lacIndex) = 0.843;
    model.ub(gluIndx) = 0;
end
cd (current)
% get measured protein exchange rxns indexes
measuredIndxs = getMeasuredProtsIndexes(model);
if ~isempty(measuredIndxs)
    abundances    = model.ub(measuredIndxs);
    objIndex      = find(model.c==1);
    sol           = solveLP(model,1);
    iteration = 1;
    while isempty(sol.x)
        model.ub(glnIndex) = model.ub(glnIndex)+0.01;
        sol                = solveLP(model,1);
        disp(['Iteration #' num2str(iteration) ' not feasible'])
        iteration = iteration+1;
    end
end    
        
        
%         
%     growth        = sol.x(objIndex);
%     % iterate while growth is underpredicted
%     while growth<0.99*gRate
%         [limIndex,flag] = findLimitingUBs(model,measuredIndxs,flexFactor,1);
%         if ~flag
%             [limIndex,~] = findLimitingUBs(model,measuredIndxs,flexFactor,2);
%         end
%         %Flexibilize the top growth limiting protein on the original ecModel
%         flexProts          = [flexProts; model.rxns(limIndex)];
%         model.ub(limIndex) = Inf;
%         sol                = solveLP(model);
%         if ~isempty(sol.x)
%            growth = sol.x(objIndex);
%            for j=1:length(limIndex)
%                 indx = limIndex(j);
%                 disp(['Modified ub for: ' model.rxnNames{indx} ' gRate: ' num2str(growth)])
%            end
%         end
%     end
%     [model,enzUsages]  = getNewBounds(model,gRate,measuredIndxs,flexProts,objIndex);
%     modifiedAbundances = model.ub(measuredIndxs);
%     exchangedProteins  = model.rxnNames(measuredIndxs);
%     modifications      = getDifferences(abundances,modifiedAbundances,exchangedProteins);
% end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffTable = getDifferences(originalBounds,newBounds,exchangedProteins)
protein_IDs     = {};
previous_values = {};
modified_values = {};
for i=1:length(originalBounds)
    if newBounds(i)>originalBounds(i)
        proteinID       = exchangedProteins{i};
        proteinID       = proteinID(1:((strfind(proteinID,'_exchange'))-1));
        protein_IDs     = [protein_IDs; proteinID];
        previous_values = [previous_values; originalBounds(i)];
        modified_values = [modified_values; newBounds(i)];
    end
end
diffTable = table(protein_IDs,previous_values,modified_values);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measuredIndxs = getMeasuredProtsIndexes(model)
measuredIndxs  = find(contains(model.rxnNames,'prot_'));
exchange_prots = find(contains(model.rxnNames(measuredIndxs),'_exchange'));
measuredIndxs  = measuredIndxs(exchange_prots(1:end-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,enzUsages] = getNewBounds(model,gRate,protIndxs,flexProts,gPos)
%Now that the model is growing at least at the specified dilution rate
%lets fix the growth rate and minimize enzymes usage
objectiveVector      = model.c;
model.lb(gPos)       = 0.99*gRate;
model.ub(gPos)       = 1.01*gRate;
model.c(:)           = 0;
protIndexes          = contains(model.rxnNames,'prot_');
model.c(protIndexes) = -1;
optSolution          = solveLP(model,1);
optSolution          = optSolution.x;
enzUsages            = zeros(length(protIndxs),1);
for i=1:length(protIndxs)
    index = protIndxs(i);
    name  = model.rxns(index);
    %If protein was flexibilized set its upper bound to the simulated
    %concentration
    if ismember(name,flexProts) && optSolution(index)>0
        model.ub(index) = optSolution(index);
    end
    enzUsages(i) = optSolution(index)/model.ub(index);
end
model.c   = objectiveVector;
enzUsages = enzUsages(enzUsages > 0);
end
    
    
