%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
% Raphael Ferreira     Last edited: 2018-03-09
% Ivan Domenzain.      Last edited: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
    range_model      = [];
    optimizedIndxs   = [];
    optimizedECIndxs = [];
    range_ecModel    = [];
    
	%if nargin>2
		Blocked_mets    = {'glucose', 'oxygen'};
		model   = block_production(model,Blocked_mets,false);
        ecModel = block_production(ecModel,Blocked_mets,true);
    %end
     
     rxnsIndxs  = find(model.c~=1);
     objIndx    = find(model.c==1);
     % Optimize for the objective function and "fix" its bounds to the
     % optimal values
     sol = solveLP(model,1);
     model.lb(objIndx) = -0.999*sol.f;
     model.ub(objIndx) = -sol.f;
     disp(-sol.f)
     % Optimize and fixes objective value for ECmodel
     objIndx    = find(ecModel.c==1);
     sol = solveLP(ecModel,1);
     ecModel.lb(objIndx) = -0.999*sol.f;
     ecModel.ub(objIndx) = -sol.f;
     disp(-sol.f)
     
     % Get the variability range for each of the rxns
     for i=1:length(rxnsIndxs)
         rangeEC = [];
         indx = rxnsIndxs(i);
         model.c = zeros(length(model.c),1);
         %Maximization
         model.c(indx) = 1;
         sol = solveLP(model);
         if ~isempty(sol.f)
             maxFlux = sol.x(indx);
             %Set objective for minimization
             model.c(indx) = -1;
             sol = solveLP(model);
             if ~isempty(sol.f)
                 minFlux = sol.x(indx);
                 range = maxFlux- minFlux;
                 rxnID          = model.rxns(indx);
    
                 % Rxn ID mapping to the EC model 
                 if model.rev(indx) == 1
                     mappedIndxs = rxnMapping(rxnID,ecModel,true);
                     %Max
                     ecModel.c = zeros(length(ecModel.c),1);
                     ecModel.c(mappedIndxs(1)) = 1;
                     sol = solveLP(ecModel);
                     % If the maximization was feasible
                     if ~isempty(sol.f)
                        maxECFlux = sol.x(indx);
                        ecModel.c = zeros(length(ecModel.c),1);
                        ecModel.c(mappedIndxs(2))  = 1;
                        % Fix forward rxn flux to "zero"
                        ecModel.ub(mappedIndxs(2)) = 0.001;
                        %Min
                        sol = solveLP(ecModel);
                        if ~isempty(sol.f)
                            minECFlux = sol.x(indx);
                            rangeEC = maxECFlux+minECFlux;
                        end                       
                     end
                 else
                     mappedIndxs = rxnMapping(rxnID,ecModel,false); 
                    %Max
                     ecModel.c = zeros(length(ecModel.c),1);
                     ecModel.c(mappedIndxs) = 1;
                     sol = solveLP(ecModel);
                     if ~isempty(sol.f)
                        maxECFlux = sol.x(indx);
                        %Min
                        ecModel.c = zeros(length(ecModel.c),1);
                        ecModel.c(mappedIndxs) = -1;
                        sol = solveLP(ecModel);
                         if ~isempty(sol.f)
                            minECFlux = sol.x(indx);
                            rangeEC = maxECFlux+minECFlux;
                         end
                     end
                 end
                 
                 if ~isempty(rangeEC)
                     range_model    = [range_model; range];
                     optimizedIndxs = [optimizedIndxs; indx];
                     range_ecModel    = [range_ecModel; rangeEC];
                 end
                 
             end
             disp(['ready with ' rxnID])
         end
     end

%
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mappedIndxs = rxnMapping(rxnID,model,revFlag)
    indexes = find(contains(model.rxns,rxnID));
    if revFlag
        backwardIndxs = find(contains(model.rxns(indexes),'_REV'));
        backwardIndxs = findArmRxns(backwardIndxs,model);
        forwardIndxs  = setdiff(indexes,backwardIndxs);
        forwardIndxs  = findArmRxns(forwardIndxs,model);  
        mappedIndxs   = vertcat(forwardIndxs,backwardIndxs);
    else
        mappedIndxs  = findArmRxns(indexes,model);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ArmIndex = findArmRxns(rxnIndexes,model)
    if length(rxnIndexes)>1
       ArmIndex = rxnIndexes(contains(model.rxns(rxnIndexes),'arm_')); 
    else
       ArmIndex = rxnIndexes;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = block_production(model,metName,revFlag)
	% COBRA function for the extraction of Exchange reaction indexes
    [~,Uptake_indxs] = findExcRxns(model);
    Uptake_indxs 	 = find(Uptake_indxs);
    % Gets the indexes of all the extracellular mets
    compIndex		 = find(strcmpi(model.compNames,'extracellular'));
    compMets		 = find(model.metComps==compIndex);
    
    %For each metabolite uptake that is going to be blocked
    for i=1:length(metName)
        % Gets the index for the metabolite locallized in the extracellular
        % compartment
        ExtMet_Index = compMets(strcmpi(model.metNames(compMets),metName(i)));
        
        if length(ExtMet_Index)>1
            disp(['warning: multiple matches were found in the compartment for ' string(metName(i))])
        else
            
            % Get uptake rxns indexes for the especified metabolite
            metRxns      = find(model.S(ExtMet_Index,:));
            if revFlag==false
                metUptakeRxn = intersect(metRxns,Uptake_indxs);
                % Block metabolite secretion
                if ~isempty(metUptakeRxn)
                    model.lb(metUptakeRxn) = 0;
                end
            else
                Uptake_indxs = find(model.S(ExtMet_Index,:)<0);
                if ~isempty(Uptake_indxs)
                    for j=1:length(Uptake_indxs)
                        product = model.metNames(find(model.S(:,Uptake_indxs(j))>0));
                        if isempty(product)
                            break
                        end
                    end
                model.ub(Uptake_indxs(j)) = 0;
                end                
            end
        end
    end
end