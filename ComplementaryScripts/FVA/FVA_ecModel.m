%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
% Raphael Ferreira     Last edited: 2018-03-14
% Ivan Domenzain.      Last edited: 2018-03-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
    range_model      = [];
    optimizedIndxs   = [];
    optimizedECIndxs = [];
    range_ecModel    = [];    
    %Block glucose and oxygen production
	%if nargin>2
		Blocked_mets    = {'D-glucose', 'oxygen'};
		model   = block_production(model,Blocked_mets,false);
        ecModel = block_production(ecModel,Blocked_mets,true);
    %end
     %Get the index for all the non-objective rxns 
     rxnsIndxs    = find(model.c~=1);     
     %Gets the optimal value for ecModel and fixes the objective value to 
     %this for both models
     [OptimalValue, ~, ecModel]       = fixObjective(ecModel);
     [OptimalValue,optFluxDist,model] = fixObjective(model,OptimalValue);
     % Get the variability range for each of non-objective reactions in the
     % original model
     for i=1:length(rxnsIndxs)
         rangeEC = [];
         indx    = rxnsIndxs(i);         
         % performs the analysis just on rxns with a non-zero flux
         %if optFluxDist(indx)~=0
             %Maximize i-th rxn
             sol = Optimizer(model,indx,1);
             %If the solution was feasible then proceed to the minimization
             %of the rxn
             if ~isempty(sol.f)
                 maxFlux = sol.x(indx);
                 %minimization
                 sol = Optimizer(model,indx,-1);
                 %If minimization was feasible then proceed to calculate
                 %the FV range
                 if ~isempty(sol.f)
                     minFlux = sol.x(indx);
                     range   = maxFlux- minFlux;
                     %Now for the EC model (irrev model)
                     rxnID       = model.rxns(indx);
                     mappedIndxs = rxnMapping(rxnID,ecModel,true);
                     if length(mappedIndxs)>1
                     %if model.lb(indx) < 0 % model.rev(indx) == 1%
                         %mappedIndxs = rxnMapping(rxnID,ecModel,true);
                         %if length(mappedIndxs)<2
                         %    disp(rxnID)
                         %end
                         %Set objective function for the maximization of the
                         %forward rxn
                         sol = Optimizer(ecModel,mappedIndxs(1,1),1);
                         % If the maximization was feasible proceed to min
                         if ~isempty(sol.f)
                            maxECFlux = sol.x(indx);
                           %Set objective function for the maximization 
                           %of thebackward rxn (minimization) and Fix 
                           %forward rxn flux to "zero"
                           sol = Optimizer(ecModel,mappedIndxs(2,1),1,...
                                                mappedIndxs(1,1),0);
                            %If max and min were feasible the FV range for
                            %the i-th rxn (mapped in the ecModel) is
                            %calculated
                            if ~isempty(sol.f)
                                minECFlux = sol.x(indx);
                                rangeEC   = maxECFlux+minECFlux;
                            end                       
                         end
                     else
                         mappedIndxs = rxnMapping(rxnID,ecModel,false); 
                        %Max
                         sol = Optimizer(ecModel,mappedIndxs,1);
                         if ~isempty(sol.f)
                            maxECFlux = sol.x(indx);
                            %Min
                            sol = Optimizer(ecModel,mappedIndxs,-1);
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
                 disp(['ready with ' string(rxnID)])
             end
         %end
     end

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, model] = fixObjective(model,priorValue)
    % Optimize and fixes objective value for GEM
     objIndx  = find(model.c~=0);
     if nargin ==2
        model.lb(objIndx) = 0.99*priorValue;
        model.ub(objIndx) = priorValue;
     else
         sol               = solveLP(model);
         model.lb(objIndx) = 0.99*sol.x(objIndx);
         model.ub(objIndx) = sol.x(objIndx);
     end
     sol = solveLP(model,1);
     if ~isempty(sol.f)
         OptimalValue = sol.x(objIndx);
         optFluxDist  = sol.x;
     end
     disp(['The optimal value is ' num2str(OptimalValue)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = Optimizer(model,indx,coeff,blockIndx,FixedFlux)
             model.c       = zeros(length(model.c),1);
             model.c(indx) = coeff;
             if nargin > 3
                 model.ub(blockIndx) = FixedFlux;
             end
             sol = solveLP(model);
end
