function [Solution,optFlag] = solveMutant(model,refModel,method,prots)
% solvemodel
%
% Function that receives a model structure and gets a flux distribution
% distribution according to the selected optimization method. 'pFBA'
% performs a first FBA optimization, then fixes the optimal value for the
% objective function and performs a minimization of the total sum of fluxes
% (minimization of total protein usage for ecModels) for getting rid of
% loops. MOMA performs a quadratic programming simulation, minimizing the
% euclidean distance from the reference flux distribution, a reference model 
% is required in this case.
% 
% model        GEM matlab structure
% refModel     A reference model for calculating minimal adjustment if MOMA
%              is selected
% method       'MOMA' or 'pFBA'
% prots        Inxedes for the protein exchange reactions (ecrefModels)
%
% Solution     Optimal flux distribution for the model refModel
% optFlag      TRUE if sumulation was succesful
%
% Created.  Ivan Domenzain 2018-10-18
%

minFlag = false;
optFlag = 0;
if strcmpi(method,'MOMA')
    [Solution,~, optFlag]=qMOMA(model,refModel,1.01);
    disp(optFlag)
elseif strcmpi(method,'pFBA')
	%Get a first optimization
    Solution = solveLP(model);
    Solution = Solution.x;
    if ~isempty(Solution) && any(Solution)
        obj = find(model.c);
        if ~isempty(prots)
            %Fix optimal value for the objective and then minimize the total
            %sum of protein usages
            model.lb(obj)  = 0.999*Solution(obj);
            model.ub(obj)  = Solution(obj);
            model.c(:)     = 0;
            model.c(prots) = -1;
            Solution       = solveLP(model,1);
            if ~isempty(Solution.x) && any(Solution.x)
                Solution = Solution.x;
                minFlag  = true;
            end
        end
        %If protein usage minimization does not work or refModel is not EC, then 
        %apply simple pFBA
        if ~minFlag
            Solution = solveLP(model,1);
            Solution = Solution.x;
        end
        optFlag = 1;
    else 
        Solution = zeros(length(model.rxns));
    end
end
end