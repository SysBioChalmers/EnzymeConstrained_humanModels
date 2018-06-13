%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results =  essentialityAnalysis(model)
%
% Script that performs a gene essentiallity analysis on a given GEM.
% INPUT:
%   model       GEM matlab structure
% OUTPUTS:
%   results.essential   cell array with the essential genes indexes
%                       (zero-Growth).
%   results.infeasible  cell array with the genes indexes that yielded
%                       infeasible solutionswhen removed.
%   results.sol_matrix  sparse matrix with the flux distributions of the
%                       optimal solution for each deletion.
%   results.wTypeSol    Optimal solution vector for the "wild_type" model
%
% Raphael Ferreira.     Created:     2018-04-13
% Ivan Domenzain.       Last edited: 2018-04-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results =  essentialityAnalysis(model)

solution_wt    = solveLP(model,1);
objIndex       = find(model.c==1);
growth_wt      = solution_wt.x(objIndex);
sol_matrix     = sparse(length(model.rxns),length(model.genes));
essential_idx  = [];
affected_grwt  = [];
gwrt_rate      = zeros(length(model.genes));


for i=1:length(model.genes)
    gene    = model.genes(i);
    model_i = model;
    model_i = removeGenes(model_i,gene,false,false);
    sol_mut = solveLP(model_i,1);
    gwrt_rate(i) = abs(sol_mut.f);
    if ~isempty(sol_mut.f)
        objMut          = find(model_i.c==1);
        growth_mut      = sol_mut.x(objMut);
        sol_matrix(:,i) = sol_mut.x;
       change           = abs(growth_wt-growth_mut)/growth_wt;
       if change >= 0.05      
           affected_grwt = [affected_grwt;i];
           disp(['Growth affected by gene: ', char(gene)])
       elseif change >= 0.50 
           essential_idx = [essential_idx;i];
           disp(['Essential gene found: ', char(gene)])
       end     
    else 
        essential_idx = [essential_idx;i];
        disp(['Essential gene found: ', char(gene)])
    end
   
    disp(['Ready with gene deletion #', num2str(i) ': ' char(gene)])
end
results.affected   = affected_grwt;
results.essential  = essential_idx;
results.sol_matrix = sol_matrix;
results.wTypeSol   = solution_wt.x;
cd ../Results
save('essentialKO_ecHepG2','results')
end
 