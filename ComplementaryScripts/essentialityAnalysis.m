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
sol_matrix     = sparse(length(model.rxns),length(model.genes));
infeasible_idx = [];
zero_grwt      = [];
gwrt_rate      = zeros(length(model.genes));

sol_matrix     = sparse(length(model.rxns),length(model.genes));
infeasible_idx = [];
zero_grwt      = [];
gwrt_rate      = zeros(length(model.genes));

for i=1:length(model.genes)
    gene    = model.genes(i);
    model_i = model;
    model_i = removeGenes(model_i,gene,false,false);
    sol_mut = solveLP(model_i,1);
    gwrt_rate(i) = abs(sol_mut.f);
    if ~isempty(sol_mut.f)
        sol_matrix(:,i) = sol_mut.x;
       change           = (abs(solution_wt.f)-abs(sol_mut.f))/abs(solution_wt.f);
       if change>0.01 && change <0.95
           disp(['Growth affected by gene', char(gene)])
       elseif change >= 0.95      
           zero_grwt      = [zero_grwt;i];
           disp(['Essential gene found: ', char(gene)])
       end     
    else 
        infeasible_idx = [infeasible_idx;i];
        disp(['Essential gene found: ', char(gene)])
    end
   
    disp(['Ready with gene deletion #', num2str(i)])
end
results.essential  = zero_grwt;
results.infeasible = infeasible_idx;
results.sol_matrix = sol_matrix;
results.wTypeSol   = solution_wt.x;
cd ../Results
save('essentialKO_HepG2','results')
end
 