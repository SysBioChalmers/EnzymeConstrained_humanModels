function results =  essentialityAnalysis(model,lethalTreshold,filename,EMEM)
% results =  essentialityAnalysis(model,lethalTreshold,filename)
%
% Script that performs a gene essentiallity analysis on a given GEM.
% INPUT:
%   model           GEM matlab structure
%   lethalTreshold  growth reduction threshold for a gene to be considered
%                   as lethal
% OUTPUTS:
%   results.essential   cell array with the essential genes indexes
%                       (zero-Growth).
%   results.infeasible  cell array with the genes indexes that yielded
%                       infeasible solutionswhen removed.
%   results.sol_matrix  sparse matrix with the flux distributions of the
%                       optimal solution for each deletion.
%   results.wTypeSol    Optimal solution vector for the "wild_type" model
%
% usage: results =  essentialityAnalysis(model,lethalTreshold,filename)
%
% Raphael Ferreira.     Created:     2018-04-13
% Ivan Domenzain.       Last edited: 2018-10-15
%

current = pwd;
%Set EMEM medium
if EMEM
     [model,~,~] = setEMEMmedium(model,0.01);
end
  solution_wt    = solveLP(model,1);
objIndex       = find(model.c==1);
growth_wt      = solution_wt.x(objIndex);
essential_idx  = [];
affected_grwt  = [];
gwrt_rate      = zeros(length(model.genes));
objective      = find(model.c==1);
%Try to find the indexes for individual protein usage reactions
  protIndexes    = find(contains(model.rxnNames,'draw_prot'));

for i=1:length(model.genes)
    gene    = model.genes(i);
    model_i = model;
    model_i = removeGenes(model_i,gene,false,false);
    sol_mut = solveLP(model_i,1);

    if ~isempty(sol_mut.f)
        gwrt_rate(i) = abs(sol_mut.f);
        growth_mut   = sol_mut.x(objective);
        %sol_matrix(:,i) = sol_mut.x;
        change       = (growth_wt-growth_mut)/growth_wt;
        %Avoid numerical inaccuracies
        if abs(change)<=1e-5
            change = 0;
        end
        %Identify genes that affect growth rate
        if change >= 0.05
            affected_grwt = [affected_grwt;i];
            disp(['Growth affected by gene #' num2str(i) ' (' char(gene) ')'  ': reduction = ' num2str(change*100) '%'])
            if change >= lethalTreshold
                essential_idx = [essential_idx;i];
                disp(['Essential gene found: ', char(gene)])
            end
        end
    %Unfeasible mutants are considered as essential
    else 
        essential_idx = [essential_idx;i];
        disp(['Essential gene found: ', char(gene)])
    end
    disp(['Ready with gene deletion #', num2str(i)])
end
%Create results table with flux distributions for the growth reducing
%genes based on the principal of protein burden minimization
sol_matrix = simulateTargets(model,affected_grwt,protIndexes);
sol_matrix = [solution_wt.x, sol_matrix];
sol_matrix = [model.rxns, num2cell(sol_matrix)];
names      = [{'Rxns'},{'Wild_Type'}, model.genes{affected_grwt}];
Results    = cell2table(sol_matrix,'VariableNames',names,'RowNames',model.rxns);
cd ../../Results
writetable(Results,filename,'Delimiter','\t')
results.affected   = affected_grwt;
results.essential  = essential_idx;
results.sol_matrix = sol_matrix;
results.wTypeSol   = solution_wt.x;
cd(current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol_matrix = simulateTargets(model,affected_grwt,prots)
sol_matrix = zeros(length(model.rxns),length(affected_grwt));
genes      = model.genes(affected_grwt);
for i=1:length(affected_grwt)
    gene = genes(i);
    disp(['Getting Flux distribution for gene: ',gene{1}])
    mutant          = removeGenes(model,gene,false,false);
    [Solution,~]    = solveMutant(mutant,[],'pFBA',prots);
    sol_matrix(:,i) = Solution;
end
end
 