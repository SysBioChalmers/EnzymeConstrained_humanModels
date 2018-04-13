%%

%load('HepG2model_modified.mat')% From models/HepG2/
model  = HepG2model_modified;
model = substituteBiomassRxns(model,false);

model        = HepG2model;
%solution_wt = solveLP(model,1);

sol_matrix     = sparse(length(model.rxns),length(model.genes));
infeasible_idx = [];
zero_grwt      = [];
gwrt_rate      = zeros(length(model.genes));

for i=1:20%length(model.genes)
    gene    = model.genes(i);
    model_i = model;
    model_i = removeGenes(model_i,gene,false,false);
    sol_mut = solveLP(model_i,1);
    gwrt_rate(i) = abs(sol_mut.f);
    if ~isempty(sol_mut.f) && (sol_mut.f < 0)
        sol_matrix(:,i) = sol_mut.x;
    
    elseif isempty(sol_mut.f)
        infeasible_idx = [infeasible_idx;i];
    elseif sol_mut.f == 0    
        zero_grwt      = [zero_grwt;i];
    end
    disp(i)
end

save('essentialKO_HepG2','zero_grwt','infeasible_idx','sol_matrix')
 