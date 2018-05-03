function [pairwiseComp,enzUsages] = fluxDistMatching(model,ecModel)
% Get flux distributions from pFBA solutions 
sol_WT     = solveLP(model,1);
sol_WT     = sol_WT.x;
sol_mutant = solveLP(ecModel,1);
sol_mutant = sol_mutant.x;
% Get size of both models
[~,n]     = size(model.S);

pairwiseComp = zeros(n,2); 

for i=1:n
    rxnID       = model.rxns(i);
    revFlag     = logical(model.rev(i));
    mappedIndxs = rxnMapping(rxnID,ecModel,revFlag);
    if ~isempty(mappedIndxs)
        if length(mappedIndxs) == 2
            revIndex = find(contains(ecModel.rxns(mappedIndxs),'_REV'));
            if ~isempty(revIndex)
                ecFlux_i = sum(sol_mutant(mappedIndxs))-2*sol_mutant(revIndex);
            else
                ecFlux_i = sum(sol_mutant(mappedIndxs));
            end
        else
            ecFlux_i = sol_mutant(mappedIndxs);
        end
    else
        disp(['Rxn #' num2str(i) ' not found in ecModel'])
        pause
    end
    Flux_i = sol_WT(i);
    
    pairwiseComp(i,1) = Flux_i; 
    pairwiseComp(i,2) = ecFlux_i; 
    disp(['Ready with rxn #' num2str(i)]) 
end

enzUsagesIndexes = contains(ecModel.rxns,'draw_prot');
enzUsages        = sol_mutant(enzUsagesIndexes);
end


            
    
    
    
    

