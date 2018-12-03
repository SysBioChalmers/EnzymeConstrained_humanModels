function [pairwiseComp,enzUsages] = fluxDistMatching(model,ecModel)
% fluxDistMatching
%
% Script that gets a GEM and its EC version, it performs a pFBA simulation
% on both and a pairwise flux distribution comparison is then obtained.
% This is done by the correct mapping of the original reactions onto the
% ecModel, taking directionality into account.
%
%   model            a Genome-scale mEtabolic Model
%   ecModel          Enzymatically-constrained version of the GEM (obtained by
%
%   pairwiseComp     Numerical matrix with the flux distribution of the
%                    original model (column 1) and the mapped ecModel 
%                    distribution (column 2)    
%   enzUsages        Cell array with enzyme usage reactions fluxes 
%                    (mmol/gDwh) for each enzyme in the model.
%
% Usage: [pairwiseComp,enzUsages] = fluxDistMatching(model,ecModel)
%
% Ivan Domenzain.      Last edited: 2018-06-12
%
% Get flux distributions from pFBA solutions 
sol_GEM      = solveLP(model,1);
sol_GEM      = sol_GEM.x;
sol_EC       = solveLP(ecModel,1);
sol_EC       = sol_EC.x;
[~,n]        = size(model.S);
pairwiseComp = zeros(n,2); 

for i=1:n
    rxnID       = model.rxns(i);
    revFlag     = logical(model.rev(i));
    mappedIndxs = rxnMapping(rxnID,ecModel,revFlag);
    if ~isempty(mappedIndxs)
        if length(mappedIndxs) == 2
            revIndex = find(contains(ecModel.rxns(mappedIndxs),'_REV'));
            if ~isempty(revIndex)
                ecFlux_i = sum(sol_EC(mappedIndxs))-2*sol_EC(revIndex);
            else
                ecFlux_i = sum(sol_EC(mappedIndxs));
            end
        else
            ecFlux_i = sol_EC(mappedIndxs);
        end
    else
        disp(['Rxn #' num2str(i) ' not found in ecModel'])
        pause
    end
    Flux_i = sol_GEM(i);
    
    pairwiseComp(i,1) = Flux_i; 
    pairwiseComp(i,2) = ecFlux_i; 
    disp(['Ready with rxn #' num2str(i)]) 
end

enzUsagesIndexes = contains(ecModel.rxns,'prot_');
enzUsagesIndexes = enzUsagesIndexes(1:end-1);
enzUsages        = sol_EC(enzUsagesIndexes);
end


            
    
    
    
    

