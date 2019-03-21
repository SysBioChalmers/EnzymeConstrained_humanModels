function model = modelModifications(model)
[grRules, rxnGeneMat] = standardizeGrRules(model);
model.grRules         = grRules;
model.rxnGeneMat      = rxnGeneMat;
model.b               = zeros(length(model.mets),1);
%model                 = substitute_grRules(model);
%disp(model.grRules)
%model                 = substituteEnsemblGeneIDs(model);
% Standardizes the metabolites names
model  = modifyMetNames(model);
% Substitute biomass reaction
model  = substituteBiomassRxns(model,false);
end