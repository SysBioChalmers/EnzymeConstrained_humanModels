function model = modelModifications(model)
[GRR, RGM]       = standardizeGrRules(model);
model.grRules    = GRR;
model.rxnGeneMat = RGM;
model.b          = zeros(length(model.mets),1);
model            = setParam(model,'obj','biomass_human',1);
end