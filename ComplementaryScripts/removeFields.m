function model_data = removeFields(model_data)
% Remove unnecessary fields for ecModels that cause conflicts with
% COBRA/RAVEN addRxn functions
%
% Ivan Domenzain.   2018-10-07
%
model = model_data.model;

if isfield(model,'rxnFrom')
    model = rmfield(model,'rxnFrom');
end

if isfield(model,'metFrom')
    model = rmfield(model,'metFrom');
end

if isfield(model,'geneFrom')
    model = rmfield(model,'geneFrom');
end

if isfield(model,'rxnReferences')
    model = rmfield(model,'rxnReferences');
end

if isfield(model,'rxnConfidenceScores')
    model = rmfield(model,'rxnConfidenceScores');
end

if isfield(model,'rxnRecon3DID')
    model = rmfield(model,'rxnRecon3DID');
end

if isfield(model,'rxnRecon3DID')
    model = rmfield(model,'rxnRecon3DID');
end

model_data.model = model;

end