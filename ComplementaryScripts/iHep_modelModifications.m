function model = iHep_modelModifications(model)
model = modifyRxnDependentFields(model);
model = substitute_grRules(model);
model = substituteEnsemblGeneIDs(model);
model.geneNames = model.genes;
% Standardizes the metabolites and rxns names
%model  = modifyMetNames(model);
model = modifyMetNames_iHep(model);
% Substitute biomass reaction
%model  = substituteBiomassRxns(model,false);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modifyRxnDependentFields(model)
n        = length(model.rxnNames);
rxnNames = cell(n,1);
rules    = model.rules;
grRules  = cell(n,1);
rev      = false(n,1);
LBs      = model.lb;
for i=1:n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  rxns  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rxnName  = model.rxnNames{i};
    position = strfind(rxnName,'R_');
    if ~isempty(position)
        position = position(1);
        rxnName  = rxnName(position+2:end);
    end
    rxnNames{i} = rxnName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% grRules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(rules{i})
        rule = rules{i};
        rule = strsplit(rule,'|');
        newRule = []; 
        for j=1:length(rule)
            gene    = strrep(rule{j},' ','');
            gene    = strrep(gene,'x','');
            gene    = strrep(gene,'(','');
            gene    = strrep(gene,')','');
            gene    = model.genes{str2double(gene)};
            newRule = [newRule gene];
            if j<length(rule)
                newRule = [newRule ' or '];
            end
        end
        disp(num2str(i))
        disp(newRule)
        grRules{i} = newRule;
    else
        grRules{i} = '';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% reversibility %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if LBs(i)<0
        rev(i) = true;
    end
end
model.rxnNames        = rxnNames;
model.grRules         = grRules;
[grRules, rxnGeneMat] = standardizeGrRules(model);
model.grRules         = grRules;
model.rxnGeneMat      = rxnGeneMat;
model.rev             = rev;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modifyMetNames_iHep(model)

    m           = length(model.metNames);
    metNames    = cell(m,1);
    mets        = cell(m,1);
    metformulas = cell(m,1);
    metComps    = cell(m,1);
    for i=1:m
        metName       = model.metNames{i};
        metabolite    = model.mets{i};
        underScorePos = strfind(metName,'_');
        bracketPos    = strfind(metabolite,'[');
        %%%%%%%%%%%%%%%%%%%%%%%%%%metNames%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(underScorePos)    
            underScorePos  = underScorePos(end);
            formula        = metName(underScorePos+1:end);
            metName        = metName(1:underScorePos-1);
            metformulas{i} = formula;
        end        
        if any(strfind(metName,'[protein]-'))
            metName = replace(metName,'[protein]-','');
            
        elseif(any(strfind(metName,'-[protein]')))
            metName = replace(metName,'-[protein]','');
            
        elseif(any(strfind(metName,'[protein C terminal]-')))
            metName = replace(metName,'[protein C terminal]-','');
        
        end
        metNames{i} = metName;
        %%%%%%%%%%%%%%%%%%%%%%%%%%metNames%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(bracketPos)   
            bracketPos  = bracketPos(end);
            compartment = metabolite(bracketPos+3);
            metabolite  = metabolite(1:bracketPos-1);
            compartment = find(strcmpi(model.comps,compartment));
            if ~isempty(compartment)
                metComps{i} = compartment;
            end
            mets{i}     = metabolite;
        end
    end
    model.mets        = mets;
    model.metNames    = metNames;
    model.metComps    = metComps;
    model.metformulas = metformulas;
end

    