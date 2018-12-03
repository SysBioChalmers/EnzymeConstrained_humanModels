function [model,components] = setEMEMmedium(model,glucBound,oxBound,stdBound,irrev)
% setEMEMmedium
%
% Set a EMEM medium for the model, if the model is overconstrained
% (according to a provided experimental gRate measurement), then an
% iterative process that allows the uptake of new metabolites until the 
% model grows properly.
%
%   model       An HMR-based GEM
%   glucBound   Absolute value for the glucose uptake rate [mmol/gDw h]
%   oxBound     Absolute value for the oxygen uptake rate [mmol/gDw h]
%   stdBound    Absolute value used as a standard uptake rate [mmol/gDw h]
%   irrev       Indicates if a model is on its irreversible format 
%               (default = true)

%   model        EMEM medium constrained model
%   components   Components that were successfully added to the medium
%
% Ivan Domenzain.      Last edited: 2018-11-29
%
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';...
                  'glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'nicotinamide';'Pantothenate';...
                  'Pyridoxine';'Riboflavin';'Thiamin';'Glucose';'O2';'H2O';...
                  'Na+';'K+';'Mg+';'PI';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
                  'H+'};
if nargin<5
    irrev = true;
end
%getExchangeRxns works with the unconstrained field if present, lets remove
%it temporarily to avoid any inconsistency with the model bounds
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
end
unconstrained = zeros(length(model.mets),1);
%Get exchange rxn indexes and block all of them, except for growth, oxygen, 
%CO2 and the protein exchanges         
[~,excRxnIndxs] = getExchangeRxns(model);
GRindex         = find(model.c);
ProtIndex = find(contains(model.rxnNames,'prot_'));
CO2Index  = find(strcmpi(model.rxns,'HMR_9058'));
excRxnIndxs = setdiff(excRxnIndxs,GRindex);
excRxnIndxs = setdiff(excRxnIndxs,ProtIndex);
excRxnIndxs = setdiff(excRxnIndxs,CO2Index);
%Get an initial solution vector
sol = solveLP(model);
%block uptakes and production of exchanged metabolites
if irrev
    model.ub(excRxnIndxs) = 0;
else
    model.lb(excRxnIndxs) = 0;
end
%The model shouldn't be able to grow
sol = solveLP(model);
components = [];

for i=1:length(exchangeMets)
    excMetabolite = exchangeMets(i);   
    if ~isempty(find(strcmpi(model.metNames,excMetabolite)))
        for j=1:length(excRxnIndxs)
            %Then search the met exchange index (one by one)
            rxnIndx = excRxnIndxs(j);
            %Get the products compartment
            rxnMets  = model.metNames(find(model.S(:,rxnIndx)));
            if irrev
                prodComp = model.metComps(find(model.S(:,rxnIndx)>0));
            else
                prodComp = model.metComps(find(model.S(:,rxnIndx)<0));
            end
            
            if isempty(prodComp)
                prodComp = 10;
            end
            %If the searched metabolite matches with one of exchange rxn 
            %mets then find its uptake rxn and enable it 
            metExcIndx = find(strcmpi(rxnMets,excMetabolite));
            if ~isempty(metExcIndx) && strcmpi(model.compNames(prodComp),'extracellular')
                components = [components; excMetabolite];
                unconstrained(find(model.S(:,rxnIndx))) = 1;
                %Allow met uptake
                if irrev
                    model.ub(rxnIndx) = stdBound;
                else
                    model.lb(rxnIndx) = -stdBound;
                end
                
                if strcmpi(excMetabolite,'glucose')
                    if irrev
                        model.ub(rxnIndx) = glucBound;
                    else
                        model.lb(rxnIndx) = -glucBound;
                    end
                end
                if strcmpi(excMetabolite,'O2')||strcmpi(excMetabolite,'H2O')||strcmpi(excMetabolite,'H+')
                    if irrev
                        model.ub(rxnIndx) = oxBound;
                    else
                        model.lb(rxnIndx) = -oxBound;
                    end
                end
                disp([excMetabolite{1} ' added to the medium'])
            end
        end
    end
end
%Allow all secretions
[~,excRxnIndxs] = getExchangeRxns(model,'out');
model.ub(excRxnIndxs) = 1000;

%Is the model growing now?
sol = solveLP(model)
model.unconstrained = unconstrained;
end