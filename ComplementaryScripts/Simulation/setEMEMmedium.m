function [model,components] = setEMEMmedium(model)
% find_EMEM_ExchangeRxns
%
% Set a EMEM medium for the model, if the model is overconstrained
% (according to a provided experimental gRate measurement), then an
% iterative process that allows the uptake of new metabolites until the 
% model grows properly.
%
% INPUT:
%   model       An HMR-based GEM
% OUTPUTS:
%   model       EMEM medium constrained model
%   essential   Essential components that were missing in the initial EMEM
%               formulation
%
% Ivan Domenzain.      Last edited: 2018-10-26
%
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';...
                  'glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'nicotinamide';'Pantothenate';...
                  'Pyridoxine';'Riboflavin';'Thiamin';'Glucose';'O2';'H2O';...
                  'Na+';'K+';'Mg+';'PI';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
                  'H+'};
%getExchangeRxns works with the unconstrained field if present, lets remove
%it temporarily to avoid any inconsistency with the model bounds
model = rmfield(model,'unconstrained');
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
model.ub(excRxnIndxs) = 0;
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
            prodComp = model.metComps(find(model.S(:,rxnIndx)>0));
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
                model.ub(rxnIndx) = 1000;
                if strcmpi(excMetabolite,'glucose')
                    model.ub(rxnIndx) = 1;
                end
                if strcmpi(excMetabolite,'O2')||strcmpi(excMetabolite,'H2O')||strcmpi(excMetabolite,'H+')
                    model.ub(rxnIndx) = 1000;
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