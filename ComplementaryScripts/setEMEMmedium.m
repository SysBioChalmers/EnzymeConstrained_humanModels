function [model,essential] = find_EMEM_ExchangeRxns(model,gRate)
% find_EMEM_ExchangeRxns
%
% Set a EMEM medium for the model, if the model is overconstrained
% (according to a provided experimental gRate measurement), then an
% iterative process that allows the uptake of new metabolites until the 
% model grows properly.
%
% INPUT:
%   model       An HMR-based GEM
%   gRate       A experimentally measured growthrate for the cell-type
% OUTPUTS:
%   model       EMEM medium constrained model
%   essential   Essential components that were missing in the initial EMEM
%               formulation
%
% Ivan Domenzain.      Last edited: 2018-09-12
%
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';...
                  'glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'nicotinamide';...
                  'Pantothenate';'Pyridoxine';'Riboflavin';'Thiamin';...
                  'Glucose';'Oxygen';'H2O';'Na+';'K+';'Mg+';...
                  'PI';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
                  'H+';'cholesterol'};
%Get exchange rxn indexes and block all of them, except for growth, oxygen, 
%CO2 and the protein exchanges         
[~,excRxnIndxs] = getExchangeRxns(model);
GRindex         = find(strcmpi(model.rxns,'biomass_components'));
ProtIndex = find(contains(model.rxnNames,'prot_'));
CO2Index  = find(strcmpi(model.rxns,'HMR_9058'));
excRxnIndxs = setdiff(excRxnIndxs,GRindex);
excRxnIndxs = setdiff(excRxnIndxs,ProtIndex);
excRxnIndxs = setdiff(excRxnIndxs,CO2Index);
%Get an initial solution vector
sol = solveLP(model)
%block uptakes and production of exchanged metabolites
model.ub(excRxnIndxs) = 0;
model.lb(excRxnIndxs) = 0;
%The model shouldn't be able to grow
sol = solveLP(model)
mediumComponents = [];

for i=1:length(exchangeMets)
    excMetabolite = exchangeMets(i);
    
    if ~isempty(find(strcmpi(model.metNames,excMetabolite)))
        for j=1:length(excRxnIndxs)
            %Then search the met exchange index (one by one)
            rxnIndx = excRxnIndxs(j);
            %Get all the mets and compartments for a given exchange rxn
            rxnMets  = model.metNames(find(model.S(:,rxnIndx)));
            subs     = model.metNames(find(model.S(:,rxnIndx)<0));
            subComp  = model.metComps(find(model.S(:,rxnIndx)<0));
            prods    = model.metNames(find(model.S(:,rxnIndx)>0));
            prodComp = model.metComps(find(model.S(:,rxnIndx)>0));
            
            metExcIndx = find(strcmpi(rxnMets,excMetabolite));
            %If the searched metabolite matches with one of exchange rxn 
            %mets then find its uptake rxn and enable it 
            if ~isempty(metExcIndx)
                mediumComponents = [mediumComponents; i];
                %if isempty(subs)
                %Allow met exchange
                    model.ub(rxnIndx) = 1000;
                    disp([excMetabolite{1} ' added to the medium'])
                %end
            end
        end
    end
end
%Is the model growing now?
sol = solveLP(model)
priorValue = abs(sol.f);
%If overconstrained then look for the missing components for the model to
%grow
essential = [];
i = 1;
while priorValue<gRate
    rxnIndx = excRxnIndxs(i);
    if ~ismember(i,rxnIndx)
        temp_model = model;
        if temp_model.ub(rxnIndx) == 0
            temp_model.ub(rxnIndx) = 1000;
            sol = solveLP(temp_model);
            if abs(sol.f)>priorValue
                disp(sol.f)
                disp(['Essential component found: ' model.rxns(rxnIndx)])
                essential = [essential;rxnIndx];
            end
        end
    end
end
end