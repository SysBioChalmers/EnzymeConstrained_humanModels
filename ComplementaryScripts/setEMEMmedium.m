function [model,essential,sensitivities] = setEMEMmedium(model,gRate)
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
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';'Cysteine';...
                  'glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'nicotinamide';...
                  'Pantothenate';'Pyridoxine';'Riboflavin';'Thiamin';...
                  'Glucose';'O2';'H2O';'Na+';'K+';'Mg+';...
                  'PI';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
                  'H+'};
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
sol = solveLP(model)
%block uptakes and production of exchanged metabolites
model.ub(excRxnIndxs) = 0;
%The model shouldn't be able to grow
sol = solveLP(model)
mediumComponents = [];

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
            metExcIndx = find(strcmpi(rxnMets,excMetabolite));
            %If the searched metabolite matches with one of exchange rxn 
            %mets then find its uptake rxn and enable it 
            if ~isempty(metExcIndx) && prodComp==1
                mediumComponents = [mediumComponents; i];
                %Allow met uptake
                    model.ub(rxnIndx) = 1000;
                    disp([excMetabolite{1} ' added to the medium'])
                %end
            end
        end
    end
end
%Allow all secretions
for i=1:length(model.rxns)
    prods = find(model.S(:,i)>0);
    if isempty(prods)
        model.ub(i) = 1000;
    end
end
     
%Is the model growing now?
sol = solveLP(model)
priorValue = abs(sol.f);
%If overconstrained then look for the missing components for the model to
%grow
essential = [];
sensitivities = [];
if priorValue<gRate 
    for i=1:length(excRxnIndxs)
        rxnIndx = excRxnIndxs(i);
        if ~ismember(i,rxnIndx)
            temp_model = model;
            if temp_model.ub(rxnIndx) == 0
                %Get exchange met name
                %prods = model.metNames{find(model.S(:,rxnIndx)>0,1)};
                %disp(prods)
                disp(model.rxns{rxnIndx})
                %Allow exchange
                temp_model.ub(rxnIndx) = 1000;
                sol = solveLP(temp_model);
                growth = sol.x(GRindex);
                sensitivity = (growth-priorValue)/gRate;                
                if abs(sol.f)>priorValue
                    disp(['Essential component found: ' model.rxns{rxnIndx}])
                    essential = [essential;rxnIndx];
                    sensitivities = [sensitivities;sensitivity];
                end
            end
        end
    end
end
end