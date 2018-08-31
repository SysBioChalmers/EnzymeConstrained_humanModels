function [model,essential] = find_EMEM_ExchangeRxns(model)
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';...
                  'glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'nicotinamide';...
                  'Pantothenate';'Pyridoxine';'Riboflavin';'Thiamin';...
                  'Glucose';'O2';'H2O';'Na+';'K+';'Mg+';...
                  'PI';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
                  'H+';'cholesterol'};
              
[excRxnIDs,excRxnIndxs] = getExchangeRxns(model);
GRindex = find(strcmpi(model.rxns,'humanGrowthOut')); 
%ProtIndex = find(strcmpi(model.rxns,'prot_pool_exchange'));
ProtIndex = find(contains(model.rxnNames,'prot_'));
CO2Index = find(strcmpi(model.rxns,'HMR_9058'));
excRxnIndxs = setdiff(excRxnIndxs,GRindex);
excRxnIndxs = setdiff(excRxnIndxs,ProtIndex);
excRxnIndxs = setdiff(excRxnIndxs,CO2Index);
sol = solveLP(model)
%block uptakes and production of exchanged metabolites
model.ub(excRxnIndxs) = 0;
model.lb(excRxnIndxs) = 0;

sol = solveLP(model)
mediumComponents = [];
for i=1:length(exchangeMets)
    excMetabolite = exchangeMets(i);
    
    if ~isempty(find(strcmpi(model.metNames,excMetabolite)))
        
        for j=1:length(excRxnIndxs)
            %disp(excMetabolite)
            rxnIndx = excRxnIndxs(j);
            
            rxnMets = model.metNames(find(model.S(:,rxnIndx)));
            subs = model.metNames(find(model.S(:,rxnIndx)<0));
            subComp=model.metComps(find(model.S(:,rxnIndx)<0));
            prods = model.metNames(find(model.S(:,rxnIndx)>0));
            prodComp=model.metComps(find(model.S(:,rxnIndx)>0));
            if ~isempty(find(strcmpi(rxnMets,excMetabolite), 1))
                mediumComponents = [mediumComponents; i];
                %if contains(model.rxnNames(rxnIndx),'reversible')    
                    model.ub(rxnIndx) = 1000;
                    
                    %disp(excMetabolite)
%                     if ~isempty(subs)
%                         disp(subs)
%                         disp(subComp)
%                     else
%                         disp('No subs')
%                     end
%                     
%                     if ~isempty(prods)
%                         disp(model.rxns(rxnIndx))
%                         disp(prodComp)
%                         disp(model.lb(rxnIndx))
%                         disp(prods)
%                     else
%                         disp('No prods')
%                         
%                     end
                %end
            end
        end
    end
end
sol = solveLP(model)
priorValue = abs(sol.f);
if priorValue<0.4
    essential = [];
    for i=1:length(excRxnIndxs)
        temp_model = model;
        rxnIndx = excRxnIndxs(i);
        disp(model.metNames(find(model.S(:,rxnIndx))))
        if temp_model.ub(rxnIndx) == 0
            
          temp_model.ub(rxnIndx) = 1000;
          sol = solveLP(temp_model);
          if abs(sol.f)>1.01*priorValue
              disp(sol.f)
              disp(model.rxns(rxnIndx))
              essential = [essential;rxnIndx];
          end
        end
    end
end
end