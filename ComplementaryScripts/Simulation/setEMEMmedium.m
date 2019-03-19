function [model,components] = setEMEMmedium(model,glucBound,oxBound,stdBound,irrev)
% setEMEMmedium
%
% Set a EMEM medium for the model, if the model is overconstrained
% (according to a provided experimental gRate measurement), then an
% iterative process that allows the uptake of new metabolites until the 
% model grows properly.
%
%   model       An HMR-based GEM
%   glucBound   UB for the glucose uptake rate [mmol/gDw h]
%   oxBound     UB for the oxygen uptake rate [mmol/gDw h]
%   stdBound    UB used as a standard uptake rate [mmol/gDw h]
%   irrev       Indicates if a model is on its irreversible format 
%               (default = true)

%   model        model 
%   components   Components that were successfully added to the medium
%
% Ivan Domenzain.      Last edited: 2019-03-14
%
exchangeMets  =  {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cystine';...
                  'Glutamate';'Glutamine';'Glycine';'Histidine';'Isoleucine';...
                  'Leucine';'Lysine';'Methionine';'Phenylalanine';'Proline';...
                  'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine';...
                  'Choline';'Folate';'Inositol';'Nicotinamide';'Pantothenate';...
                  'pyridoxine';'Riboflavin';'Thiamin';'Glucose';'o2';'h2o';...
                  'na+';'k+';'mg+';'PI';'Sulfate';'ca2+';'fe2+';'fe3+';'hco3-';...
                  'h+'};
if nargin<5
    irrev = true;
end
%getExchangeRxns works with the unconstrained field if present, lets remove
%it temporarily to avoid any inconsistency with the model bounds
 if isfield(model,'unconstrained')
     model = rmfield(model,'unconstrained');
 end
unconstrained = zeros(length(model.mets),1);
%Allow all secretions and block all uptakes
[~,exchIndxs]       = getExchangeRxns(model);
%model.ub(exchIndxs) = 1000;
if irrev
    upTakeIndxs = exchIndxs(find(contains(model.rxns(exchIndxs),'_REV')));
    model.ub(upTakeIndxs) = 0;
else
    model.lb(exchIndxs) = 0;
    upTakeIndxs         = exchIndxs;
end

components = [];

for i=1:length(exchangeMets)
    excMetabolite = exchangeMets(i);   
    if ~ismember(excMetabolite,model.metNames)
        for j=1:length(upTakeIndxs)
            %Then search the met exchange index (one by one)
            rxnIndx = upTakeIndxs(j);
            %Get the products compartment
            metIndex = find(model.S(:,rxnIndx));
            excMet   = model.metNames(metIndex);
            prodComp = model.metComps(metIndex); 
            if isempty(prodComp)
                prodComp = 10;
            end
            %If the searched metabolite matches with one of exchange rxn 
            %mets then find its uptake rxn and enable it 
            if strcmpi(model.compNames(prodComp),'extracellular') & strcmpi(excMet,excMetabolite)
                components              = [components; excMetabolite];
                unconstrained(metIndex) = 1;
                %Allow met uptake
                if strcmpi(excMetabolite,'glucose')
                    Bound = glucBound;
                else
                    if strcmpi(excMetabolite,'O2')||strcmpi(excMetabolite,'H2O')||strcmpi(excMetabolite,'H+')
                       Bound = oxBound; 
                    else
                       Bound = stdBound;
                    end
                end
                        
                if irrev
                    model.ub(rxnIndx) = Bound;
                else
                    model.lb(rxnIndx) = -Bound;
                end
            end
        end
    end
end
%Allow all secretions
%[~,excRxnIndxs]       = getExchangeRxns(model,'out');
%model.ub(excRxnIndxs) = 1000;

model.unconstrained = unconstrained;

%If DM_rxns are present (coming from RECON3D) block them
DMrxns = find(contains(model.rxns,'DM_'));
model.ub(DMrxns) = 0;
model.lb(DMrxns) = 0;
Sink_rxns = find(contains(model.rxns,'sink_'));
model.ub(Sink_rxns) = 0;
model.lb(Sink_rxns) = 0;
%Allow biomass production
bioIndex = find(contains(model.rxns,'humanGrowthOut'));
if ~isempty(bioIndex)
    model.ub(bioIndex) = 1000;
    model.lb(bioIndex) = 0;
end

end