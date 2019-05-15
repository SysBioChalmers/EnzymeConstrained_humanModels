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
% Ivan Domenzain.      Last edited: 2019-03-21
%
exchangeMets  =  {'alanine';'arginine';'asparagine';'aspartate';'cysteine';...
                  'glutamate';'glutamine';'glycine';'histidine';'isoleucine';...
                  'leucine';'lysine';'methionine';'phenylalanine';'proline';...
                  'serine';'threonine';'tryptophan';'tyrosine';'valine';...
                  'choline';'folate';'inositol';'nicotinamide';'pantothenate';...
                  'pyridoxine';'riboflavin';'thiamin';'glucose';'O2';'H2O';...
                  'Na+';'K+';'Mg+';'Pi';'sulfate';'Ca2+';'Fe2+';'Fe3+';'HCO3-';...
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
%Allow all secretions and block all uptakes
[~,exchIndxs] = getExchangeRxns(model);
%Discard protein pool exchange pseudo-reaction
exchIndxs     = exchIndxs(find(~strcmpi(model.rxns(exchIndxs),'prot_pool_exchange')));
model.ub(exchIndxs) = 1000;
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
    index         = find(strcmp(model.metNames,excMetabolite));
    if ~isempty(index)
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

model.unconstrained = unconstrained;
%If DM_rxns are present (coming from RECON3D) block them
DMrxns              = find(contains(model.rxns,'DM_'));
model.ub(DMrxns)    = 0;
model.lb(DMrxns)    = 0;
Sink_rxns           = find(contains(model.rxns,'sink_'));
model.ub(Sink_rxns) = 0;
model.lb(Sink_rxns) = 0;
%Allow biomass production
bioIndex = find(contains(model.rxns,'humanGrowthOut'));
if ~isempty(bioIndex)
    model.ub(bioIndex) = 1000;
    model.lb(bioIndex) = 0;
end
%Block bicarbonate production to allow CO2 production
index           = find(strcmpi(model.rxns,'HMR_9078'));
model.ub(index) = 0;
end