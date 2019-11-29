function exchModel = setHamsMedium(model,Csource,irrev,measuredMets,fluxes)
% setHamsMedium
%
% Set a Hams culture medium for a humanGEM based model. This function works
% for either standard or enzyme constrained-GEMs
%
%   model           An ihuman-based GEM
%   Csource         (string) metName for the main carbon source
%   irrev           (logical) TRUE if model comes in an irreversible format.
%                   Default = false
%   measuredMets    (cell) met IDs for measured compounds. Optional
%   fluxes          (vector) Measured fluxes [mmol/gDw h], always as positive 
%                   quantities. Optional
%
%   exchModel       (struct) Model with Ham's media constraints
%
% Ivan Domenzain.      Last edited: 2019-11-29

if nargin<5
	fluxes = [];
    if nargin<4
    	measuredMets = [];
        if nargin<3 
        	irrev = false;
        end
    end
end
%Remove unconstrained field, if available
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
end
%Check if boundary metabolites are present, if so then remove them
boundaryIndx = find(strcmpi(model.compNames,'Boundary'));
boundary     = find(model.metComps==boundaryIndx);
if ~isempty(boundary)
    model = removeMets(model,boundary,false,false,false,true);
end
%Locate carbon source uptake reaction
Cindxs   = find(strcmpi(model.metNames,Csource));
Cindxs   = Cindxs(find(model.metComps(Cindxs)==1));
Csource  = model.mets(Cindxs);
%Ham's media composition
mediaComps ={ Csource{1} ...
             'm01365s' ...	%arginine[s]
             'm02125s' ...	%histidine[s]
             'm02471s' ...	%methionine[s]
             'm02724s' ...	phenylalanine[s]
             'm03089s' ...	tryptophan[s]
             'm03101s' ...	tyrosine[s]
             'm01307s' ...	alanine[s]
             'm01986s' ...	glycine[s]
             'm02896s' ...	serine[s]
             'm02993s' ...	threonine[s]
             'm01370s' ...	aspartate[s]
             'm01974s' ...	glutamate[s]
             'm01369s' ...	asparagine[s]
             'm01975s' ...	glutamine[s]
             'm02184s' ...	isoleucine[s]
             'm02360s' ...	leucine[s]
             'm02770s' ...	proline[s]
             'm03135s' ...	valine[s]
             'm01628s' ...	cysteine[s]
             'm02982s' ...	thiamin[s]
             'm02159s' ...	hypoxanthine[s]
             'm01830s' ...	folate[s]
             'm01401s' ...	biotin[s]
             'm02680s' ...	pantothenate[s]
             'm01513s' ...	choline[s]
             'm02171s' ...	inositol[s]
             'm02583s' ...	nicotinamide[s]
             'm02817s' ...	pyridoxine[s]
             'm02842s' ...	riboflavin[s]
             'm02996s' ...  %thymidine[s]
             'm02630s' ... %Oxygen
             'm02040s' ... %Water
             'm02519s' ... %sodium
             'm02200s' ... %K+
             'm02751s' ... %Pi
             'm02946s' ... %Sulfate
             'm01413s' ... % calcium
             'm01821s' ... %Fe2+
             'm01822s' ... %Fe3+
             'm02046s' ... %HCO3-
             'm01450s' ... %cholesterol
             'm01596s' ... %CO2
             'm02039s'};   %H+	
%Default flux bounds
fluxBounds = ones(1,length(mediaComps))*1000;
%Check if provided mets are part of media's formulation
if ~isempty(measuredMets)
    [iA,iB] = ismember(measuredMets,mediaComps);
else 
    iA = 0;
end

if ~irrev
    %Modify fluxBounds with the correspondent provided flux measurements
    if sum(iA)>0
        fluxBounds(iB(iA)) = fluxes(iA);
    end
    %Set the right flux direction
    fluxBounds = -1*fluxBounds;
    %Allow all exchanges (Makes sure that all secretions are open
    [exchModel,~] = setExchangeBounds(model);
    %Close uptakes for mets not present in the media
    [exchModel,~] = setExchangeBounds(exchModel,mediaComps,fluxBounds,1000,true);
else
    if sum(iA)>0
        fluxBounds(iB(iA)) = fluxes(iA);
    end
    [exchRxns,exchIndxs] = getExchangeRxns(model);
    %Exclude protein pool exchange
    exchIndxs = exchIndxs(1:end-1);
    exchRxns  = exchRxns(1:end-1);
    %Differentiate between uptakes and production reactions
    uptkIndxs = exchIndxs(find(contains(exchRxns,'_REV')));
    prodIndxs = exchIndxs(find(~contains(exchRxns,'_REV')));
    %Open all production reactions
    exchModel = setParam(model,'ub',prodIndxs,1000);
    %close all uptakes
    exchModel = setParam(exchModel,'ub',uptkIndxs,0);
    %Open uptake of media components one by one
    for i=1:length(mediaComps)
        %Get metabolite indx
        metIndx = find(strcmpi(model.mets,mediaComps{i}));
        %Get rxns for metabolite
        metRxns = find(model.S(metIndx,:));
        %Get the uptake reaction for the metabolite
        metExchRxn = intersect(metRxns,uptkIndxs);
        %Open it!
        exchModel.ub(metExchRxn) = fluxBounds(i);
    end
end   
%Check if model is feasible
sol = solveLP(exchModel); 
if ~isempty(sol.x)
	disp('Constrained model is feasible')
else
    disp('Constrained model is unfeasible')
	exchModel = [];
end
end   