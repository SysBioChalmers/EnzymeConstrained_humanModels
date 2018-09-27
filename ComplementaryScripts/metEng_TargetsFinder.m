function [resultsMat,WT_yields] = metEng_TargetsFinder(model,c_sourceID,metList,direction,compartment,action)
% metEng_TargetsFinder
%
% Function that computes the effect of all the possible gene deletions or 
% overexpressions in the total production or consumption yields for a given
% set of metabolites in the network.
% 
% model       a model structure
% c_sourceID  Rxn ID for the carbon source uptake reaction
% metList     Cell array containing the metName for all the metabolites to
%             analyze.
% direction   String 'production' when the total flux that produces a given
%             metabolite is required, 'consumption' when the total flux that 
%             consumes a metabolite is required instead.
% compartment String that indicates the model compartment in which the total 
%             flux calculations should be done
% action      ?deletion?for deleting all the genes in the model (one at a
%             time) and 'OE' for overexpressions with a default OE factor of 20.
%
% resultsMat  Matrix with the same number of rows as genes in the model
%             and each column corresponding to a metabolite in metList, the 
%             last column corresponds to changes in the growth rate. The 
%             values in this matrix are the fold changes in the total calculated 
%             flux yield (in relation to glucose consumption) for the metabolite 
%             of interest when modifying the i-th gene compared to a
%             wild-type flux distribution.
% WT_yields   Vector containing the value for the total flux yield for the
%             specified metabolites in [mol metabolite/mol glucose] and 
%             also for the growth rate in [g biomass/g glucose]
%
% Created.  Ivan Domenzain 2018-09-27
%

nMets = length(metList);
%Identify target metabolites and the rxns that produce them on the desired
%compartment
metIndexes = cell(0,nMets);
for i=1:nMets
    [~,metIndexes{i}] = getMetProdIndexes(model,metList{i},true,compartment,direction);
    disp([metList{i} ' ' direction ' reactions:'])
    disp(model.rxns(metIndexes{i}))
end

%preallocate results matrices
m = length(model.genes)+1;
resultsMat = zeros(m,nMets+1);

%Find the glucose uptake rate
glucUptkIndx = find(strcmpi(model.rxns,c_sourceID));
%Get Wild type solution 
objIndex   = find(model.c);
base_sol   = solveLP(model,1);
glucUptake = base_sol.x(glucUptkIndx);
WTgrowth   = base_sol.x(objIndex);
%Get all protein usage indexes
prot_Indxs = find(contains(model.rxnNames,'prot'));
%Exclude the protein pool
prot_Indxs = prot_Indxs(1:end-1);

%Calculate WT yield for selected metabolites
indexes = {};
WT_yields = zeros(1,nMets);
for i=1:nMets
    %get the total flux towards/from the i-th metabolite
    totalFlux    = sum(base_sol.x(metIndexes{i}));
    indexes      = [indexes;metIndexes(i)];
    WT_yields(i) = totalFlux/glucUptake;
    disp(['WT ' metList{i} ' ' direction ' yield: ' num2str(WT_yields(i)) ' [mol/mol glucose]'])
end
WT_yields(WT_yields==0) = 1E-6;
indexes = [indexes;{glucUptkIndx};{objIndex}];
growthYield = WTgrowth/(glucUptake*0.180);
disp(['WT Growth yield: ' num2str(growthYield) ' [gBiomass/gGlucose]'])
fprintf('\n')
%Loop through all the original genes
for j=1:length(model.genes)    
    gene = model.genes(j);
    %Get met production yields for every  mutant
    [Fchanges,successWT,mutGrowth] = getResultsForGene(model,model,gene,'pFBA',WT_yields,indexes,action,prot_Indxs);
    %save results
    gRateFC         = mutGrowth/WTgrowth;
    resultsMat(j,:) = [Fchanges,gRateFC];
    disp(['Ready with gene ' num2str(j) ' feasible: ' num2str(successWT) ' average FC: ' num2str(sum(Fchanges)/nMets)])
end
%Rearrange results WT yields in the first row and each mutant in the rest
%of the rows
resultsMat(2:end,:)    = resultsMat(1:end-1,:);
[Fchanges,~,mutGrowth] = getResultsForGene(model,model,'','pFBA',WT_yields,indexes,action,prot_Indxs);
gRateFC                = mutGrowth/WTgrowth;
resultsMat(1,:)        = [WT_yields,growthYield];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results,success,gRate] = getResultsForGene(strainModel,model,gene,method,WT_yields,indexes,action,prot_Indxs)
success = 0;
nMets   = length(WT_yields);
objIndex  = indexes{end};
glucUptkIndx  = indexes{end-1};

results = zeros(1,nMets); 
gRate   = 0;
%find gene in model's genes
if ~isempty(find(strcmpi(strainModel.genes,gene), 1)) || strcmpi(gene,'') 
    if strcmpi(action,'deletion')
        %Get mutant (single deletion)
         mutant = removeGenes(model,gene);     
    elseif strcmpi(action,'OE')
        %Get mutant (single OE)
        mutant = getMutant(strainModel,{gene,1,20});
    end
    %optimize
    [mutSolution,flag] = solveMutant(mutant,model,method,prot_Indxs);
    
    %if the optimization was feasible, then get results
    if flag ==1 
        %glucose uptake rate
        mutUptake = mutSolution(glucUptkIndx);
        if mutUptake>0
            success = 1;
            %get yield fold-changes for each precursor
            for i=1:nMets
                results(1,i) = getYieldFoldChange(indexes{i},mutSolution,mutUptake,WT_yields(i));
            end
            gRate = mutSolution(objIndex);
        end
    end
else
    [mutSolution,~] = solveMutant(model,model,method,prot_Indxs);
    results = ones(1,nMets);
    gRate   = mutSolution(objIndex);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YieldFC = getYieldFoldChange(prodIndex,mutSolution,Cuptake,WT_Yield)
mutProdYield = sum(mutSolution(prodIndex))/Cuptake;
%Get metabolite production yields ratio betweeen WT and
%double mutants
YieldFC = mutProdYield/WT_Yield;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mutSolution,flag] = solveMutant(mutant,model,method,prot_Indxs)
if strcmpi(method,'MOMA')
    [mutSolution,~, flag]=qMOMA(mutant,model,1.01);
    disp(flag)
elseif strcmpi(method,'pFBA')
    mutSolution = solveLP(mutant);
    if ~isempty(mutSolution.x) && any(mutSolution.x)
        %Get a first optimization
        mutSolution = mutSolution.x;
        index = find(mutant.c);
        %Fix optimal value for the objective and then minimize the total
        %sum of protein usages
        mutant.lb(index) = 0.999*mutSolution(index);
        mutant.ub(index) = mutSolution(index);
        mutant.c(:) = 0;
        mutant.c(prot_Indxs) = -1;
        mutSolution = solveLP(mutant,1);
        mutSolution = mutSolution.x;
        flag = 1;
    else 
        mutSolution = [];
        flag = 0;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [metindx,metProduction] = getMetProdIndexes(model,metName,compFlag,compName,direction)
metProduction  = [];
metindx        = find(strcmpi(model.metNames,metName));
%metindx        = metindx(find(model.metComps(metindx)==compIndex));
if strcmpi(direction,'production')
    coeff = 1;
elseif strcmpi(direction,'consumption')
    coeff = -1;
end

if ~compFlag
    for i=1:length(metindx)
        row = model.S(metindx(i),:);
        metProduction  = [metProduction,find(coeff*row>0)];
    end
else
    compIndex      = find(strcmpi(model.compNames,compName));
    metindx        = metindx(find(model.metComps(metindx)==compIndex));
    row            = model.S(metindx,:);
    metProduction  = find(coeff*row>0);
end
end


            


