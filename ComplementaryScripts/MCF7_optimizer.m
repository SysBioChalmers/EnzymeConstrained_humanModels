%function [resultsMat,targets,WT_yields] = MCF7_optimizer(model)
%
%
% Created.  Ivan Domenzain 2018-09-18
%
function [resultsMat,WT_yields] = MCF7_optimizer(model,action)

%Identify target metabolites and the rxns that produce them on the desired
%compartment
[~,AcCproduction] = getMetProdIndexes(model,'Acetyl-CoA',true,'cytosol');
[~,MalProduction] = getMetProdIndexes(model,'malonyl-CoA',true,'cytosol');
[~,CoAproduction] = getMetProdIndexes(model,'CoA',true,'cytosol');

%preallocate results matrices
m = length(model.genes)+1;
n = 3;
resultsMat = zeros(m,n+1);

%Find the glucose uptake rate
glucUptkIndx = find(strcmpi(model.rxns,'HMR_9034'));
%Get Wild type solution 
gIndex = find(strcmpi(model.rxnNames,'humanGrowthOut'));
% if nargin>2
%     model.ub(gIndex) = gRate;
%     model.lb(gIndex) = 0.99*gRate;
%     model.c(:)       = 0;
%     model.c(glucUptkIndx) = -1;
% end

base_sol   = solveLP(model,1);
gRate      = base_sol.x(gIndex);
%Calculate WT yield for selected metabolites
AcCprod    = sum(base_sol.x(AcCproduction));
MalProd    = sum(base_sol.x(MalProduction));
CoAprod    = sum(base_sol.x(CoAproduction));
glucUptake = base_sol.x(glucUptkIndx);
WTgrowth   = base_sol.x(gIndex);
indexes    = [{AcCproduction};{MalProduction};{CoAproduction};{glucUptkIndx};{gIndex}];

WT_AcCYield = AcCprod/glucUptake;
WT_MalYield = MalProd/glucUptake;
WT_CoAyield = CoAprod/glucUptake;
WT_yields   = [WT_AcCYield;WT_MalYield;WT_CoAyield];
disp(WT_yields)
growthYield = WTgrowth/(glucUptake*0.180);
disp(growthYield)

%Loop through all the original genes
for j=1:length(model.genes)    
    gene = model.genes(j);
    %Get met production yields for every  mutant
    [Fchanges,successWT,mutGrowth]  = getResultsForGene(model,model,gene,'pFBA',WT_yields,indexes,action);
    %save results
    gRateFC         = mutGrowth/WTgrowth;
    resultsMat(j,:) = [Fchanges,gRateFC];%Fchanges;
    disp(['Ready with gene ' num2str(j) ' feasible: ' num2str(successWT) ' FC: ' num2str(sum(Fchanges))])
end
%Rearrange results WT yields in the first row and each mutant in the rest
%of the rows
resultsMat(2:end,:)    = resultsMat(1:end-1,:);
[Fchanges,~,mutGrowth] = getResultsForGene(model,model,'','pFBA',WT_yields,indexes,action);
gRateFC                = mutGrowth/WTgrowth;
resultsMat(1,:)        = [Fchanges,gRateFC];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results,success,gRate] = getResultsForGene(strainModel,model,gene,method,WT_yields,indexes,action)

success = 0;
WT_AcCYield = WT_yields(1);
WT_MalYield = WT_yields(2);
WT_CoAyield = WT_yields(3);


AcCproduction = indexes{1}; 
MalProduction = indexes{2};
CoAproduction = indexes{3};
glucUptkIndx  = indexes{4};
gIndex        = indexes{5};

results = zeros(1,3); 
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
    [mutSolution,flag] = solveMutant(mutant,model,method);
    
    %if the optimization was feasible, then get results
    if flag ==1 
        %glucose uptake rate
        mutUptake = mutSolution(glucUptkIndx);
        if mutUptake>0
            success = 1;
            %get yield fold-changes for each precursor
            AcCyieldFC = getYieldFoldChange(AcCproduction,mutSolution,mutUptake,WT_AcCYield);
            MalYieldFC = getYieldFoldChange(MalProduction,mutSolution,mutUptake,WT_MalYield);
            CoAyieldFC = getYieldFoldChange(CoAproduction,mutSolution,mutUptake,WT_CoAyield);
            %save results
            results(1,1) = AcCyieldFC;
            results(1,2) = MalYieldFC;
            results(1,3) = CoAyieldFC;
            gRate = mutSolution(gIndex);
        end
    end
else
    [mutSolution,~] = solveMutant(model,model,method);
    results = ones(1,3);
    gRate   = mutSolution(gIndex);
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
function [mutSolution,flag] = solveMutant(mutant,model,method)
if strcmpi(method,'MOMA')
    [mutSolution,~, flag]=qMOMA(mutant,model,1.01);
    disp(flag)
elseif strcmpi(method,'pFBA')
    mutSolution = solveLP(mutant,1);
    if ~isempty(mutSolution.x) && any(mutSolution.x)
        mutSolution = mutSolution.x;
        flag = 1;
    else 
        mutSolution = [];
        flag = 0;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [metindx,metProduction] = getMetProdIndexes(model,metName,compFlag,compName)
metProduction  = [];
metindx        = find(strcmpi(model.metNames,metName));
%metindx        = metindx(find(model.metComps(metindx)==compIndex));
if ~compFlag
    for i=1:length(metindx)
        row = model.S(metindx(i),:);
        %metProduction  = find(row>0);
        %disp(model.rxns(find(row>0)))
        metProduction  = [metProduction,find(row>0)];
    end
else
    compIndex      = find(strcmpi(model.compNames,compName));
    metindx        = metindx(find(model.metComps(metindx)==compIndex));
    row            = model.S(metindx,:);
    metProduction  = find(row>0);
end
end


            


