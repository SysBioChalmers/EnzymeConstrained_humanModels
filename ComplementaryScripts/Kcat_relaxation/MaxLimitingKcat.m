%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [model,] = modifyKcats_sce(model,Protdatabase,GECKO_path)
%
% Function that gets the limiting Kcat values in an EC model (according to
% a sensitivity analysis) according to a sensitiviy analysis based on a
% lower limit of 0.01 for the Growth control coefficients.
%
% Ivan Domenzain    Last edited. 2017-12-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [limKcats,limRxns,breakFlag] = MaxLimitingKcat(model,prevUnicodes)

    %First, find the growth limiting Kcats in the Batch model
    [limKcats,limRxns] = findLimitingKcat(model);
    if ~isempty(limKcats{1})
        %Sort the limiting Kcats cell according to their control coefficient 
        limKcats = sortByCCoeff(limKcats);
        j        = 1;
        pos      = 0;
        % Choose the most controlling Kcat that is not included in the cell
        % prevUnicodes for avoiding infinite loops
        breakFlag = false;
        while j<=length(limKcats{1})
            %The top limiting coefficient is identified with its corresponding
            %uniprot code and the rxn index in the model because there might 
            %be promiscous enzymes can diferent reactions at different rates
            Protname = [limKcats{1}{j}(6:end) '_' num2str(limKcats{3}(j))];
            if ~ismember(Protname,prevUnicodes)
                pos = j;
                break
            else
                j = j+1;
            end
        end

        for i=1:length(limKcats)
            if pos ~=0
                limKcats{i} = limKcats{1,i}(pos);
            else
                % If all the significant Limiting Kcats have been previously
                % modified then the iterations are stopped
                limKcats{i} = [];
                breakFlag   = true;
            end
        end
    end
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kcat_output = sortByCCoeff(kcat_input)
    len = length(kcat_input);
    for i=1:len
        kcat_output{i} = [];
    end
    i = 0;
    
    while i<length(kcat_input{5})
        [M, indxs] = max(kcat_input{5});
        if ~isempty(indxs)
            for j=1:length(indxs)
                for k=1:len
                    kcat_output{k} = [kcat_output{k}; kcat_input{k}(indxs(j))];
                end
            end
        end
        kcat_input{5}(indxs) = 0;
        i = i+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kcatIndxs, limRxns] = findLimitingKcat(model)
   %Get an initial simulation optimizing for biomass production
    base_sol      = solveLP(model);

    %Extract the metabolite indexes corresponding to the model enzymes
    %(avoiding the protein pool)
    enz_met_indxs = find(~cellfun(@isempty,strfind(model.metNames,'prot_')));
    enz_met_indxs = enz_met_indxs(1:end-1);
    
    %Extract the enzyme usage pseudoreactions indexes
    enzUsageIndxs = find(~cellfun(@isempty,strfind(model.rxns,...
                                                           'prot_'))); 
    enzUsageIndxs = enzUsageIndxs(1:end-1);                                             
       
    kcatIndxs = cell(1,6);
    limRxns   = cell(1,3);     
    
    if base_sol.f ~= 0     
       %Get the original rxn indexes of the flux carrier enzymes (enzyme usages)   
        enzUsageFlux  = base_sol.x(enzUsageIndxs);
        fluxCarriers  = find(enzUsageFlux);
        fluxCarriers  = enzUsageIndxs(fluxCarriers);
       %Extract the metabolite index for the flux carrying enzymes
        usedEnzIndxs  = zeros(size(fluxCarriers));
        for i=1:length(usedEnzIndxs)
            usedEnzIndxs(i) = find(model.S(:,fluxCarriers(i))==1);
        end    
       %For each flux carrying enzyme calculate the growthRate control
       %coefficient (gRCC) on each of its kinetic coefficients
        for i = 1:length(usedEnzIndxs)
            enzPos        = usedEnzIndxs(i);
            enzName       = model.mets(enzPos);
           %Get the nonzero kinetic coefficients
            enz_rxnsCoeffs = find(model.S(enzPos,:));
           %For each enzymatic reaction (avoid enzyme usage pseudoRxn)
            if ~isempty(enz_rxnsCoeffs)
                for j=1:length(enz_rxnsCoeffs)-1
                    temp_model = model;
                    coeffPos   = enz_rxnsCoeffs(j);
                    temp_model.S(enzPos,coeffPos) = ...
                                       model.S(enzPos,coeffPos)/1000;
                    new_sol   = solveLP(temp_model);
                    gRCC      = (new_sol.f - base_sol.f)/base_sol.f;
                    
                    if abs(gRCC) > 1e-2
                        Kcat         = (-1/model.S(enzPos,coeffPos))/3600;
                        kcatIndxs{1} = [kcatIndxs{1}; enzName];
                        kcatIndxs{2} = [kcatIndxs{2}; enzPos];
                        kcatIndxs{3} = [kcatIndxs{3}; coeffPos];
                        kcatIndxs{4} = [kcatIndxs{4}; Kcat];
                        kcatIndxs{5} = [kcatIndxs{5}; gRCC];
                        kcatIndxs{6} = [kcatIndxs{6}; model.rxnNames(coeffPos)];
                    end
                end
            end
        end
    else
        for i=1:(enzUsageIndxs(1)-1)
            disp(['Relaxing Kcat coefficients for Rxn#' string(i)])
            temp_model = model;
            rxnPos     = enzUsageIndxs(i);
            temp_model.S(enz_met_indxs(1:end-1),rxnPos) = ...
                          temp_model.S(enz_met_indxs(1:end-1),rxnPos)/1000;  
            new_sol    = solveLP(temp_model);
            deltaGR    = abs(new_sol.f - base_sol.f);
            if deltaGR>0
                limRxns{1} = [limRxns{1}; model.rxnNames(rxnPos)];
                limRxns{2} = [limRxns{2}; rxnPos];
                limRxns{3} = [limRxns{3}; deltaGR];
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag
        matching = find(~cellfun(@isempty,matching),1);
    else
        matching = find(~cellfun(@isempty,matching));
    end
end


