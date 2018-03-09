%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
%  
% This function goes through to each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares the cumulative the cumulative distributions of both.
%
% Raphael Ferreira     Last edited: 2018-03-09
% Ivan Domenzain.      Last edited: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel)
	if nargin>2
		[model,ecModel] = block_production(model,ecModel);
	end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,ecModel] = block_production(model,ecModel,metName)
 	[~,Uptake_indxs] = findExcRxns(model);
 	Uptake_indxs 	 = find(Uptake_indxs);
 	% Gets the indexes of all the extracellular mets
 	compIndex		= find(strcmpi(model.compNames,'extracellular'));
 	compMets		= find(model.metComps==compIndex);
 	%For each metabolite uptake that is going to be blocked
 	for i=1:length(metName)
 		ExtMet_Index    = find(strcmpi(model.metNames(compMets),metName(i)));

 		if length(ExtMet_Index)>1
 			disp(['warning: multiple matches were found in the compartment for ' string(metName(i))])
		else
			uptakeIndex = find(model.S(ExtMet_Index,:));
			product_Indxfind(model.S(ExtMet_Index,uptakeIndex)>0)
			else

				model.lb(uptakeIndex) = 0;
			end
		end
 	end
end