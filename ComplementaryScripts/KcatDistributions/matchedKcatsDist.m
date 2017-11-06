%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function kcat_values = matchedKcatsDist(model_kcats)
% 
% Gets all the non repeated Kcat values matched to an enzyme constrained
% GEM for the forward and backward reactions. If the model contains the
% field subsystems then the distribution of values can be splitted into
% different metabolic blocks.
%
% Ivan Domenzain.  Last edited: 2017-11-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kcat_values = matchedKcatsDist(model,model_kcats, pathGroups)
    matrix      = model_kcats.forw.kcats; 
    for i=1:4
        kcat_values{i} = [];
    end
    flag        = isfield(model,'subSystems');
    kcat_values = getValues(matrix/3600,flag,pathGroups,kcat_values,model);
    matrix      = model_kcats.back.kcats; 
    kcat_values = getValues(matrix/3600,flag,pathGroups,kcat_values,model);
    for i=1:4
        kcat_values{i} = unique(kcat_values{i});
    end
end
    
function kcat_values = getValues(matrix,flag,pathGroups,kcat_values,model)  
    [m,n] = size(matrix);
    for i=1:m
        for j=1:n
            if matrix(i,j)~=0
                kcat_values{1} = [kcat_values{1}; matrix(i,j)];
                if flag
                    
                    for k=1:3
                        if any(strcmpi(model.subSystems(i),pathGroups{2}{k}))
                            kcat_values{k} = [kcat_values{k}; matrix(i,j)];
                        end
                    end
                end
            end
        end
    end
end