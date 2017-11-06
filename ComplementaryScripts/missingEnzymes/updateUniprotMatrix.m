%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model] = updateUniprotMatrix(model)
%
% Receives the HMR model structure 
%
% INPUT:
%   model    GEM matlab structure
%
% OUTPUT:
%   model    GEM matlab structure with the updated Uniprot codes matrix
%
% Ivan Domenzain.      Last edited: 2017-10-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function model = updateUniprotMatrix(model,model_data)
    
    cd ../Databases
    load('hsa_ProtDatabase.mat')
    file_name = 'mart_export.txt';
    fID       = fopen(file_name);
    data      = textscan(fID,'%s %s %s %s','delimiter',',');
    fclose('all');
    
    for i=1:length(model_data.uniprots)
        uniprotRow    = model_data.uniprots(i,:);
        grRule        = strsplit(model.grRules{i},' or ');
        nonEmptyIndxs = find(~cellfun(@isempty,uniprotRow));
        gRllength     = length(grRule);
        Unilength     = length(grRule);
        if isempty(grRule{1})
            gRllength = 0;
        end
        if Unilength < gRllength
            for j=1:length(nonEmptyIndxs)
                
        end
            
        
        %display(['i' i])
        %display(length(grRule))
        %pause
        for j=1:length(uniprotRow)
            if ~isempty(uniprotRow{j})
            end
        end
    end
%end
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag == true
        matching = find(~cellfun(@isempty,matching),1);
    else
        matching = find(~cellfun(@isempty,matching));
    end
end
    