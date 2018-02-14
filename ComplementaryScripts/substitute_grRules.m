%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = substitute_grRules(original_model)

% Receives the HMR model .mat structure and adds a new grRules field with 
% correspondent equivalences between ENSEMBL gene codes (original HMR) and 
% gene short names in order to provide compatibility with the Swissprot and
% KEGG protein databases.
%
% Ivan Domenzain. Last edited: 2017-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = substitute_grRules(model)
    current = pwd;
    %Equivalences table between ENSEMBL and HGNC gene IDs for the organism
    %previously downloaded from:   http://www.ensembl.org/biomart,
    %including the fields [Gene stable ID/Gene name/HGNC ID/ Swissprot ID]
    file_name = 'mart_export.txt';
    cd ../Databases/ENSEMBL
    fID       = fopen(file_name);
    data      = textscan(fID,'%s %s %s %s','delimiter',',');
    fclose('all');    
    model.grRules_ENSEMBL = model.grRules;
    
    for i=1:length(model.grRules)
        model.grRules(i) = strrep(model.grRules(i),' or',' OR ');
        if ~isempty(model.grRules{i})
            str_cells = strsplit(model.grRules{i},' ');
            grRule    = [];
            %Decomposes the grRule into individual genes
            for j=1:length(str_cells)
                str = (str_cells{j});
                if ~strcmpi(str,'or') && ~strcmpi(str,'and')
                    ensemblID = replace(str,'(','');
                    ensemblID = replace(ensemblID,')','');
                    index     = indexes_string(data{1},ensemblID);
                    %If the gene is found in the ENSEMBL iDs it is replaced
                    %by its corresponding gene name
                    if ~isempty(index)
                         gene         = data{2}(index);
                    end
                end
                grRule = strcat(grRule,str_cells{j});
                if j<length(str_cells)
                    grRule = strcat(grRule,{' '});
                end
            end
            model.grRules{i} = char(grRule);
        end
        disp(strcat('ready with grRule #',string(i)))
    end
    cd (current)
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str)
    matching  = strfind(cell_array,str);
    matching = find(~cellfun(@isempty,matching),1);
end
