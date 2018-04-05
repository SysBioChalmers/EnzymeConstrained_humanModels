%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [grRules,grRules_ENSEMBL] = substitute_grRules(model)

% Receives the HMR model .mat structure and adds a new grRules field with 
% correspondent equivalences between ENSEMBL gene codes (original HMR) and 
% gene short names in order to provide compatibility with the Swissprot and
% KEGG protein databases.
%
% Ivan Domenzain. Last edited: 2018-04-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified_model = substitute_grRules(model)
    current = pwd;
    %Equivalences table between ENSEMBL and HGNC gene IDs for the organism
    %previously downloaded from:   http://www.ensembl.org/biomart,
    %including the fields [Gene stable ID/Gene name/HGNC ID/ Swissprot ID]
    file_name = 'mart_export.txt';
    cd ../Databases/ENSEMBL
    fID       = fopen(file_name);
    data      = textscan(fID,'%s %s %s %s','delimiter',',');
    fclose('all');    
    modified_model  = model;
    grRules_ENSEMBL = model.grRules;
    grRules         = grRules_ENSEMBL;
    for i=1:length(grRules)
        if ~isempty(grRules{i})
            str_cells = strsplit(grRules{i},' ');
            grRule    = [];
            %Decomposes the grRule into individual genes
            for j=1:length(str_cells)
                str = (str_cells{j});
                if ~strcmpi(str,'or') && ~strcmpi(str,'and')
                    ensemblID = replace(str,'(','');
                    ensemblID = replace(ensemblID,')','');
                    index     = find(contains(data{1},ensemblID),1);
                    %If the gene is found in the ENSEMBL iDs it is replaced
                    %by its corresponding gene name
                    if ~isempty(index)
                         gene = data{2}(index);
                         str  = gene;
                    end
                end
                grRule = strcat(grRule,str);
                if j<length(str_cells)
                    grRule = strcat(grRule,{' '});
                end
            end
            grRules{i} = char(grRule);
        end
        disp(strcat('ready with grRule #',num2str(i)))
    end
    modified_model.grRules = grRules;
    [grRules,~]            = standardizeGrRules(modified_model);
    modified_model.grRules = grRules;
    cd (current)
end    
