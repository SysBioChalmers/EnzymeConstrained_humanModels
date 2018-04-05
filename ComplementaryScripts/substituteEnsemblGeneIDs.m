%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = substituteEnsemblGeneIDs(original_model)

% Receives the HMR model .mat structure and adds a new grRules field with 
% correspondent equivalences between ENSEMBL gene codes (original HMR) and 
% gene short names in order to provide compatibility with the Swissprot and
% KEGG protein databases.
%
% Ivan Domenzain. Last edited: 2017-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = substituteEnsemblGeneIDs(model)
    %Equivalences table between ENSEMBL and HGNC gene IDs for the organism
    %previously downloaded from:   http://www.ensembl.org/biomart,
    %including the fields [Gene stable ID/Gene name/HGNC ID/ Swissprot ID]
    file_name = 'mart_export.txt';
    cd ../Databases/ENSEMBL
    fID       = fopen(file_name);
    data      = textscan(fID,'%s %s %s %s','delimiter',',');
    fclose('all');    
    model.genes_ENSEMBL = model.genes;

    for i=1:length(model.genes)
        index     = find(strcmpi(data{1},model.genes(i)),1);
        %If the gene is found in the ENSEMBL iDs it is replaced
        %by its corresponding gene name
        if ~isempty(index)
            model.genes(i) = data{2}(index);
            
        end
        disp(strcat('ready with gene #',num2str(i)))
    end
 end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
