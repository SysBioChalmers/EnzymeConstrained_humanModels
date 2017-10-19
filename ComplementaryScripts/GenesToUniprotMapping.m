%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [] = GenesToENSEMBL_UniprotMapping(model)
%
% Receives the HMR model structure and the list of human genes in
% http://www.ensembl.org/biomart/martview/27e01d07c7284a13ed26b212c23e6d3dadds
% the downloaded file contains the columns: <Gene (ensembl) / Gene name /
% HGNC ID / ENSEMBL_UniprotKB/Swiss-Prot ID>. The algorithm searches for the
% ENSEMBL_Uniprot code associated to every gene of the model in the ENSEMBL database.
%
% INPUT:
%   model           GEM matlab structure
%
% OUTPUT:
%   genesCounter    Cell containing (1) number of genes in ENSEMBL, (2) number
%                   of genes in the model, (3) number of genes in the model 
%                   with an associated ENSEMBL_Uniprot code.
%
%  Ivan Domenzain.      Last edited: 2017-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ENSEMBL_Uniprot = GenesToUniprotMapping(model)

    cd ../Databases
    file_name = 'mart_export.txt';
    fID       = fopen(file_name);
    data      = textscan(fID,'%s %s %s %s','delimiter',',');
    fclose('all');
    
    ENSEMBL_Uniprot{1}   = []; ENSEMBL_Uniprot{2}   = [];
    for i=1:length(model.genes)
        % Looks for each model gene in the ENSEMBL file
        indxs = indexes_string(data{1},model.genes{i});
        if ~isempty(indxs)
            matching = find(~cellfun(@isempty,data{4}(indxs)),1);
            % If the gene was found on the ENSEMBL list and it has an
            % associated ENSEMBL_Uniprot code then both (ENSEMBL_Uniprot and gene name) are
            % saved as output
            if ~isempty(matching)
                ENSEMBL_Uniprot{1} = [ENSEMBL_Uniprot{1}; data{1}(indxs(matching))];
                ENSEMBL_Uniprot{2} = [ENSEMBL_Uniprot{2}; data{4}(indxs(matching))]; 
            end
           
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.

function matching = indexes_string(cell_array,str)
    matching  = strfind(cell_array,str);
    matching = find(~cellfun(@isempty,matching));
end
