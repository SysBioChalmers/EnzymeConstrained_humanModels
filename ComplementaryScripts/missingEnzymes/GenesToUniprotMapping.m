%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [] = GenesToUniprotMapping(model)
%
% Receives the HMR model structure and the list of human genes in
% http://www.ensembl.org/biomart/martview/27e01d07c7284a13ed26b212c23e6d3dadds
% the downloaded file contains the columns: <Gene (ensembl) / Gene name /
% HGNC ID/ENSEMBL_UniprotKB/Swiss-Prot ID>. The algorithm searches for the
% ENSEMBL_Uniprot code associated to every gene of the model in the ENSEMBL 
% database.
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

%function ENSEMBL_Uniprot = GenesToUniprotMapping(model)

%    cd ../../Databases/ENSEMBL
%     file_name = 'mart_export.txt';
%     fID       = fopen(file_name);
%     data      = textscan(fID,'%s %s %s %s','delimiter',',');
%     fclose('all');
    
    %[ENSEMBL_Uniprot, NonUni_ENSEMBL] = findInEnsembl(data,model);
    %[NonUni_Swissprot, NonUni_kegg]   = findInSwissprot_KEGG(data,model);
    %[missingGenes,matchedGenes]       = missingUniprots(data,model);
    %missingGenes                      = findECnumbers(missingGenes);
    missingGenes                      = findMissingGenesIndxs(model,missingGenes);
        
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function missingGenes = findMissingGenesIndxs(model,missingGenes)
    missingGenes{4} = [];
    grRules         = model.grRules_ENSEMBL;
    for i=1:length(missingGenes{1})
        if ~isempty(missingGenes{3}{i})
           matching = indexes_string(grRules,missingGenes{1}(i),false);
           if ~isempty(matching)
               missingGenes{4}{i} = matching;
           end
        else 
            missingGenes{4}{i} = {};
        end
    end
    missingGenes{4} = transpose(missingGenes{4});
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receive the genes that are included in the original model
% but were not matched to any uniprot code after the EC-GEM extension and
% looks for an EC number for them in the Swissprot and KEGG databases.

function missingGenes = findECnumbers(missingGenes)
    %cd ..
    cd ../../Databases
    load('hsa_ProtDatabase.mat');
    Swissprot_Uni{1} = [swissprot(:,1)]; Swissprot_Uni{2} = [swissprot(:,4)];
    kegg_Uni{1}      = [kegg(:,1)];      kegg_Uni{2}      = [kegg(:,3)];
    missingGenes{3}  = [];
    
    
    for i=1:length(missingGenes{1})
        uniprot = missingGenes{2}{i};
        % If the gene has an associated Uniprot code then it looks in
        % Swissprot for its EC number.
        if ~isempty(uniprot)
           indxs = indexes_string(Swissprot_Uni{1},uniprot,true);  
           if ~isempty(indxs) && ~isempty(Swissprot_Uni{2}(indxs))
               missingGenes{3} = [missingGenes{3}; Swissprot_Uni{2}(indxs)];
           else
           % If it wasn't found in Swissprot then it is searched in the
           % KEGG database cell array.
               indxs = indexes_string(kegg_Uni{1},uniprot,true);
               if ~isempty(indxs) && ~isempty(kegg_Uni{2}(indxs))
                   missingGenes{3} = [missingGenes{3}; kegg_Uni{2}(indxs)];
               else
                   missingGenes{3} = [missingGenes{3}; {''}];
               end
           end
        else
            missingGenes{3} = [missingGenes{3}; {''}];
        end
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [missingGenes,matchedGenes] = missingUniprots(ENSEMBL,model)

    enzymes         = model.enzymes;
    matchedGenes{1} = []; matchedGenes{2} = [];
    for i=1:length(enzymes)
        % Looks for each model enzyme in the ENSEMBL cell
        indxs = indexes_string(ENSEMBL{4},enzymes{i},true); 
        if ~isempty(indxs)
            matchedGenes{1} = [matchedGenes{1}; ENSEMBL{1}(indxs)];
            matchedGenes{2} = [matchedGenes{2}; indxs];
        end  
    end
    missingGenes{1} = []; missingGenes{2} = [];
    for i=1:length(model.genes)
        % Identify the model genes that were not assigned with any uniprot
        % code
        if ~ismember(model.genes{i},matchedGenes{1})
            missingGenes{1} = [missingGenes{1}; model.genes(i)];
            % looks for the "missing gene" in the ENSEMBL data
            indxs           = indexes_string(ENSEMBL{1},model.genes{i},true);
            if ~isempty(indxs) && ~isempty(ENSEMBL{4}{indxs})
                %disp(ENSEMBL{4}(indxs))
                %missingGenes{2}{length(missingGenes{1})} = ENSEMBL{4}(indxs);
                missingGenes{2} = [missingGenes{2}; ENSEMBL{4}(indxs)];
            else
                %missingGenes{2}{length(missingGenes{1})} = [];
                missingGenes{2} = [missingGenes{2}; {''}];
            end
        end
    end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NonUni_Swissprot, NonUni_kegg] = findInSwissprot_KEGG(data,model)
%   Now each gene its translated to its correspondant gene short name and
%   searched in the Uniprot and KEGG databases for H. sapiens
    cd ..
    load('hsa_ProtDatabase.mat');
    Swissprot_Uni{1}    = [swissprot(:,1)]; 
    Swissprot_Uni{2}    = [swissprot(:,3)];
    NonUni_Swissprot{1} = []; 
    NonUni_Swissprot{2} = [];   
    kegg_Uni{1}         = [kegg(:,1)]; kegg_Uni{2} = [kegg(:,3)];
    NonUni_kegg{1}      = []; NonUni_kegg{2} = [];
    
    for i=1:length(model.genes)
        indx = indexes_string(data{1},model.genes{i},true);
        if ~isempty(indx)
            % Gets the translated gene name
            gene_name = data{2}(indx);
            indx      = indexes_string(Swissprot_Uni{2},gene_name,true);
            if isempty(Swissprot_Uni{2}(indx))
                NonUni_Swissprot{1} =[NonUni_Swissprot{1}; model.genes(i)];
                matching            = indexes_string(model.grRules_ENSEMBL,...
                                                     model.genes{i},false);
                NonUni_Swissprot{2} = [NonUni_Swissprot{2}; {matching}];        
            end
        end
        % Gets the translated gene name from ENSEMBL to HGNC
        gene_name      = data{3}{i};
        gene_name      = gene_name(strfind(gene_name,'HGNC')+5:end);
        % Looks for the HGNC in the KEGG database
        %indx           = indexes_string(kegg_Uni{2},gene_name,true);
        indx            = find(strcmpi(kegg_Uni{2},gene_name));
        if isempty(indx) %&& isempty(kegg_Uni{1}(indx))
            NonUni_kegg{1} =[NonUni_kegg{1}; model.genes(i)];
            matching       = indexes_string(model.grRules_ENSEMBL,...
                                            model.genes{i},false);
            NonUni_kegg{2} = [NonUni_kegg{2}; {matching}];        
        end    
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ENSEMBL_Uniprot, NonUni_ENSEMBL] = findInEnsembl(data,model)
    ENSEMBL_Uniprot{1}   = []; ENSEMBL_Uniprot{2}   = [];
    for i=1:length(model.genes)
        % Looks for each model gene in the ENSEMBL file
        indxs = indexes_string(data{1},model.genes{i},false);
        if ~isempty(indxs)
            matching = find(~cellfun(@isempty,data{4}(indxs)),1);
            % If the gene was found on the ENSEMBL list and it has an
            % associated ENSEMBL_Uniprot code then both (ENSEMBL_Uniprot 
            % and gene name) are saved as output
            if ~isempty(matching)
                ENSEMBL_Uniprot{1} = [ENSEMBL_Uniprot{1}; data{1}(indxs(matching))];
                ENSEMBL_Uniprot{2} = [ENSEMBL_Uniprot{2}; data{4}(indxs(matching))]; 
            end
           
        end
    end
    
    % Loop that looks for every gene in the model on the ENSEMBL_Uniprot
    % cell.
   NonUni_ENSEMBL{1} =[];NonUni_ENSEMBL{2} =[];
    for i=1:length(model.genes)
        if ~ismember(model.genes{i},ENSEMBL_Uniprot{1})
           NonUni_ENSEMBL{1} = [NonUni_ENSEMBL{1}; model.genes(i)];
           matching = indexes_string(model.grRules_ENSEMBL,model.genes{i},true);
           NonUni_ENSEMBL{2} = [NonUni_ENSEMBL{2}; {matching}];           
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag == false
        matching = find(~cellfun(@isempty,matching));
    else
        matching = find(~cellfun(@isempty,matching),1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      U_12=[]; U_13 = []; U_14=[];U_23=[]; U_24=[]; U_34=[];U_123=[];U_124=[];...
%      U_134=[]; U_234=[];U_1234=[];
%      nonUniGenes{1} = model.genes; nonUniGenes{2} = NonUni_ENSEMBL{1};
%      nonUniGenes{3} = NonUni_Swissprot{1}; nonUniGenes{4} = NonUni_kegg{1};
% 
%     for i=1:length(nonUniGenes{1})
% 
%         U_12  = union_2(nonUniGenes{1}(i),nonUniGenes{2},U_12);  
%         U_13  = union_2(nonUniGenes{1}(i),nonUniGenes{3},U_13);     
%         U_14  = union_2(nonUniGenes{1}(i),nonUniGenes{4},U_14);  
%         U_123 = union_3(nonUniGenes{1}(i),nonUniGenes{2},nonUniGenes{3},U_123);
%         U_124 = union_3(nonUniGenes{1}(i),nonUniGenes{2},nonUniGenes{4},U_124);
%         U_134 = union_3(nonUniGenes{1}(i),nonUniGenes{3},nonUniGenes{4},U_134);
%         
%         if (ismember(nonUniGenes{1}(i),nonUniGenes{2}) && ...
%             ismember(nonUniGenes{1}(i),nonUniGenes{3}) && ...
%             ismember(nonUniGenes{1}(i),nonUniGenes{4}))   
%             U_1234  = [U_1234; nonUniGenes{1}(i)];
%         end
%     end
%     
%     for i=1:length(nonUniGenes{2})
%         U_23  = union_2(nonUniGenes{2}(i),nonUniGenes{3},U_23);
%         U_24  = union_2(nonUniGenes{2}(i),nonUniGenes{4},U_24);
%         U_234 = union_3(nonUniGenes{2}(i),nonUniGenes{3},nonUniGenes{4},U_234);
%     end
%     
%     U_34  = union_2(nonUniGenes{2}(i),nonUniGenes{4},U_34);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function U_ij = union_2(cell_i,cellArray_j,U_ij)  
%     if ismember(cell_i,cellArray_j)
%         U_ij  = [U_ij; cell_i];
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function U_ijl = union_3(cell_i,cellArray_j,cellArray_l,U_ijl)  
%     if (ismember(cell_i,cellArray_j) && ...
%         ismember(cell_i,cellArray_l))   
%         U_ijl  = [U_ijl; cell_i];
%     end            
% end
