% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [new_database] = modifyKEGGProtDB(database,name)
% 
% Ivan Domenzain. Last edited: 2017-09-20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
url  = 'http://rest.kegg.jp/list/hsa';
data = textscan(webread(url),'%s %s','delimiter','\t');
data = data{1};
cd ../Databases
fileID = fopen('KEGGgenes.txt','w');
for i=1:length(data)
    fprintf(fileID,'%s\n',data{i});
end
fclose(fileID);

% Now upload the generated file (KEGGgenes.txt) on 
% http://www.uniprot.org/uploadlists/ and select the option from KEGG to
% UniProtKB, generate and download a file containing the next columns
% Entry / KEGG ID / Gene name / Protein name / EC number / Isoform Map
% The generated file should be saved on the .../Databases folder with the
% name KEGG-uniprot_hsa.tab
data = [];
fileID = fopen('KEGG-uniprot_hsa.tab','r');
data   = textscan(fileID,'%s %s %s %s %s','delimiter','\t');
load('hsa_ProtDatabase.mat')
KEGGIDs = kegg(:,3);
for i=1:length(KEGGIDs)
    str  = strcat('hsa:', KEGGIDs{i});
    indx = find(strcmpi(str,data{2}));
    if length(indx)>1
        ecCells = data{3}(indx);
        ecCells = find(~cellfun(@isempty,ecCells));
        if ~isempty(find(~cellfun(@isempty,ecCells)))
            disp(ecCells)
        end
        
    end
%     if ~isempty(indx) && ~isempty(data{3}(indx))
%        kegg(i,3) =  data{3}(indx);
%     end
end
i = unique(data{2});