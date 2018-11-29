function convertDatasetIDs(cellLine)
    %Open dataset
    cd (['../../models/' cellLine '/Data']) 
    fileID = fopen(['proteomics_TCGA_' cellLine '.txt']);
    data   = textscan(fileID,'%s','delimiter','\n','HeaderLines',1);
    data   = data{1};
    fclose('all');
    %Count the number of columns and reload the data in the properly
    %allocated format
    cellStr = strsplit(data{1},'\t');
    [~,nCols] = size(cellStr);
    format  = ['%s ' repmat('%f ',1,nCols-1)];
    format  = format(1:end-1);
    fileID  = fopen(['proteomics_TCGA_' cellLine '.txt']);
    data    = textscan(fileID,format,'delimiter','\t','HeaderLines',1);
    fclose('all');
    proteins = data{1};
    proteins = strrep(proteins,'"','');
    %Load swissprot data
    load('../../../Databases/swissprot_shortNames.mat')
    load('../../../Databases/ProtDatabase.mat')
    shortNames = swissprot_shortNames(:,3);
    shortNames = strrep(shortNames,'"','');
    for i=1:length(proteins)
        index = find(contains(shortNames,proteins{i}));
        L = length(index);
        if L>1 
            disp(L)
        end
    end
end