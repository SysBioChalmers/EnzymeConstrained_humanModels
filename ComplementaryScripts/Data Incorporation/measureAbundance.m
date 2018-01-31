%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,count] = measureAbundance(model,abundance_file)
% 
%
% Benjam??n J. S??nchez. Last edited: 2016-03-18
% Ivan Domenzain.      Last edited: 2018-01-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,count] = measureAbundanceHepG2(enzymes)
    cd ../../Databases/paxDB
    %Read downloaded data of abundance:
    fID        = fopen('prot_abundance.txt');
    data1       = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',12);
    fclose(fID);
    %Read equivalence table for ENSP to uniprot codes
    cd ../ENSEMBL
    fID       = fopen('ENSEMBL_2_UNIPROT.txt');
    data2     = textscan(fID,'%s %s %s','delimiter',',');
    fclose(fID);
    [genes,abundance] = translateProteinIDs(data1,data2);

    Load KEGG data:
    cd ..
    data      = load('hsa_ProtDatabase.mat');
    swissprot = data.swissprot;

    Main loop:
    MW_ave  = mean(cell2mat(swissprot(:,5)));
    Pmodel  = 0;
    Ptot    = 0;
    counter = zeros(size(genes));
    for i = 1:length(genes)
        MW = MW_ave;
        gene_name = genes{i};
        %Find gene in swissprot database:
        swiss = swissprot(:,1);
        index = find(strcmpi(gene_name,swiss));
        %for j = 1:length(swissprot)
        %    if ~isempty(strfind(swissprot{j,i},gene_name))
        if ~isempty(index)   
                MW = swissprot{index,5};
                %Check if uniprot is in model:
                for k = 1:length(enzymes)
                    if strcmpi(enzymes{k},swissprot{index,1})
                        counter(i) = 1;
                        Pmodel = Pmodel + MW*abundance(i)/1000000;
                    end
                end
            %end
        %end
        end
        Ptot = Ptot + MW*abundance(i)/1000000;
        disp(['Calculating total abundance: Ready with gene ' gene_name])
    end

    f     = Pmodel/Ptot;
    count = [length(counter);sum(counter(:,1))];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genes,abundance] = translateProteinIDs(Data1,Data2)

    genes     = [];
    abundance = [];

    for i=1:length(Data1{2})
        %Trim gene name:
        ENSP_name = Data1{2}{i};
        ENSP_name = ENSP_name(strfind(ENSP_name,'.')+1:end);
        % Look for the ENSP 
        index     = find(strcmpi(ENSP_name,Data2{2}),1);
        
        if ~isempty(index)
            % Look for the equivalent UNIPROT code
            if ~isempty(Data2{2}(index))
                genes     = [genes; Data2{3}(index)];
                abundance = [abundance; Data1{3}(i)];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




