%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateDatabases
% Updates all databases for protein matching (KEGG and Swiss-Prot).
%
% Note: Before using this script, one should manually download from 
%       http://www.uniprot.org/uniprot a tab delimited file for the
%       desired organism with the following format:
%       Entry - Protein names - Gene names - EC number - Sequence
%       OBS: filter with the Swiss-Prot option
% 
% Benjam?n S?nchez. Last edited: 2017-04-18
% Ivan Domenzain.   Last edited: 2017-10-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateDatabases(GECKO_path,DB_path,organism_code)
    %File name of the Uniprot organism data:
    cd ([DB_path '/UNIPROT']) 
    fileID_uni     = fopen('uniprot_hsa.tab');
    swissprot      = textscan(fileID_uni,'%s %s %s %s %s %s',...
                                 'delimiter','\t');
    %Looks just for the entries with associated genes
    nonEmptyGenes  = find(~cellfun(@isempty,swissprot{3}));
    nonEmptyGenes  = find(~cellfun(@isempty,swissprot{3}));
    swissprot{4}   = strrep(swissprot{4},'"','');
    swissprot      = [swissprot{1}(nonEmptyGenes) swissprot{2}(nonEmptyGenes)...
                      swissprot{3}(nonEmptyGenes) swissprot{4}(nonEmptyGenes)...
                      swissprot{5}(nonEmptyGenes)];

    swissprot(1,:) = [];
    fclose(fileID_uni);
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
    for i = 1:length(swissprot)
        %Leave protein name as lower case, remove ';' from ECs & calculate MW:
        prot_name      = lower(swissprot{i,2});
        uni            = swissprot{i,1};
        sequence       = swissprot{i,5};
        MW             = calculateMW(sequence);
        swissprot{i,2} = prot_name;
        swissprot{i,4} = strrep(swissprot{i,4},';','');
        swissprot{i,5} = MW;
        swissprot{i,6} = sequence;
        disp(['Updating Swiss-Prot database: Ready with protein ' uni])
    end

    %Retrieve KEGG info (uniprot code - protein name - systematic gene name...
    %                    - EC number - MW - pathway - sequence):
    cd ([DB_path '/KEGG/' organism_code])
    file_names      = dir();
    file_names(1:2) = [];
    %wouldn?t be better to allocate kegg structure to the dimension of file_names 
    kegg            = cell(100000,7);
    n               = 0;
    for i = 1:length(file_names)
        %i=1;
        file_name = file_names(i).name;
        %3rd column: systematic gene name
        gene_name = file_name(1:end-4);
        %Retrieve all data as a cell with all rows:
        fID  = fopen(file_name);%
        text = textscan(fID,'%s','delimiter','\t');
        fclose(fID);
        text = text{1};
        cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])

        uni       = '';
        EC_names  = '';
        sequence  = '';
        MW        = 0;
        pathway   = '';
        for j = 1:length(text)
            %j=6;
            line = text{j};
            if length(line) > 10
                %1st column: uniprot code
                if strcmp(line(1:8),'UniProt:')
                    uni = line(10:end);

                %2nd column: protein name 
                elseif strcmp(line(1:10),'DEFINITION')
                    pos_RefSeq  = strfind(line,'(RefSeq)');
                    if isempty(pos_RefSeq)
                        prot_name = lower(line(13:end));
                        disp([gene_name ': no RefSeq'])
                    else
                        prot_name = lower(line(pos_RefSeq+9:end));
                    end

                %4th column: EC number
                elseif strcmp(line(1:9),'ORTHOLOGY')
                    pos_EC    = strfind(line,'[EC:');
                    if ~isempty(pos_EC)
                        EC_names  = line(pos_EC+4:end-1);
                    end

                %5th column and 7th column: MW & sequence
                elseif strcmp(line(1:5),'AASEQ')
                    end_seq  = false;
                    for k = j+1:length(text)
                        if length(text{k}) > 10
                            if strcmp(text{k}(1:5),'NTSEQ')
                                end_seq = true;
                            end
                        end
                        if ~end_seq
                            sequence = [sequence text{k}];
                        end
                    end
                    MW = calculateMW(sequence);

                %6th column: pathway
                elseif strcmp(line(1:7),'PATHWAY')
                    start    = strfind(line,organism_code);
                    pathway  = line(start(1):end);
                    end_path = false; 
                    met_path_string = [organism_code ...
                                       '01100  Metabolic pathways'];                      
                    while ~end_path
                        for k = j+1:length(text) 
                            nospace  = strrep(text{k},met_path_string,'');
                            nospace  = strrep(nospace,' ','');
                            if length(nospace) > 10
                                if strcmp(nospace(1:3),organism_code) && ...
                                   ~end_path
                                   start    = strfind(text{k},organism_code);
                                   pathway  = [pathway ' ' ...
                                               text{k}(start(1):end)];
                                else
                                    end_path = true;
                                end
                            end
                        end
                    end
                end
            end
        end
        %Create one aditional association in kegg:
        n         = n+1;
        kegg{n,1} = uni;
        kegg{n,2} = prot_name;
        kegg{n,3} = gene_name;
        kegg{n,4} = EC_names;
        kegg{n,5} = MW;
        kegg{n,6} = pathway;
        kegg{n,7} = sequence;
        cd ([DB_path '/KEGG/' organism_code])
        disp(['Updating KEGG database: Ready with gene ' gene_name])
    end
    kegg(n+1:end,:)         = [];

    %Save all databases as .mat files:
    cd ../..
    save(['hsa_ProtDatabase' '.mat'],'kegg','swissprot');
    cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
