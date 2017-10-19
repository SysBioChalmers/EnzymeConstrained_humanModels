%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEGGEcRetrieve
% Access the KEGG API and retrieves all data available for each enzyme 
% related to an EC number.
%
% Relates each EC number in KEGG to a metabolic subgroup according to the 
% next classification:
%
%   Carbon and Energy Metabolism: 
%                               - carbohydrate metabolism
%                               - Energy metabolism
%   Amino Acids, Fatty acids and Nucleotides metabolism:
%                               - Lipid metabolism
%                               - Nucleotide metabolism
%                               - Amino Acid metabolism
%   Intermediate and Secondary metabolism:
%                               - Metabolism of other Amino acids
%                               - Glycan biosynthesis and metabolism
%                               - Metabolism of cofactors and vitamins
%
%   INPUTS:     
%         - None
%   OUTPUTS:
%         - EC_cells:   Array with substrates/products, pathways, genes and 
%                       metabolic subgroups related to each EC number found
%                       in KEGG enzyme database.
%         - pathGroups: Array containing the name and the pathways of each
%                       metabolic subgroup.
%
%
% Ivan Domenzain. created:     2017-07-12
% Ivan Domenzain. Last edited: 2017-10-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EC_cells, pathGroups] = retrieveKEGGEcs

    % Reads the previously downloaded enzymes list from the KEGG FTP
    %path = 'write your repo path here'
    path = '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM';
    cd ([path '/Databases']);
    fileID_uni     = fopen('enzyme.txt');
    enzymes_data  = textscan(fileID_uni,'%s','delimiter','\n');
    fclose('all');
    enzymes_data = enzymes_data{1}; 
    %Data for all the EC numbers are lumped together on the file, being
    %delimited by '///'
    delimiters = find(~cellfun(@isempty,strfind(enzymes_data, '///')));
    pos = 0;
    EC_cells = [];

    %Main loop that goes through all the different EC numbers in the enzyme
    %file
    for i=1:length(delimiters)
        
        %Extract all data available in the enzyme file for the EC number
        ECdata         = enzymes_data(pos+1:delimiters(i)-1);
        %Extract EC numbers
        row            = (strsplit(ECdata{1},' '));
        EC_cells{1}{i} = row{3};
        %Extract substrates
        cells=[];
        substrate = find(~cellfun(@isempty,strfind(ECdata,'SUBSTRATE')));
        product   = find(~cellfun(@isempty,strfind(ECdata,'PRODUCT')));
        comment   = find(~cellfun(@isempty,strfind(ECdata,'COMMENT')));
        history   = find(~cellfun(@isempty,strfind(ECdata,'HISTORY')));
        
        %cells     = ECdata(substrate:min([comment,history])-1);
        cells     = ECdata(substrate:product-1);
        if ~isempty(substrate)
            %Omits the word 'SUBSTRATE' of the substrates list
            cells{1} = cells{1}(13:strfind(cells{1},'[')-2);
            if length(cells)>1
                for j=2:length(cells)
                   %if strfind(cells{j},'PRODUCT')
                   %    cells{j} = cells{j}(13:strfind(cells{j},'[')-2);
                   %else
                        cells{j} = cells{j}(1:strfind(cells{j},'[')-2);
                   %end
                end
            end
        end
        EC_cells{2}{i} = lower(cells);        

        %Extract Pathways
        pathway   = find(~cellfun(@isempty,strfind(ECdata,'PATHWAY')));
        orthology = find(~cellfun(@isempty,strfind(ECdata,'ORTHOLOGY')));  
        genes     = find(~cellfun(@isempty,strfind(ECdata,'GENES')));
        links     = find(~cellfun(@isempty,strfind(ECdata,'DBLINKS')));
        cells=[];
        if ~isempty(pathway) 
            %Saves all the pathways rows between 'PATHWAY' and the next
            %field and omits the word 'PATHWAY' from the first row
            cells     = ECdata(pathway:min([orthology, genes, links])-1);
            cells{1} = cells{1}(13:end);
        end
        EC_cells{3}{i} = cells; 
        
        %Extract genes, omitting the word 'GENES'
        cells=[];
        if ~isempty(genes)      
            cells     = ECdata(genes:links-1);
            cells{1} = cells{1}(13:end);
        end
        EC_cells{4}{i} = cells; 
        
        %Updates the position in the enzyme file
        pos = delimiters(i);         
    end
    cd ../ComplementaryScripts/KcatDistributions
    [EC_cells, pathGroups] = KEGGPathsRetrieve(EC_cells);
end

    
   
    
    


