%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEGGPathsRetrieve
% Extracts KEGG pathways, sorts them and relate the different EC numbers
% (also in KEGG) to a metabolic pathways category.

% Ivan Domenzain. Last edited: 2017-10-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EC_cells, pathGroups] = KEGGPathsRetrieve(EC_cells)
    %path = 'write your repo path here'
    path = '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM';
    %Read and save the data in pathways list file (downloaded  from the 
    %KEGG API.
    cd ([path '/Databases']);
    fileID_uni     = fopen('pathway.txt');   
    path_data  = textscan(fileID_uni,'%s','delimiter','\n');
    fclose('all');
    path_data = path_data{1}; 
    cat = 0; pathways=[];
    for i=2:length(path_data)
        marker =find(strfind(path_data{i},'##'));
         if ~isempty(marker)
             cat = cat+1;
             pathways{1}{cat} =  path_data{i};
             pathways{2}{cat} = [];
         elseif isempty(find(strfind(path_data{i},'#'), 1))
             letters = find(isletter(path_data{i}));
             pathways{2}{cat} =  [pathways{2}{cat};...
                                  string(path_data{i}(letters(1):end))];
         end
    end
    %Sort selected pathways into three main metabolic groups
    pathGroups{1} = {'Carbohydrate and Energy primary metabolism',...
       'Aminoacid, Fatty acid and Nucleotide primary metabolism',...
        'Intermediate and Secondary metabolism'};
    pathGroups{2}{1} = [pathways{2}{2};pathways{2}{3}];
    pathGroups{2}{2} = [pathways{2}{4};pathways{2}{5};pathways{2}{6}];
    pathGroups{2}{3} = [pathways{2}{7};pathways{2}{8};pathways{2}{9}];
    
    %Creates a new dimension in EC_cells which contains the metabolic 
    %subgroup(s) of each EC number.
    for i=1:length(EC_cells{3})
        EC_cells{5}{i} = [];
        for j=1:length(EC_cells{3}{i})
            %EC_cells{5}{i}{j} = [];
            str = EC_cells{3}{i}{j}(min(strfind(EC_cells{3}{i}{j},' '))+2:end);
            for k=1:length(pathGroups{2})
                if any(strcmpi(pathGroups{2}{k},str))
                    if ~ismember(pathGroups{1}{k},EC_cells{5}{i})
                        EC_cells{5}{i} = [EC_cells{5}{i};pathGroups{1}(k)];
                    end
                end
            end
        end
    end
end
        
    
    
    