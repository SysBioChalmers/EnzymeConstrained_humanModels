%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kcat_distributions
% 
% Shows the cumulative distribution for enzymatic parameters from the BRENDA 
% database for the organisms specified by the user. The overall distributions 
% will be splitted into distributions for different metabolic subgroups. 
%
% Ivan Domenzain.  Last edited: 2017-11-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %path = 'write your repo path here'
    path = '/Users/ivand/Documents/EnzymeConstrained-HMR-GEM';
    
    %Extract the Data available in KEGG enzyme database
     %Structure with all the organisms in KEGG and its corresponging
    %phylogenetic domain
    cd ([path '/ComplementaryScripts/KcatDistributions'])
    load('EC_path_taxonomy.mat')

    %Extracts all the enzymatic data (Kcat,Km, SA, Mw) queried from BRENDA
    %each value in these files corresponds to the maximum value found for
    %the especific EC number / substrate / organism triplet.
    BRENDAfiles_path=[path '/Databases/BRENDA_data'];
    %Choose one parameter to analyze {'KCAT', 'KM', 'SA','MW'}
    parameter_name   = 'KCAT';
    %Extract info
    [parameter_data] = enzymes_data(BRENDAfiles_path,parameter_name);
    %Organisms to be searched in the BRENDA database
    organisms = string({'homo sapiens'});
   %Initialize cell arrays size: (number of organisms)x(number of metabolic
   %subgroups+1) 
    for i=1:length(organisms)
        for j=1:length(pathGroups{1})+1
            E_parameter{i}{j} = [];
        end
    end
   %Parses the file for the selected parameter and gets the distribution
    %of values
    E_parameter = getDistributions(parameter_data, EC_cells, organisms, ...
                                   pathGroups, E_parameter);
    % Extract the Kcat values of the EC model sorted by metabolic blocks
    cd ([path '/ComplementaryScripts/KcatDistributions'])
     kcat_values{1} = matchedKcatsDist(model,model_kcats, pathGroups);
    
    met_group   = string({'All enzymes',...
                          'Central carbon and energy metabolism',...
                          'Amino acids, Fatty acids and Nucleotides metabolism',...
                          'Intermediate and Secondary metabolism'});
   % Plots the cumulative distribution of the selected parameter in BRENDA
   titleStr    = [organisms 'K_{cat} in BRENDA']; 
   [y_param, stats_param] = plotSplittedCumDist(E_parameter,met_group,titleStr);
   % Plots the cumulative distribution of the Kcats matched to the ECmodel
   titleStr               = [organisms 'ecModel K_{cat}?s'];                                        
   [y_param, stats_param] = plotSplittedCumDist(kcat_values,met_group,titleStr);
    % Shows a comparison of the indicated distributions
    Distcell{1} = {E_parameter{1}{1} kcat_values{1}{1}};
    names    = string({'BRENDA all enzymes', 'hsa_ecModel all enzymes'});
    titleStr               = ['K_{cat} distributions comparison'];                                        
    [y_param, stats_param] = plotSplittedCumDist(Distcell,names,titleStr);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_param, stats_param] = plotSplittedCumDist(E_parameter,legends,...
                                                      titlestr)
   figure
   for i=1:length(E_parameter)
    for j=1:length(E_parameter{i})
        if ~isempty(E_parameter{i}{j})
        	str(j) = legends(j)+' ('+string(length(E_parameter{i}{j}))+...
                     ' / ' + string(median(E_parameter{i}{j}))+ ')';
            [y_param(j), stats_param(j)] = cdfplot(E_parameter{i}{j});
            title(titlestr)
            ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
            xlabel('K_{cat} [s^{-1}]','FontSize',30,'FontWeight','bold');
            set(gca, 'XScale', 'log')
            hold on
        end
    end
   end
   legend(y_param,str);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_parameter= getDistributions(parameter_data, EC_cells, organisms,...
                                       pathGroups, E_parameter)

    for i=1:length(parameter_data{1})
    	%Filter EC number match just if it is also present in KEGG (enzyme)
    	ec_indx = find(strcmpi(EC_cells{1},parameter_data{1}{i}(3:end))~=0);
        if ~isempty(ec_indx)
        	%Looks into the KEGG data if the substrate of the entry appears 
            %on the list of natural substrates/products for the EC number
            subs_indx = find(strcmpi(EC_cells{2}{ec_indx},...
                                     parameter_data{2}(i))~=0);
            if ~isempty(subs_indx)
            	name   = extract_orgName(parameter_data{3}{i});

                %For each organisms especified by the user
                for j=1:length(organisms)
                    if strcmpi(organisms{j},name)
                    	E_parameter{j}{1}    = [E_parameter{j}{1};...
                                                parameter_data{4}(i)];
                        %Looks for the pathway associated to each 
                        %element on the ec_indx subset on the metabolic
                        %subgroups and if found appends its K value to 
                        %the corresponding subgroup distribution
                        if ~isempty(EC_cells{5}{ec_indx})
                        	for k = 1:length(EC_cells{5}{ec_indx})
                                subDist = find(strcmpi(EC_cells{5}{ec_indx}(k),...
                                                       pathGroups{1}));
                                E_parameter{j}{subDist+1} =...
                                  [E_parameter{j}{subDist+1}; parameter_data{4}(i)];  
                            end
                        end
                    end
                end
            end
        end

   end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function orgName = extract_orgName(str)

    pos1 = strfind(str,'/');
    pos2 = strfind(str,' ');
    if length(pos2)>1
    	pos  = min(pos1(1),pos2(2));
    else
    	pos = pos1(1);
    end
    orgName = str(1:pos-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that queries a metabolic pathway from KEGG and looks if the EC #
% is related to it.

function presence = ECnum_presence(pathway,ec)
   
    URL       = 'http://rest.kegg.jp/get/'+string(pathway);
    path_data = webread(URL);
    path_data = textscan(path_data,'%s','delimiter','\n');
    %Extracts just the lines that relate genes and proteins information
    gene_row = indexes_string(path_data{1},'GENE ',false);
    % Looks if either compound or reference appears first on the data in
    % order to extract just the GENE related part.
    ending   = [indexes_string(path_data{1},'COMPOUND ',true),...
                indexes_string(path_data{1},'REFERENCE ',true)];
    ending   = min(ending);
    extract  = path_data{1}(gene_row:ending-1);
    presence = indexes_string(extract,ec,true);
    if isempty(presence)
        presence = false;
    else
        presence = true;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag
        matching = find(~cellfun(@isempty,matching),1);
    else
        matching = find(~cellfun(@isempty,matching));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that queries the full KEGG organism list and returns a cell array
% containing KEGG codes and names 
function orgList = KEGG_orgList
    URL     = 'http://rest.kegg.jp/list/organism';
    data    = webread(URL);
    data    = textscan(data,'%s %s %s %s','delimiter','\t');
    orgList{1} = data{2};
    orgList{2} = cellfun(@lower,data{3},'UniformOutput',false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract all the enzymatic data (Kcat,Km, SA, Mw) queried from BRENDA
%each value in these files corresponds to the maximum value found for the
%especific EC number / substrate / organism triplet.
function parameter_data = enzymes_data(BRENDAfiles_path,parameter_name)
     cd (BRENDAfiles_path)
     if strcmpi(parameter_name,'KM')
     	file_name = ['min_'  parameter_name  '.txt'];
     else
        file_name = ['max_'  parameter_name  '.txt'];
     end
     disp(['Parsing' parameter_name 'values from BRENDA'])
     fID       = fopen(file_name);
     parameter_data = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
     fclose('all');    

     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
