%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ecModel = modifyKcats(ecModel,ecModel_batch,PDB,limKcats, gRexp)
%
% Function that gets the limiting Kcat values in an EC model (according to
% a sensitivity analysis), then it modifies each of those values according to 
% the maximal available values in the BRENDA files (Kcats and SA*Mw) when a
% manual curated option is not specified.
%
% Ivan Domenzain    Last edited. 2017-01-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ecModel = modifyKcats(ecModel,ecModelBatch,PDB,gRexp,batch_flag)
    current = pwd;
    T  = []; error   = []; changed_kcats = [];
    e  = -100; i=1; 
    %Load BRENDA data:
    KCAT_file          = 'max_KCAT.txt';
    SA_file            = 'max_SA.txt';
    MW_file            = 'max_MW.txt';
    [BRENDA,SA_cell]   = loadBRENDAdata(KCAT_file,SA_file,MW_file);

    %Iterates while growth rate is being underpredicted
    disp('********************Limiting Kcats curation********************')
    while e<=0
        cd (current)
        %Get the top growth rate-limiting enzyme (uniprot code basis)
        [limKcat,limRxns,breakFlag] = MaxLimitingKcat(model_min,changed_kcats);
        
        if breakFlag == false
            [model_min,data] = changeKcat(model_min,limKcat,gRexp,...
                                                       PDB,BRENDA,SA_cell);
                                          
            cd (current)
            T              = [T; data];
            e              = data{1,7};
            %Add a string with the uniprot code and the rxn number in order
            %to keep track of the modified coefficients
            str            = {horzcat(data{1},'_',num2str(limKcat{3}))};
            changed_kcats  = [changed_kcats; str];
            disp(str)
            disp(limKcat{6}(1))
            error          = [error; e];           
            str            = ['#' num2str(i) ' prev:' ...
                              num2str(data{1,5}) ' new:' num2str(data{1,6}) ...
                              ' CC:' num2str(limKcat{1,5}) ' Err:' num2str(e)];
            disp(str)          
            i = i+1;            
        else
            break
        end
    end
    
    cd (current)
    [m,n]     = size(ecModel.S);
    ecModel.S = ecModelBatch.S(1:m,1:n);
    T         = cell2table(T,'VariableNames',{'Unicode','Organism',...
                                              'Modified','Parameter',...
                                              'oldValue','newValue',...
                                              'error','gRControlCoeff'});
    writetable(T, 'KcatModifications.txt');  
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,data] = changeKcat(model,limKcats,gR_exp,PDB,BRENDA,SA_cell)
   % Gets the Unicode 
    UniCode  = limKcats{1}{1}(strfind(limKcats{1}{1},'_')+1:end);
    % Map the UNIPROT code (kcat)
    [ECnumber, Mw] = findECnumber(UniCode,PDB);
    enzIndx = limKcats{2}(1);
    rxnIndx = limKcats{3}(1);
    data    = {UniCode,[],[],[],[],[],[],[]};
    e       = 0;
    
    if ~isempty(ECnumber) && ~isempty(Mw)
    	flag           = false;
    	previous_value = -1/(3600*model.S(enzIndx,rxnIndx)); %[1/s]
        
            disp(['Automatic search // ' 'EC#: ' ECnumber])
           %Looks for the maximal value available for the respective EC
           %number (Kcats and SA*Mw if available)
            [Kcat,org, match] = findMaxValue(ECnumber,BRENDA,SA_cell,Mw);
            coeff             = -1/(Kcat);  
           %Change the kinetic coefficient just if a higher value was found
            if coeff > model.S(enzIndx,rxnIndx)
            	flag = true;
                model.S(enzIndx,rxnIndx) = coeff;
            end
        end
        new_value = -1/(3600*model.S(enzIndx,rxnIndx));
           
        % After changing the i-th kcat limiting value a simulation is
        % performed and the growth rate and absolute error are saved 
        model_sim            = model;
        gR_pos               = find(strcmpi(model_sim.rxnNames,'growth'));
        model_sim.c          = zeros(size(model_sim.c));
        model_sim.c(gR_pos)  = 1;
        solution             = solveLP(model_sim);
        model_sim.lb(gR_pos) = 0.999*solution.x(gR_pos);
        model_sim.ub(gR_pos) = solution.x(gR_pos);
        solution             = solveLP(model_sim,1);
        e                    = ((solution.x(gR_pos)-gR_exp)/gR_exp)*100;
        
        data                 = {UniCode,org,flag,match,previous_value,...
                                new_value,e,limKcats{5}(1)};
    end 
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [value, organism, parameter] = findMaxValue(EC_cell,BRENDA,...
                                                                SA_cell,Mw)
    %Looks for the maximum turnover number available for the EC# associated
    %with the uniprot code
    EC_cell    = strsplit(EC_cell,' ');
    value      = [];
    organism   = [];
    parameter  = [];
    for i=1:length(EC_cell)
        find_flag  = false;
        %In case that wild cards are present in the EC number the search on
        %the BRENDA file will be relaxed.
        if ~isempty(strfind(EC_cell{i},'-'))
             EC_cell{i} = EC_cell{i}(strfind(EC_cell{i},'-')-1:end);
             find_flag  = true;
        end    
        ECnumber = ['EC' EC_cell{i}];
        Kcat     = 0; orgK = '';
        if find_flag == true
            matching = indexes_string(BRENDA{1},ECnumber,false);
        else
            % If no wild cards are present the EC number search in the
            % BRENDA file will look for an exact match
            matching = find(strcmpi(ECnumber,BRENDA{1}));
        end
        %Gets the maximum Kcat value for the queried EC#
        if ~isempty(matching)
            [Kcat, maxIndx] = max(BRENDA{4}(matching));
            orgK            = BRENDA{3}(matching(maxIndx));
        end        
        % Looks for the maximum SA*Mw value available for the EC number
        SA_Mw = 0; orgS = '';
        if find_flag == true
            matching = indexes_string(SA_cell{1},ECnumber,false);
        else
            matching = find(strcmpi(ECnumber,SA_cell{1}));
        end
        %Gets the maximum SA*Mw value for the queried EC#
        if ~isempty(matching)
            [SA_Mw, maxIndx] = max(SA_cell{4}(matching));
            SA_Mw            = SA_Mw; %[1/hr]
            orgS             = SA_cell{3}(matching(maxIndx));
        end        
        %Choose the maximal available value as a turnover number for the EC
        value  = [value; max(Kcat,SA_Mw)]; 

        if Kcat>SA_Mw
            organism  = [organism; {orgK}];
            parameter = [parameter; 'K_cat'];
        else
        	organism  = [organism; {orgS}];
            parameter = [parameter; 'SA*Mw'];
        end   
    end
    [value, maxIndx] = max(value);
    organism         = organism(maxIndx);
    parameter        = parameter(maxIndx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ECnumber, Mw] = findECnumber(Unicode,ProtDatabase)
    cd ../../../Databases
    load (ProtDatabase)
    DB1{1} = swissprot(:,1);DB1{2} = swissprot(:,4);DB1{3} = swissprot(:,5);
    DB2{1} = kegg(:,1);     DB2{2} = kegg(:,4);     DB2{3} = kegg(:,5);
    ECnumber = {};
    % First look for the UNIPROT ID in the swissprot DB structure
    matching = indexes_string(DB1{1},Unicode,true);  
    if ~isempty(matching)
        ECnumber = DB1{2}{matching};
        Mw       = DB1{3}{matching};
    end
    % If nothing comes up then look into the KEGG DB structure
    if isempty(ECnumber)
        matching = indexes_string(DB2{1},Unicode,true);
        if ~isempty(matching)
            ECnumber = DB2{2}{matching};
            Mw       = DB1{3}{matching};
        end
    end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [KCATcell, SAcell] = loadBRENDAdata(KCAT_file,SA_file,MW_file )
     current = pwd;
     cd ../../../Databases
     %Extract BRENDA DATA from files information
     scallingFactor = 3600;   %[1/s] -> [1/h]
     KCATcell      = openDataFile(KCAT_file,scallingFactor); 
     scallingFactor = 60;     %[umol/min/mg] -> [mmol/h/g]
     SA             = openDataFile(SA_file,scallingFactor); 
     scallingFactor = 1/1000; %[g/mol] -> [g/mmol]
     MW             = openDataFile(MW_file,scallingFactor); 
     
     for i=1:4
         SAcell{i} = [];
     end
     previousEC = []; EC_indexes = [];
     for i=1:length(SA{1})
         %Gets the indexes of the EC repetitions in the MW cell for every
         %new (different) EC
         if ~strcmpi(SA{1}(i), previousEC)
             EC_indexes = find(strcmpi(SA{1}(i),MW{1}));
         end
         mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
         % just looks for the first match because just the maximal value for
         % each EC# / Orgaism is reported on the file
         org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
         if ~isempty(org_index)
             SAcell{1} = [SAcell{1};SA{1}(i)];
             SAcell{2} = [SAcell{2};SA{3}(i)];
             SAcell{3} = [SAcell{3}; SA{4}(i)* mwEC{2}(org_index)]; %[1/hr]
             SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];
         end
         previousEC = SA{1}(i);
     end
     cd (current)
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function data_cell = openDataFile(fileName,scallingFactor)
     fID          = fopen(fileName);
     data_cell    = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
     fclose(fID);
     data_cell{4} = data_cell{4}*scallingFactor;
     %Split string for each organism in the BRENDA data 
     %{name, taxonomy, KEGG code}
     data_cell{3}  = cellfun(@stringSplit, data_cell{3});
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
function string_cells = stringSplit(cell_array)
        string_cells = {strsplit(cell_array,'//')};
        string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
