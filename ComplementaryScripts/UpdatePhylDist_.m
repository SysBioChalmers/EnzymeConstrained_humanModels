%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phylDistStruct] = UpdatePhylDist_(keggPath,unicellular)
%   
%   Calculates distance between species in KEGG based on systematic name.
%
%   REQUIREMENTS
%     A valid username and password for the KEGG FTP
%     ftp://ftp.bioinformatics.jp/
%
%   INPUT
%     keggPath        if keggPhylDist.mat is not in the RAVEN\external\kegg
%                     directory, this function will attempt to read data<<
%                     from a local FTP dump of the KEGG database. taxPath
%                     is the path to the root of this database.
%
%     unicellular     if true it generates a matrix with distance = Inf to
%                     all the multicellular organisms (Animals, Plants, 
%                     Protists)
%
%   OUTPUT   
%     phylDistStruct  a structure with a list of organism ids, names  and a 
%                     matrix that specifies their pairwise distances
%
%   This simple metric is based on the number of nodes that two organisms 
%   are away from each other in the KEGG taxonomical tree.
%
%   Usage: phylDistStruct = getPhylDist(keggPath,unicellular)
%
%   Created:           Eduard Kerkhoven, 2017-02-28
%   Last modification: Ivan Domenzain,   2017-06-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phylDistStruct] = UpdatePhylDist_(keggPath,unicellular)

    % Sets the path for the Taxonomy file that will be parsed
    taxPath = [keggPath '/genes/misc'];   
    if nargin<2
        unicellular=true;
    end
    % First it looks for the PhylDist.mat structure in the current folder
    % if it does not exist then the structure will be created
    if exist('PhylDist.mat')
        fprintf(['NOTE: Importing KEGG phylogenetic distance matrix from '...
                 'keggPhylDist.mat' '.\n']);
        load('PhylDist.mat');
    else
       fprintf(['Cannot locate ' 'PhylDist.mat ' ...
           ' and will try to generate it from the local KEGG database.\n']);
       
       % It accesses the KEGG ftp directory and looks for the taxonomy file     
       if ~exist(fullfile(taxPath,'taxonomy'),'file')   
           EM = fprintf(['The file ''taxonomy'' cannot be located at '...
                        strrep(taxPath,'\','/')...
                        '/ and should be downloaded from the KEGG FTP.\n']);
           dispEM(EM);
       else
           phylDistStruct.ids   = {};
           phylDistStruct.names = {};
           %Keeps the categories for each organism
           orgCat               = {};
           currentCat           = {}; %Keeps track of the current category
           depth                = 0; %Keeps track of the current tree depth
           %Loop through the file
           orgCounter           = 0;
           fid                  = fopen(fullfile(taxPath,'taxonomy'), 'r');
           %fileOrgs             = fopen('KEGGorganisms.txt','w');
                         
           while 1
             %Get the next line
             tline = fgetl(fid);
             %Abort at end of file
             if ~ischar(tline)
                 break;
             end

             if any(tline)
                %Check if it's a new category
                if tline(1)=='#'
                    %Find the first space (=depth +1)
                    sPos=strfind(tline,' ')-1; %Should always exist
                    sPos=sPos(1);
                    %If we have stepped back one step in the tree
                    if sPos<depth
                        currentCat=currentCat(1:sPos);
                    end
                    depth=sPos;
                    currentCat{depth}=tline(sPos+2:end);
                else
                    orgCounter=orgCounter+1;
                    %It is an organism
                    %Gets the organism name and KEGG ID
                    sPos=find(isstrprop(tline, 'wspace')); %Should always exist
                    phylDistStruct.ids{orgCounter}   = tline(sPos(1)+1:sPos(2)-1);
                    phylDistStruct.names{orgCounter} = tline(sPos(3)+1:end);
                    orgCat{orgCounter}               = currentCat;
                end
             end
          end

          %Generate a distance matrix 
          phylDistStruct.distMat = zeros(numel(phylDistStruct.ids));
          for i=1:numel(phylDistStruct.ids)
              % Creates a .txt file with KEGG IDs and names
              %nbytes = fprintf(fileOrgs,'%i\t%s\t%s\n',i,phylDistStruct.ids{i},...
              %                 phylDistStruct.names{i});
                         
              for j=1:numel(phylDistStruct.ids)
                  %If organism is unicellular then the distance to 
                  %multicellular organisms is set to Inf
                  if unicellular == true && (strcmpi('animals',orgCat{j}(2))||...
                    strcmpi('plants',orgCat{j}(2)) ||strcmpi('protists',orgCat{j}(2))) 
                    phylDistStruct.distMat(i,j)=Inf;
                  else
                  	% Distance calculation
                    k = 1;
                    %Counts the common categories between organisms.
                    while strcmpi(orgCat{i}(k),orgCat{j}(k)) 
                    	k=k+1;
                        if k>min(length(orgCat{i}),length(orgCat{j}))
                            break;
                        end
                    end
                    % The distance is the total depth of both organisms
                    % minus twice the common depth in the tree.
                    phylDistStruct.distMat(i,j) =  numel(orgCat{i})+...
                                                   numel(orgCat{j})- 2*(k-1);
                  end
              end
          end
       end    
   end
   %fclose(fileOrgs);
   %dlmwrite('phylDistMatrix.txt',phylDistStruct.distMat,'delimiter','\t')
   %Saves the structure
   currentFolder = pwd;
   save('PhylDist.mat','phylDistStruct');

end
