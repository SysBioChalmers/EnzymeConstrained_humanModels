%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modelProteinContent(model)
%
% Receive an ecModel and an absolute proteomics dataset and calculate the
% total protein concentration (gProt/gDw) that was included in the model. 
% 
% Ivan Domenzain. Last edited: 2018-01-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totalProtModel, f] = modelProteinContent(model)
	%Get total protein content according to the biomass reaction
	Ptotal    = getPtotal(model); %[g/gDw] 
	[f,count] = measureAbundance(model.enzymes);% [g Pmodel/g Ptotal(DB)] 
	%Load Absolute proteomics data
	cd ../Databases/HepG2proteomics
	file_name = 'HepG2Proteomics_qualified.txt';
	fID       = fopen(file_name);
    data      = textscan(fID,'%s %f %f','delimiter',',');
    fclose('all'); 

    data{2}        = data{2}*(1/1000); %[pmol/mgDw] -> [mmol/gDw]
    meanConc       = mean(data{2});
    enzymes   	   = model.enzymes;
	MWs		       = model.MWs;
	totalProtModel = 0;
	for i=1:length(enzymes)
		
		index = find(strcmpi(enzymes(i),data{1}),1);
		if ~isempty(index)
			enzConcentration = data{2}(index)*MWs(i);
			totalProtModel	 = totalProtModel+enzConcentration;
		end
	end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function totalMass = getPtotal(model) 

%'human_proteinPool', ('0.0937 alanine[c] + 0.0507 arginine[c] + ...
%                       0.04 asparagine[c] + 0.04 aspartate[c] + ...
%                       0.0142 cysteine[c] + 0.1118 glutamine[c] + ...
%                       0.1831 glycine[c] + 0.0198 histidine[c] + ...
%                       0.0309 isoleucine[c] + 0.0664 leucine[c] + ...
%                       0.0571 lysine[c] + 0.0156 methionine[c] + ...
%                       0.029 phenylalanine[c] + 0.0853 proline[c] + ...
%                       0.0491 serine[c] + 0.0402 threonine[c] + ...
%                       0.019 tyrosine[c] + 0.0471 valine[c] + ...
%                       0.0072 tryptophan[c]) => proteinPool[c]'
%'HumanGrowth', Other precurssors + 6 proteinPool[c]  => biomass[c]' [g/gDw h]

% AMINOACIDS 				  MWs [g/mol]
% A	Alanine 					 71.08
% B	Aspartic acid or Asparagine 114.60
% C	Cysteine  				    103.14
% D	Aspartic acid 				115.09
% E	Glutamic acid 				129.11
% F	Phenylalanine 				147.17
% G	Glycine  					 57.05  
% H	Histidine 					137.14
% I	Isoleucine 					113.16
% J	Leucine or Isoleucine 		113.16
% K	Lysine 						128.17
% L	Leucine 					113.16
% M	Methionine 					131.20
% N	Asparagine 					114.10
% O	Pyrrolysine 				255.31
% P	Proline 				     97.12
% Q	Glutamine 					128.13
% R	Arginine 					156.19
% S	Serine 						 87.08
% T	Threonine 					101.10 
% U	Selenocysteine  			150.04 
% V	Valine 						 99.13 
% W	Tryptophan 					186.21
% X	any 						126.50
% Y	Tyrosine 					163.17
% Z	Glutamic acid or Glutamine  128.62

totalMass = 6*(18+ 0.0937*71.08 + 0.0507*156.19 + 0.04*114.10+ 0.04*114.60+ ...
            0.0142*103.14 + 0.1118*128.13 + 0.1831*57.05 + 0.0198*137.14 + ...
            0.0309*113.16 + 0.0664*113.16 + 0.0571*128.17 + 0.0156*131.20+ ...
			0.029*147.17 + 0.0853*97.12 + 0.0491*87.08 + 0.0402*101.10 + ...
			0.019*163.17 + 0.0471*99.13 + 0.0072*186.21)/1000;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

