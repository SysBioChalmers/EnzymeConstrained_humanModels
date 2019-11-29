%generate_human_ecModels_NCI60
% 
% Main script for enhancement of 11 cell-line specific GEMs (from the
% NCI60 cell-lines) with enzyme-constraints with the use of the GECKO
% pipeline.
%
% Ivan Domenzain.      Last edited: 2019-11-20

%Clone GECKO and substitute human-specific scripts
git('clone https://github.com/SysBioChalmers/GECKO.git')
%Replace scripts in GECKO:
fileNames = dir('GECKO_humanFiles');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
        fullName   = ['GECKO_humanFiles/' fileName];
        GECKOpath = dir(['GECKO/**/' fileName]);
        GECKOpath = GECKOpath.folder;
        copyfile(fullName,GECKOpath)
    end
end
clear
clc
load('../models/humanGEM_cellLines/11models.mat')
modelNames = who;
current    = pwd;
GECKO_path =  [current '/GECKO'];
%Generate enzyme-constrained models
for i=1:length(modelNames)
    cd (current)
    cellName = modelNames{i};
    mkdir (['../models/' cellName])
    cd (['../models/' cellName])
    save([cellName '.mat'],cellName)
    mkdir Data
    cd (current)
    %Models mat files are saved in their respective folder by enhanceGEM_cellLine
    [ecModel,ecModel_batch] = enhanceGEM_cellLine(cellName);
end
cd (current)
rmdir('GECKO', 's')