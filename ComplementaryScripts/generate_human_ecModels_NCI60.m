load('../models/humanGEM_cellLines/11models.mat')
modelNames = who;
current    = pwd;
%Clone GECKO and substitute human-specific scripts
git('clone https://github.com/SysBioChalmers/GECKO.git')
GECKO_path =  [current '/GECKO'];
%Replace scripts in GECKO:
fileNames = dir('GECKO_humanFiles');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..')
        fullName   = ['GECKO_humanFiles/' fileName];
        GECKOpath = dir(['GECKO/**/' fileName]);
        GECKOpath = GECKOpath.folder;
        copyfile(fullName,GECKOpath)
    end
end
for i=1:length(modelNames)
    cd (current)
    cellName = modelNames{i};
    mkdir (['../models/' cellName])
    cd (['../models/' cellName])
    save([cellName '.mat'])
    mkdir Data
    cd (current)
    [ecModel,model_data,kcats] = enhanceGEM_cellLine(cellName);
end
cd (current)
rmdir('GECKO', 's')