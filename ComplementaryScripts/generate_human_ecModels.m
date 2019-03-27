current = pwd;
load('../models/humanGEM_cellLines/11models.mat')
modelNames = who;
modelNames = modelNames(~strcmpi(modelNames,'current'));
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