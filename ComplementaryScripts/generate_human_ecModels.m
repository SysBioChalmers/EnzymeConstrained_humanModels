current = pwd;
load('../models/humanGEM_cellLines/11models.mat')
modelNames = who;
modelNames = modelNames(~strcmpi(modelNames,'current'));
for i=1:length(modelNames)
    cellName = modelNames{i};
    [ecModel,model_data,kcats] = enhanceGEM_cellLine(cellName);
end
cd (current)
rmdir('GECKO', 's')