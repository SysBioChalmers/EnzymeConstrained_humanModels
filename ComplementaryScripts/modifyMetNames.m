%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modifyMetNames(model)
% Receives the HMR model .mat structure and modify the metNames field in 
% order to make it compatible with the substrate names of the max_Kcats
% file created with the kinetic data from BRENDA.
%
% Ivan Domenzain. Last edited: 2017-10-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modifyMetNames(model)

    model.OriginalmetNames = model.metNames;
    metList                = model.metNames;
    
    for i=1:length(metList)
        met = model.metNames{i};
        % Removes the '_termination code' of each metabolite
        pos = strfind(met,'_');
        if ~isempty(pos)
            met = met(1:pos(end)-1);
        end
        
        if any(strfind(met,'[protein]-'))
            met = replace(met,'[protein]-','');
            
        elseif(any(strfind(met,'-[protein]')))
            met = replace(met,'-[protein]','');
            
        elseif(any(strfind(met,'[protein C terminal]-')))
            met = replace(met,'[protein C terminal]-','');
        
        end
        model.metNames{i} = met;
    end
    
end
