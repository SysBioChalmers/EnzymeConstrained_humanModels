%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sol = MAXmin_Optimizer(model,indx,coeff,blockIndx,FixedFlux)
%  
%
%
% Ivan Domenzain.      Last edited: 2018-03-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = MAXmin_Optimizer(model,indx,coeff,blockIndx,FixedFlux)
function sol = MAXmin_Optimizer(model,indx,coeff,blockIndx,FixedFlux)
             model.c       = zeros(length(model.c),1);
             model.c(indx) = coeff;
             if nargin > 3
                 model.ub(blockIndx) = FixedFlux;
             end
             sol = solveLP(model);
end