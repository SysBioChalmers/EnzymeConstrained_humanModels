function model = manualKcatCuration(model)
%Growth limiting enzyme Q16850 // EC#1.14.13.70
%The automatic method chose the highest value available, for Bos taurus,
%however additional information was found on BRENDA references.
%Ivan Domenzain 2018-09-13
Kcat = 9.8*60; % Aoyama, Y.,J Biochem. 1996 May;119(5):926-33. [1/min] -> [1/h] 
protIndex = find(contains(model.metNames,'Q16850'));
protRxns  = find(model.S(protIndex,:));
%Assign curated value
model.S(protIndex,protRxns) = -1/Kcat;
end