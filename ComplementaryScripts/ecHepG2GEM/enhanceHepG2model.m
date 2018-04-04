function ecHepG2Model = enhanceHepG2model
% enhanceHepG2model
%   Converts the HepG2 specific GEM to an Enzyme constrained version
%
%
%   Ivan Domenzain, 2018-04-04
%
    current = pwd;
    cd ../models/HepG2
    %Import HepG2 model xml file and save it as a matlab structure
    HepG2model = importModel('Hep-G2.xml');
    save('HepG2model.mat','HepG2model')
end