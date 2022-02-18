function provePAcostante (nameFile)
%PROVEPACOSTANTE
% Il file provePAcostante1.txt contiene tre colonne. La prima corrisponde
% al volume (D = 3 dimensioni), la seconda a (deltaTau)^4 e la terza a
% PACC. Si vuole verificare che PACC rimanga circa costante.

if nargin < 1
    nameFile = 'provePAcostante1.txt';
end

data = load(nameFile,'\t');
figure()
scatter3(data(:,1),data(:,2),data(:,3));
title('Probabilità di accettanza costante');
xlabel('Volume');
ylabel('(\delta\tau)^4');
zlabel('P_{ACC}');