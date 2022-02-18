function autocorrPlot (nameFile,data1)

if nargin < 1 
    nameFile = 'provaHMC.txt';
    data1 = load(nameFile,'\t');
end
if nargin < 2 data1 = load(nameFile,'\t'); end
% Import data
%counter = data1(:,1);
m1 = data1(:,1);
%absMag = data1(:,3);
%M2 = data1(:,4);
%M4 = data1(:,5);

meanMag1 = mean(m1);

% Set iniziale grafico
figure()
axis square
title("autocorrPlot");
hold on;

N = length(m1);

nMax = N/10;
points = ones(nMax,1);

for t = 1:nMax
    gamma = 0;
    for i = 1:N-t
        gamma = gamma + (m1(i)-meanMag1)*(m1(i+t)-meanMag1);
    end
    points(t) = gamma/(N-t);
end

t = 1:nMax;
scatter(t,points,'.','r');
hold off;

end