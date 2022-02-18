%% ANALYZE BINS

% Import data
data1 = load("provaHMC.txt",'\t');
% Facciamo l'esempio con la magnetizzazione
V = 64;
%counter = data(:,1);
m1 = data1(:,2)/V; % magnetizzazione
naccu = 100;
%abs_m = data(:,3);
%m2 = data(:,4);
%m4 = data(:,5);

% Set iniziale grafico
shg
clf
set(gcf,'numbertitle','off','name','Plot errore vs bin');

axis square
title("naccu = 100, osservabile: magnetizzazione");

hold on;

for t = 1:100
    % Make averages every n values and add them to the plot
    n = t*10; % average every n values
    % a = reshape(mag,[],1); inutile qui, poichè abbiamo già un vettore: a=mag
    meanVec = arrayfun(@(i) mean(m1(i:i+n-1)),1:n:length(m1)-n+1)'; % the averaged vector
    
    % Variante 1
    errVec = sqrt(var(meanVec)/size(meanVec,1));
    scatter(n*naccu,mean(errVec),'.');
    
    % Variante 2: stesso risultato della variante 1
    %res = (meanVec-mean(mag)).^2;
    %res = sqrt(mean(res)/(size(meanVec,1)-1));
    %scatter(n*naccu,res,'.');
end
hold off;
clear;

% h.XData = data(:,1);
% h.YData = data(:,2);
% drawnow

%% AUTO-CORRELATION FUNCTIONS

% Import data
data1 = load("provaHMC.txt",'\t');
%counter = data1(:,1);
V = 64;
m1 = data1(:,2)/V;
%absMag = data1(:,3);
%M2 = data1(:,4);
%M4 = data1(:,5);

meanMag1 = mean(m1);

% Set iniziale grafico
shg
clf
set(gcf,'numbertitle','off','name','Funzione di auto-correlazione');

axis square
title("naccu = 100, osservabile: magnetizzazione");

hold on;

N = length(m1);

for t = 1:100
    gamma = 0;
    for i = 1:N-t
        gamma = gamma + (m1(i)-meanMag1)*(m1(i+t)-meanMag1);
    end
    scatter(t,gamma/(N-t),'.');
end

hold off;
clear;

% h.XData = x(:,1);
% h.YData = x(:,2);
% drawnow

%% AUTO-CORRELATION FUNCTIONS FOR DIFFERENT DELTAs

% Set iniziale grafico
shg
clf
set(gcf,'numbertitle','off','name','Funzioni di auto-correlazione');
axis square
ttl = title("naccu = 100, osservabile: magnetizzazione");
hold on;


% Import data
V = 64;
data1 = load("provaBinningD025.txt",'\t');
m1 = data1(:,2)/V;
meanMag1 = mean(m1);
data2 = load("provaBinningD022.txt",'\t');
mag2 = data2(:,2)/V;
meanMag2 = mean(mag2);
data3 = load("provaBinningD040.txt",'\t');
mag3 = data3(:,2)/V;
meanMag3 = mean(mag3);
data4 = load("provaBinningD010B1.txt",'\t');
mag4 = data4(:,2)/V;
meanMag4 = mean(mag4);
% data4 e data5 sono gli stessi ma in data5 binniamo di 10 in 10 e poi
% plottiamo, osservando che la funzione di auto-correlazione va a 0 più
% velocemente
data5 = load("provaBinningD010B10.txt",'\t');
mag5 = data5(:,2)/V;
n = 10;
mag5 = arrayfun(@(i) mean(mag5(i:i+n-1)),1:n:length(mag5)-n+1)';
meanMag5 = mean(mag5);


N = 10000;
points = ones(100,1);
points2 = ones(100,1);
points3 = ones(100,1);
points4 = ones(100,1);
points5 = ones(100,1);

for t = 1:100
    gamma = 0;
    gamma2 = 0;
    gamma3 = 0;
    gamma4 = 0;
    gamma5 = 0;
    for i = 1:N-t
        gamma = gamma + (m1(i)-meanMag1)*(m1(i+t)-meanMag1);
        gamma2 = gamma2 + (mag2(i)-meanMag2)*(mag2(i+t)-meanMag2);
        gamma3 = gamma3 + (mag3(i)-meanMag3)*(mag3(i+t)-meanMag3);
        gamma4 = gamma4 + (mag4(i)-meanMag4)*(mag4(i+t)-meanMag4);
        gamma5 = gamma5 + (mag5(i)-meanMag5)*(mag5(i+t)-meanMag5);
    end
    points(t) = gamma/(N-t);
    points2(t) = gamma2/(N-t);
    points3(t) = gamma3/(N-t);
    points4(t) = gamma4/(N-t);
    points5(t) = gamma5/(N-t);
end
t = 1:100;
scatter(t,points,'.','r');
scatter(t,points2,'.','g');
scatter(t,points3,'.','b');
scatter(t,points4,'.','c');
scatter(t,points5,'.','k');
legend("D = 0.25","D = 0.22","D = 0.40","D = 0.10, B = 1", "D = 0.10, B = 10");

hold off;
clear;
% h.XData = x(:,1);
% h.YData = x(:,2);
% drawnow

%% TRYING TO COMPUTE INTEGRATED AUTO-CORRELATION TIME (TAU)
%tic
% Import data
data1 = load("provaBinningD025.txt",'\t');
%V = 64;
m1 = data1(:,2);%/V;
meanMag1 = mean(m1);

N = 10000;
endpoint = 110; %massimo N-1

gamma0 = 0;
for i = 1:N
    gamma0 = gamma0 + (m1(i)-meanMag1)*(m1(i)-meanMag1);
end
gamma0 = gamma0/N;

tau = 0;
points = ones(1,endpoint+1);
points(1) = gamma0;
for t = 1:endpoint
    gamma = 0;
    for i = 1:N-t
        gamma = gamma + (m1(i)-meanMag1)*(m1(i+t)-meanMag1);
    end
    points(t+1) = gamma/(N-t);
    gamma = gamma/(N-t);
    tau = tau + gamma/gamma0;
end
tau = 0.5*(1+2*tau);
fprintf("endpoint t = %d, tau = %f\n",endpoint,tau);

% Set grafico
shg
clf
set(gcf,'numbertitle','off','name','Funzione di auto-correlazione');
hold on;
axis square
title(sprintf("naccu = 100 || endpoint: t = %d || osservabile: magnetizzazione",endpoint));
plot(1:endpoint+1,points,'.');
hold off;
%clear;
%toc
[~,~,~,tauint,dtauint,~,wopt] = ...
    UWerrTexp(data1,[],[],'TestM',2);
fprintf("tauint = %f, dtauint = %f, wopt = %f\n",tauint,dtauint,wopt);

%% EXAMPLES USING UWerrTexp

data = load("provaBinningD025.txt",'\t');
data(:,1) = [];
[value,dvalue,ddvalue,tauint,dtauint,qval,wopt,gammaFbb,drho] = ...
    UWerrTexp(data,[],[],"Test magnetizzazione",1);
hold on;
title("Confronto UWerrTexp funzione di auto-correlazione");
plot(gammaFbb,'r');
N = 10000;
points = ones(100,1);
m1 = data(:,1);
meanMag1 = mean(m1);
for t = 1:100
    gamma = 0;
    for i = 1:N-t
        gamma = gamma + (m1(i)-meanMag1)*(m1(i+t)-meanMag1);
    end
    points(t) = gamma/(N-t);
end
t = 1:100;
scatter(t,points,'.','b');
legend("UWerrTexp","Analitica");
hold off


%% CHECK JACKKNIFE

% Import data
data1 = load("evoluzioneOsservabili1.txt",'\t');
data1(:,1) = []; % la prima colonna è il #passi, non mi serve qui
data1(:,1) = []; % la seconda colonna contiene i valori medi di m, che qui non mi servono esplicitamente
N = size(data1,1);
nObs = size(data1,2);
V = 64;

% Calcolo medie (valori di aspettazione) di |M|, M^2, M^4
means = mean(data1(:,1:nObs));

% jackknife vectors for primary variables: |M|, M^2, M^4
primary_jacks = means(1:3)-(data1(1:N,1:3)-means(1:3))/(N-1);
primary_errors = sqrt(var(primary_jacks)*((N-1)^2)/N);

% Theoretical jackknife for derived variables: susceptibility, binder cumulant
susc_jack = susceptibility(primary_jacks(:,1),primary_jacks(:,2),V);
susc_error = sqrt(var(susc_jack)*((N-1)^2)/N);
bc_jack = binderCumulant(primary_jacks(:,2),primary_jacks(:,3));
bc_error = sqrt(var(bc_jack)*((N-1)^2)/N);

% Explicit calculation of errors
cov_susc = mean((data1(1:N,1)-means(1)).*(data1(1:N,2)-means(2)));
cov_bc = mean((data1(1:N,2)-means(2)).*(data1(1:N,3)-means(3)));
primary_stdv = primary_errors*sqrt(N);
susc_error2 = sqrt((primary_stdv(2)^2 + 4*(means(1)*primary_stdv(1))^2 - 4*means(1)*cov_susc)/(V^2)/N);
bc_error2 = sqrt((primary_stdv(3)^2 + 4*(means(3)*primary_stdv(2)/means(2))^2 - 4*means(3)*cov_bc/means(2))/(means(2)^4)/N);

% Più rapidamente
c_susc = cov(data1(:,1),data1(:,2));
susc_error3 = sqrt((c_susc(2,2) + 4*(means(1)^2)*c_susc(1,1) - 4*means(1)*c_susc(1,2))/(V^2)/N);
c_bc = cov(data1(:,2),data1(:,3));
bc_error3 = sqrt((c_bc(2,2) + 4*c_bc(1,1)*(means(3)/means(2))^2 - 4*means(3)*c_bc(1,2)/means(2))/(means(2)^4)/N);

% Ora, check con UWerrTexp
% Cambia le "Run" con 'Nome Plot' per generare i plot
[UWabs,UWabsErr,UWabsErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(1)/V);
[UWm2,UWm2Err,UWm2ErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(2)/V);
[UWsusc,UWsuscErr,UWsuscErrErr] = UWerrTexp(data1,[],[],"Run",@(x) (x(2)-x(1)^2)/V);
[UWbinder,UWbinderErr,UWbinderErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(3)/x(2)^2);
UWvalues = [UWabs,UWm2,UWsusc,UWbinder];
UWerrors = [UWabsErr,UWm2Err,UWsuscErr,UWbinderErr];
UWerrorsErrors = [UWabsErrErr,UWm2ErrErr,UWsuscErrErr,UWbinderErrErr];

% Print
format long g;
fprintf("Controllo errori variabili primarie:\n");
display([primary_errors;sqrt(var(data1(:,1:3))/N)]);
fprintf("Differenze:\n");
display(primary_errors-sqrt(var(data1(:,1:3))/N));
fprintf("UW con errori:\n");
display([UWerrors(1:2);UWerrorsErrors(1:2)]);
fprintf("Controllo errori variabili secondarie:\n");
display([susc_error,bc_error,;susc_error2,bc_error2]);
fprintf("Differenze:\n");
display([susc_error-susc_error2,bc_error-bc_error2]);
fprintf("UW con errori:\n");
display([UWerrors(3:4);UWerrorsErrors(3:4)]);


%% ANALISI DATI

% Import data
data1 = load("provaHMC.txt",'\t');
data1(:,1) = []; % la prima colonna è il #passi, non mi serve qui
data1(:,1) = []; % la seconda colonna contiene i valori medi di m, che qui non mi servono esplicitamente
N = size(data1,1);
nObs = size(data1,2);
V = 64;

% Calcolo medie (valori di aspettazione) di |M|, M^2, M^4
means = mean(data1(:,1:nObs));

% jackknife vectors for primary variables: |M|, M^2, M^4
primary_jacks = means(1:3)-(data1(1:N,1:3)-means(1:3))/(N-1);
primary_errors = sqrt(var(primary_jacks)*((N-1)^2)/N);

% jackknife for derived variables: susceptibility, binder cumulant
susc_jack = susceptibility(primary_jacks(:,1),primary_jacks(:,2),V);
susc_error = sqrt(var(susc_jack)*((N-1)^2)/N);
bc_jack = binderCumulant(primary_jacks(:,2),primary_jacks(:,3));
bc_error = sqrt(var(bc_jack)*((N-1)^2)/N);

% Results for observables

% Summary
observables = [means(1:2)/V,susceptibility(means(1),means(2),V),binderCumulant(means(2),means(3))];
errors = [primary_errors(1:2)/V,susc_error,bc_error];
format long g;
display([observables;errors]);

% Ora, check con UWerrTexp
% Cambia le "Run" con 'Nome Plot' per generare i plot
[UWabs,UWabsErr,UWabsErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(1)/V);
[UWm2,UWm2Err,UWm2ErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(2)/V);
[UWsusc,UWsuscErr,UWsuscErrErr] = UWerrTexp(data1,[],[],"Run",@(x) (x(2)-x(1)^2)/V);
[UWbinder,UWbinderErr,UWbinderErrErr] = UWerrTexp(data1,[],[],"Run",@(x) x(3)/x(2)^2);
UWvalues = [UWabs,UWm2,UWsusc,UWbinder];
UWerrors = [UWabsErr,UWm2Err,UWsuscErr,UWbinderErr];
UWerrorsErrors = [UWabsErrErr,UWm2ErrErr,UWsuscErrErr,UWbinderErrErr];
fprintf("CONTROLLO UWERRTEXP\n");
display([UWvalues;UWerrors;UWerrorsErrors]);

clear;

% Functions definitions at the end of scripts
function s = susceptibility (absM,M2,V)
s = (M2-absM.^2)/V;
end
function b = binderCumulant (M2,M4)
b = M4./(M2.^2);
end

