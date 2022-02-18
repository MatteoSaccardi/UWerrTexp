function PAerfc (nameFile,T,PACCnameFile,wantchecks)
%PAerfc con PAerfc.txt generato da PAerfc.c
%
% PAerfc.txt (nameFile) è un file che contiene 4 colonne di dati,
% corrispondenti a un contatore e ai valori di <deltaH>, <|deltaH|>,
% <exp(-deltaH)> ogni naccu step a deltaTau fissato. Ogni T (100) valori,
% viene cambiato deltaTau.
% Qui, calcoliamo le medie e i loro errori sugli intervalli da 100 e
% facciamo il plot della probabilità di accettanza rispetto a tali valori.
% Vogliamo verificare l'andameno di PA = erfct[0.5sqrt(<dH>)].
%
% Esempi di chiamate:
% PAerfc
% PAerfc('PAerfc.txt',100,'PACCDATAPAerfc.txt')
% PAerfc('PAerfc.txt',100,'PACCDATAPAerfc.txt',1)
%

if (nargin < 4)
    wantchecks = 0;
end
if nargin < 3
    PACCnameFile = 'PACCDATAPAerfc.txt';
end
if (nargin < 2)
    T = 100; %periodo
end
if (nargin < 1)
    nameFile = 'PAerfc.txt';
end

data = load(nameFile,'\t');
data(:,1) = []; % elimino la colonna dei conteggi
data(:,2) = []; % elimino la colonna di <|dH|>
data(:,2) = []; % elimino la colonna di <exp(-dH)>
N = length(data(:,1));

% calcoliamo le medie sui blocchi da T (100) con i rispettivi errori
dHvec = arrayfun(@(i) mean(data(i:i+T-1)),1:T:length(data(:)-T+1));
%errorsVec = arrayfun(@(i) sqrt(var(data(i:i+T-1))),1:T:length(data(:)-T+1));
errorsVec = arrayfun(@(i) sqrt(var(data(i:i+T-1))/T),1:T:length(data(:)-T+1));

if wantchecks
    % I primi dati hanno probabilità di accettanza bassa, possono essere
    % correlati; ne terremo conto
    checkFile = fopen("PAerfc_checkFile.txt","w");
    for i = 1:N/T
        [~,~,~,tauint,dtauint] = UWerrTexp(data(T*(i-1)+1:T*i),[],[],0);
        % Ci attendiamo che sia intorno a 0.5 se non c'è autocorrelazione
        fprintf(checkFile,"tauint da %f a %f:\n%f +- %f\n",T*(i-1)+1,T*i,tauint,dtauint);
    end
    
    autocorrPlot([],data(1:100,1));

    % Controlliamo che gli altri errori non differiscano da quelli
    % calcolati di troppo
    dHvecCheck = zeros(N/T,1);
    errorsVecCheck = zeros(N/T,1);
    errorsErrors = zeros(N/T,1);
    for j = 1:N/T
        [dHvecCheck(j),errorsVecCheck(j),errorsErrors(j)] = ...
            UWerrTexp(data(T*(j-1)+1:T*j),[],[],1,1);
    end
    fprintf(checkFile,"Differenza medie: %f\n",norm(dHvec'-dHvecCheck,'fro')); % dovrebbe venire 0
    nSigma = 2;
    A = [abs(errorsVec'-errorsVecCheck)-nSigma*errorsErrors abs(errorsVec'-errorsVecCheck)+nSigma*errorsErrors];
    fprintf(checkFile,"\n");
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            fprintf(checkFile,"%f\t",A(i,j));
        end
        fprintf(checkFile,"\n");
    end
    % dovrebbero essere a sx negativi e a dx positivi affinché 0 sia
    % compatibile entro 2 sigma con errorsVec-errorsVecCheck i.e. affinché
    % questi due valori siano compatibili i.e. no autocorrelazione da
    % tenere in considerazione.
end

PACC = load(PACCnameFile,'\t');
expectedPACC = erfc(0.5*sqrt(dHvec));
% Se
%           P = erfc[0.5*sqrt(dH)]
% allora
%           P'(dH) = dP/d(dH) = 1/(4*sqrt(dH))*[-2/sqrt(pi)*e^(-dH/4)]
% e l'errore su P dovuto all'errore su dH (jackknife) è
%           sP(dH) = |P'(dH)] * errH
PACCprime = (-1*(2/sqrt(pi)/4)./sqrt(dHvec)).*exp(-dHvec/4);
errPACC = abs(PACCprime).*errorsVec;
figure()
hold on
xlabel("\langle\DeltaH\rangle");
ylabel("P_{ACC}");
title("P_{ACC} = erfc[1/2 sqrt{\langle\DeltaH\rangle}]");
errorbar(dHvec,expectedPACC,errorsVec,'horizontal','.r','Markersize',10);
errorbar(dHvec,PACC,errPACC,errPACC,errorsVec,errorsVec,'.b','Markersize',10);
legend('Expected values','Data');

% Fit con svdfit
% Ci attendiamo un comportamento lineare per PACC vicino a 1 i.e. per <dH>
% basso
% Dall'interpolazione escludiamo i due con deltaH maggiore, siccome qui
% l'approssimazione non è più buona.
% Il fit non viene comunque buono, però si nota che diminuendo i punti
% (i.e. con un'approssimazione lineare teoricamente sempre migliore) il chi
% quadro dimiuisce, fino a diventare circa 1, e l'intercetta diventa sempre
% più compatibile con 1, come atteso.
% E' inoltre evidente la deviazione dal comportamento lineare a valori di
% <dH> grandi i.e. PACC molto diversi da 1-(piccolo).
figure()
hold on
xlabel("(\delta\tau)^2");
ylabel("P_{ACC}");
title("Linear Fit with 95% Prediction Interval");
inizio = 2;
i = 0:9;
deltaTau = 1./(5+2.*i); % vedi hmc.c
x = deltaTau(inizio:end).^2;
xq = linspace(0,deltaTau(inizio)^2);
y = PACC(inizio:end);
erry = errPACC(inizio:end);
[a,err_a,covmat,chisqr] = svdfit(x,y,erry,1);
str = sprintf("Chi^2 = %f\nIntercetta: %e +- %e",chisqr,a(1),err_a(1));
dim = [.2 .2 .1 .1];
annotation("textbox",dim,"String",str,"FitBoxToText","on");
% Plottiamo (estrapoliamo)
errorbar(deltaTau.^2,PACC,errPACC,'.r','Markersize',10);
[yq,err_yq] = svdpolyval(a,xq,covmat);
plot(xq,yq,'b-');
plot(xq,yq+2*err_yq,'m--',xq,yq-2*err_yq,'m--')

legend('Data','Linear Fit','95% Prediction Interval');

hold off


end
