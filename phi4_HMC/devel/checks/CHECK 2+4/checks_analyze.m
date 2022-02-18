function checks_analyze (nameFile,T,wantchecks)
%ANALISI CHECK 2 E CHECK 4 CON checks_analyze.txt
%
% checks_analyze.txt (nameFile) è un file che contiene 4 colonne di dati,
% corrispondenti a un contatore e ai valori di <deltaH>, <|deltaH|>,
% <exp(-deltaH)> ogni naccu step a deltaTau fissato. Ogni T (100) valori,
% viene cambiato deltaTau.
% Qui, calcoliamo le medie e i loro errori sugli intervalli da 100 e
% facciamo il plot di tau VS questi valori. Per <|deltaH|> ci attendiamo un
% andamento quadratico, per <exp(-deltaH)> ci attendiamo che tale valore
% sia compatibile con 1.
% I risultati possono essere inseriti in un file .txt per essere
% eventualmente usati altrove e.g. su gnuplot.
%
% Esempi di chiamate:
% checks_analyze
% checks_analyze('checks_analyze_default.txt',100,1)
% checks_analyze('checks_analyzeT10000.txt',10000,0)
% checks_analyze('checks_analyzeT10000.txt',10000,1)
%

if (nargin < 3)
    wantchecks = 0;
end
if (nargin < 2)
    T = 100; %periodo
end
if (nargin < 1)
    nameFile = 'checks_analyze_default.txt';
end

resultsFile = fopen("results_check2e4.txt","w");

data = load(nameFile,'\t');
data(:,1) = []; % elimino la colonna dei conteggi
N = length(data(:,1));

meansVec = zeros(N/T,3);
errorsVec = zeros(N/T,3);
for j = 1:3
    % calcoliamo le medie sui blocchi da T (100) con i rispettivi errori
    meansVec(:,j) = arrayfun(@(i) mean(data(i:i+T-1,j)),1:T:length(data(:,j)-T+1));
    errorsVec(:,j) = arrayfun(@(i) sqrt(var(data(i:i+T-1,j))/T),1:T:length(data(:,j)-T+1));
end

% NON calcoliamo il valore medio sulle medie col rispettivo errore poiché
% sono simulazioni con deltaTau diverso
% NO meanVec = mean(meansVec);
% NO errorVec = sqrt(cov(meansVec)/length(meanVec));

i = 0:9;
deltaTau = 1./(5+2.*i)'; % confronta con hmc.c

if wantchecks
    % I primi dati hanno probabilità di accettanza bassa, possono essere
    % correlati; ne terremo conto
    checkFile = fopen("checks_analyze_checkFile.txt","w");
    fprintf(checkFile,"CHECK TAUINT\n");
    for i = 1:N/T
        [~,~,~,tauint1,dtauint1] = UWerrTexp(data(T*(i-1)+1:T*i,1),[],[],0);
        % Ci attendiamo che sia intorno a 0.5 se non c'è autocorrelazione
        [~,~,~,tauint2,dtauint2] = UWerrTexp(data(T*(i-1)+1:T*i,2),[],[],0);
        [~,~,~,tauint3,dtauint3] = UWerrTexp(data(T*(i-1)+1:T*i,3),[],[],0);
        fprintf(checkFile,"tauint da %f a %f:\n%f +- %f\n%f +- %f\n%f +- %f\n"...
            ,T*(i-1)+1,T*i,tauint1,dtauint1,tauint2,dtauint2,tauint3,dtauint3);
    end
    fprintf(checkFile,"FINE CHECK TAUINT\n\n");
    
    autocorrPlot([],data(1:100,1));
    
    figure()
    errorbar(deltaTau,meansVec(:,1),errorsVec(:,1),'.r','Markersize',14);
    xlabel('\delta\tau');
    ylabel('\langle\DeltaH\rangle');
    title('\langle\DeltaH\rangle');
    
    figure()
    errorbar(deltaTau,meansVec(:,2),errorsVec(:,2),'.r','Markersize',14);
    xlabel('\delta\tau');
    ylabel('\langle|\DeltaH|\rangle');
    title('\langle|\DeltaH|\rangle');
    
    figure()
    errorbar(deltaTau,meansVec(:,3),errorsVec(:,3),'.r','Markersize',14);
    xlabel('\delta\tau');
    ylabel('\langleexp(-\DeltaH)\rangle');
    title('\langleexp(-\DeltaH)\rangle');
    hold on
    line([deltaTau(1),deltaTau(end)],[1,1]);
    hold off
    
    % Stampiamo i risultati
    fprintf(checkFile,"STAMPA RISULTATI\n");
    fprintf(checkFile,"deltaTau, <dH> +- err, <|dH|> +- err, <exp(-dH)> +- err\n");
    for i = 1:length(deltaTau)
        fprintf(checkFile,'%f \t %f \t %f \t %f \t %f \t %f \t %f \n',...
        deltaTau(i),meansVec(i,1),errorsVec(i,1),meansVec(i,2),...
        errorsVec(i,2),meansVec(i,3),errorsVec(i,3));
    end
    fprintf(checkFile,"FINE STAMPA RISULTATI\n\n");

    % Controlliamo che gli altri errori non differiscano da quelli 
    % calcolati di troppo
    fprintf(checkFile,"INIZIO CHECK RISULTATI CON UWerrTexp\n");
    meansVecCheck = zeros(N/T,3);
    errorsVecCheck = zeros(N/T,3);
    errorsErrors = zeros(N/T,3);
    for i = 1:size(data,2)
        for j = 1:N/T
            [meansVecCheck(j,i),errorsVecCheck(j,i),errorsErrors(j,i)] = ...
                UWerrTexp(data(T*(j-1)+1:T*j,i),[],[],1,1);
        end
    end
    fprintf(checkFile,"Differenza medie: %f\n",norm(meansVec-meansVecCheck,'fro')); % dovrebbe venire 0
    nSigma = 2;
    A = [abs(errorsVec-errorsVecCheck)-nSigma*errorsErrors abs(errorsVec-errorsVecCheck)+nSigma*errorsErrors];
    aux = A(:,4); A(:,4) = A(:,2); A(:,2) = aux;
    aux = A(:,3); A(:,3) = A(:,5); A(:,5) = aux;
    aux = A(:,3); A(:,3) = A(:,4); A(:,4) = aux;
    fprintf(checkFile,"\n");
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            fprintf(checkFile,"%f\t",A(i,j));
        end
        fprintf(checkFile,"\n");
    end
    fprintf(checkFile,"FINE CHECK RISULTATI CON UWerrTexp\n\n");
    % 10 righe, 3 coppie di colonne adiacenti; ci si attende che
    % il primo valore sia negativo e il secondo positivo, dato che questi
    % sono gli estremi dell'intervallo di confidenza entro cui i valori
    % degli errori sono compatibili tra loro (i.e. la loro differenza è
    % compatibile con zero).
    
    % esempio con checks_analyze_default:
%{
A =
    -0.000558	0.000692	-0.000534	0.000729	-0.000700	0.000736
    -0.000234	0.000258	-0.000039	0.000335	-0.000244	0.000435
    -0.000043	0.000184	-0.000063	0.000211	0.000011	0.000214
    -0.000082	0.000086	-0.000048	0.000106	-0.000073	0.000091
    -0.000061	0.000124	-0.000075	0.000078	-0.000064	0.000117
    -0.000024	0.000059	-0.000041	0.000046	-0.000027	0.000057
    -0.000036	0.000038	-0.000043	0.000056	-0.000033	0.000039
    -0.000023	0.000026	-0.000033	0.000034	-0.000025	0.000043
    -0.000020	0.000022	-0.000013	0.000026	-0.000020	0.000022
    -0.000021	0.000023	-0.000013	0.000024	-0.000017	0.000018
%}
    % qui, l'unico valore che devia è nella terza riga. Andando a vedere il
    % tau0int, troviamo che vale circa 0.35 in tutto il terzo blocco di
    % misure. Questo è l'unico blocco che crea problemi. A 3 sigma, non ci
    % sono più problemi, mentre a 1 sigma tutto il terzo blocco ne dà.
    
    % esempio con checks_analyze (T = 10000, più statistica)
%{
    -0.000002	0.000011	-0.000005	0.000007	-0.000007	0.000012	
    -0.000003	0.000004	-0.000002	0.000003	-0.000003	0.000004	
    -0.000000	0.000003	-0.000001	0.000003	-0.000001	0.000002	
    -0.000001	0.000001	-0.000001	0.000001	-0.000001	0.000001	
    -0.000000	0.000001	-0.000000	0.000001	-0.000000	0.000001	
    -0.000000	0.000001	-0.000000	0.000000	-0.000000	0.000001	
    -0.000000	0.000000	-0.000000	0.000000	-0.000000	0.000000	
    -0.000000	0.000000	-0.000000	0.000000	-0.000000	0.000000	
    -0.000000	0.000000	-0.000000	0.000001	-0.000000	0.000000	
    -0.000000	0.000000	-0.000000	0.000000	-0.000000	0.000000	
%}
    % Come si vede, aumentando la statistica non si hanno più problemi, nè
    % nel tauint nè nella compatibilità degli errori!
end

% Dall'analisi della funzione di autocorrelazione per le varie osservabili
% e notando che la probabilità di accettanza per il primo valore di
% deltaTau è 0.68 (molto minore rispetto a tutte le altre: la seconda è di
% 0.85, la terza 0.92, per arrivare a valori di circa 0.96-0.98), si nota
% che l'autocorrelazione non è esattamente nulla per tutte le osservabili.
% Tuttavia, gli errori che commettiamo nel calcolare gli errori sono
% trascurabili, ANCHE PER IL TERZO PUNTO (entro nSigma).
% Potremmo comunque, per questo valore, stimare l'errore in maniera esatta.
% Per fare ciò, potremmo utilizzare la funzione UWerrexp.

%{
for i = 1:size(data,2)
    for j = 3 %j = 1:N/T per cambiarli tutti
        [meansVec(j,i),errorsVec(j,i)] = ...
                UWerrTexp(data(T*(j-1)+1:T*j,i),[],[],1,1);
        errorsVec(j,i) = errorsVec(j,i)*sqrt(T);
        %errore come radice della varianza
    end
    %[meansVec(1,i),errorsVec(1,i)] = UWerrTexp(data(1:T,i),[],[],1,1);
    %[meansVec(2,i),errorsVec(2,i)] = UWerrTexp(data(T+1:2*T,i),[],[],1,1);
end
%}

% Possiamo anche controllare che l'algoritmo svdfit funzioni confrontando
% i risultati con quelli ottenuti da svdfitLinear, implementazione
% esplicita del fit lineare con gli stessi output
if wantchecks
    errMax = 1e5*eps; %1e-11
    inizio = 1;
    j = 2;
    x = deltaTau(inizio:end).^2;
    y = meansVec(inizio:end,j);
    erry = errorsVec(inizio:end,j);
    [a,err_a,covmat,chisqr] = svdfit(x,y,erry,1);
    [a_check,err_a_check,covmat_check,chisqr_check] = svdfitLinear(x,y,erry);
    if norm(abs(a-a_check),'fro') > errMax || norm(abs(err_a-err_a_check),'fro') > errMax || ...
            norm(abs(covmat-covmat_check),'fro') > errMax || abs(chisqr-chisqr_check) > errMax
        warning('check svdfitLinear non passato');
        fprintf(checkFile,"check svdfitLinear non passato con errore di soglia %e\n",errMax);
    else
        fprintf(checkFile,"check svdfitLinear passato con errore di soglia %e\n",errMax);
    end
end

% Fit DeltaH con svdfit
figure()
hold on
xlabel("(\delta\tau)^4");
ylabel("\langle\DeltaH\rangle");
title("Linear Fit of \langle\DeltaH\rangle with 95% Prediction Interval");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<dH>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
inizio = 3; %CAMBIARE QUI PER INIZIARE IL FIT DA UN PUNTO DIVERSO
x = deltaTau(inizio:end).^4;
xq = linspace(0,deltaTau(inizio)^4);
y = meansVec(inizio:end,1);
erry = errorsVec(inizio:end,1);
[a,err_a,covmat,chisqr] = svdfit(x,y,erry,1);% 2); per includere ordini superiori nel fit
str = sprintf("Chi^2 without %d point(s) = %f\nIntercetta: %e +- %e",...
            inizio-1,chisqr,a(1),err_a(1));
fprintf(resultsFile,sprintf("Without %d point(s)\n",inizio-1));
fprintf(resultsFile,"Intercetta: %f +- %f\nPendenza: %f +- %f\n",a(1),err_a(1),a(2),err_a(2));
dim = [.2 .5 .3 0.4];
annotation("textbox",dim,"String",str,"FitBoxToText","on");
% Plottiamo (estrapoliamo)
[yq,err_yq] = svdpolyval(a,xq,covmat);
errorbar(x,y,erry,'.','Markersize',10);
plot(xq,yq,'r-');
plot(xq,yq+2*err_yq,'m--',xq,yq-2*err_yq,'m--')
legend('Data','Linear Fit','95% Prediction Interval','Location','Southeast');
hold off


% Fit |DeltaH| con svdfit
figure()
hold on
xlabel("(\delta\tau)^2");
ylabel("\langle|\DeltaH|\rangle");
title("Linear Fit of \langle|\DeltaH|\rangle with 95% Prediction Interval");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<|dH|>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
inizio = 2; %CAMBIARE QUI PER INIZIARE IL FIT DA UN PUNTO DIVERSO
x = deltaTau(inizio:end).^2;
xq = linspace(0,deltaTau(inizio)^2);
y = meansVec(inizio:end,2);
erry = errorsVec(inizio:end,2);
%[a,err_a,covmat,chisqr] = svdfit(x,y,erry,1);
% OSS: se includiamo anche i termini di ordine (deltaTau)^4, interpolazione
% decisamente migliore (anche con inizio = 1)!!!
[a,err_a,covmat,chisqr] = svdfit(x,y,erry,2);% 2); per includere ordini superiori nel fit
str = sprintf("Chi^2 without %d point(s) = %f\nIntercetta: %e +- %e",...
            inizio-1,chisqr,a(1),err_a(1));
fprintf(resultsFile,sprintf("Without %d point(s)\n",inizio-1));
fprintf(resultsFile,"Intercetta: %f +- %f\nPendenza: %f +- %f\n",a(1),err_a(1),a(2),err_a(2));
dim = [.2 .5 .3 0.4];
annotation("textbox",dim,"String",str,"FitBoxToText","on");
% Plottiamo (estrapoliamo)
[yq,err_yq] = svdpolyval(a,xq,covmat);
errorbar(x,y,erry,'.','Markersize',10);
plot(xq,yq,'r-');
plot(xq,yq+2*err_yq,'m--',xq,yq-2*err_yq,'m--')
legend('Data','Linear Fit','95% Prediction Interval','Location','Southeast');
hold off

% Fit exp(-DeltaH) con svdfit
figure()
hold on
xlabel("\delta\tau");
ylabel("\langleexp(-\DeltaH)\rangle");
title("Linear Fit of \langleexp(-\DeltaH)\rangle with 95% Prediction Interval");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<exp(-dH)>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
inizio = 2; %CAMBIARE QUI PER INIZIARE IL FIT DA UN PUNTO DIVERSO
x = deltaTau(inizio:end);
xq = linspace(0,deltaTau(inizio));
y = meansVec(inizio:end,3);
erry = errorsVec(inizio:end,3);
[a,err_a,covmat,chisqr] = svdfit(x,y,erry,1);% 2); per includere ordini superiori nel fit
str = sprintf("Chi^2 without %d point(s) = %f\nIntercetta: %e +- %e\nPendenza: %e +- %e",...
            inizio-1,chisqr,a(1),err_a(1),a(2),err_a(2));
fprintf(resultsFile,sprintf("Without %d point(s)\n",inizio-1));
fprintf(resultsFile,"Intercetta: %f +- %f\nPendenza: %f +- %f\n",a(1),err_a(1),a(2),err_a(2));
dim = [.2 .1 .3 0.4];
annotation("textbox",dim,"String",str,"FitBoxToText","on");
% Plottiamo (estrapoliamo)
[yq,err_yq] = svdpolyval(a,xq,covmat);
errorbar(x,y,erry,'.','Markersize',10);
plot(xq,yq,'r-');
plot(xq,yq+2*err_yq,'m--',xq,yq-2*err_yq,'m--')
legend('Data','Linear Fit','95% Prediction Interval');
hold off


fclose(resultsFile);

end

%{
Fit |DeltaH|

% Adesso, facciamo una prova di interpolazione
% Vediamo ad occhio che il valore 1 è compatibile con il valore ottenuto
% da <exp(-deltaH)> per ogni valore di deltaTau.
% Interpoliamo linearmente <|deltaH|> rispetto a deltaTau^2.
% Dopodiché, estrapoliamo per deltaTau = 0.
% Possiamo escludere il primo punto, che ha un errore molto alto (anche
% tenendolo viene comunque; per tenerlo, mettere inizio = 1)

% SENZA ERRORI

inizio = 1;
x = deltaTau(inizio:end).^2;
y = meansVec(inizio:end,2);
erry = errorsVec(inizio:end,2);
[yp,S] = polyfit(x,y,1);
xq = linspace(0,deltaTau(inizio)^2);
[yq,delta] = polyval(yp,xq,S);
figure()
plot(x,y,'o');
hold on
plot(xq,yq,'r-');
plot(xq,yq+2*delta,'m--',xq,yq-2*delta,'m--')
xlabel('\delta\tau^2');
ylabel('\langle|\DeltaH|\rangle');
title('Linear Fit of Data with 95% Prediction Interval without weights');
legend('Data','Linear Fit','95% Prediction Interval');

%}