function inexact_vs_exact (V,nameFile,N,wantchecks)
%INEXACT_VS_EXACT
% inexact_vs_exact.txt è un file che contiene dati divisi a blocchi e in
% 5 colonne. Le colonne indicano a un contatore, <M>, <|M|>, <M^2>, <M^4>.
% Le righe si riferiscono a misure indipendenti ogni naccu. Ogni file
% è costituito da due blocchi consecutivi. Il primo presenta le misure
% effettuate con lo step di accept_reject, nel secondo lo si trascura.
% In ogni blocco, si effettuano N misure consecutive mantenendo deltaTau
% (eps) costante. Si effettuano misure per nMeasures valori di deltaTau
% diversi.
% Qui vogliamo calcolare la media delle osservabili nei due metodi per ogni
% valore di deltaTau, confrontarli e plottarli rispetto a questo.
% Ci aspettiamo che i valori ottenuti nelle due sezioni siano sempre più
% vicini per deltaTau --> 0 e che il valore medio non dipenda da deltaTau
% nel caso in cui si consideri lo step di accept_reject.
%
% Esempi di chiamate:
% inexact_vs_exact
% inexact_vs_exact(64,'inexact_vs_exact_default.txt',100,1)
% inexact_vs_exact(64,'inexact_vs_exactT10000.txt',10000,1)
% inexact_vs_exact(64,'inexact_vs_exactT10000.txt',10000) %Usa questo per un buon grafico
%

% Osserviamo che gli errori con T = 10000 si riducono molto, rendendo
% inoltre tutti i tauint in checkTauInt.txt 0.5 (entro l'errore)!!!

if nargin < 4
    wantchecks = 0;
end
if nargin < 3
    N = 100; % guarda la prima riga del file inexact_vs_exact.txt
end
if nargin < 2
    nameFile = 'inexact_vs_exact_default.txt';
end
if nargin < 1
    V = 64;
end

resultsFile = fopen("results_check3.txt","w");

data1 = load(nameFile,'\t');
data1(:,1) = []; %non ci interessa il conteggio
data1(:,1) = []; %non ci interessa <m>

% guarda nella funzione inexact_vs_exact in hmc.c per il range di deltaTau
% I valori di deltaTau sono stampati all'inizio di ogni misura nel file
% inexact_vs_exact.txt
deltaTau = 1./(5+3.*(0:9))';

nMeasures = length(data1)/N/2;
% Creiamo un vettore di dati da analizzare singolarmente. Il primo indice
% fa riferimento alla sezione, il secondo all'uso o meno del passo di
% accept_reject. Le componenti saranno 100x3
data = zeros(nMeasures,2,N,size(data1,2));
for i = 1:nMeasures
    data(i,1,:,:) = data1(2*(i-1)*N+1:2*(i-0.5)*N,:);
    data(i,2,:,:) = data1(2*(i-0.5)*N+1:2*i*N,:);
end

if wantchecks
    checkTauInt = fopen("checkTauInt.txt","w");
    for i2 = 1:size(data,2)
        for i1 = 1:size(data,1)
            A = reshape(data(i1,i2,:,:),size(data,3),size(data,4));
            for i4 = 1:size(data,4)
                [~,~,~,tauint,dtauint] = UWerrTexp(A(:,i4),[],[],0);
                fprintf(checkTauInt,"tauint at %d,%d for obs %d is = %f +- %f\n",i1,i2,i4,tauint,dtauint);
                % Ci attendiamo che sia intorno a 0.5 se non c'è autocorrelazione
            end
        end
    end
    fclose(checkTauInt);
end

means1 = zeros(nMeasures,4);
means2 = zeros(nMeasures,4);
errors1 = zeros(nMeasures,4);
errors2 = zeros(nMeasures,4);

for i = 1:nMeasures
    [means1(i,:),errors1(i,:)] = analyze(reshape(data(i,1,:,:),size(data,3),size(data,4)),V);
    [means2(i,:),errors2(i,:)] = analyze(reshape(data(i,2,:,:),size(data,3),size(data,4)),V);
end
% Errori come incertezze statistiche, non come errori sulle medie
%errors1 = error1*sqrt(size(data,3));

% confronto con metropolis
data_m = load('metropolis_inexact_vs_exact.txt','\t');
data_m(:,1:2) = [];
[means_m,errors_m] = analyze(data_m,V);

%esempio di plot: <m^2>
figure()
hold on
for i = 1:nMeasures
    errorbar(deltaTau.^2,means1(:,2),errors1(:,2),'ob');
    errorbar(deltaTau.^2,means2(:,2),errors2(:,2),'or');
end
yline(means_m(2),'k','metropolis','Linewidth',1);
yline(means_m(2)+errors_m(2),':k','Linewidth',1);
yline(means_m(2)-errors_m(2),':k','Linewidth',1);
title('\langlem^2\rangle');
xlabel('(\delta\tau)^2');
ylabel('\langlem^2\rangle');

% Fit con svdfit
%fit a&r
inizio = 1;
x = deltaTau(inizio:end).^2;
y1 = means1(inizio:end,2);
err_y1 = errors1(inizio:end,2);
[a,err_a,covmat,chisqr] = svdfit(x,y1,err_y1,1);
%text(0,0,"Intercetta: %f \pm %f",a(1),err_a(1));
xq = linspace(0,deltaTau(inizio)^2);
[yq,err_yq] = svdpolyval(a,xq,covmat);
plot(xq,yq,'b-');
plot(xq,yq+2*err_yq,'b--',xq,yq-2*err_yq,'b--')
%fit without a&r
y2 = means2(inizio:end,2);
err_y2 = errors2(inizio:end,2);
[a2,err_a2,covmat2,chisqr2] = svdfit(x,y2,err_y2,1);
str = sprintf("chi^2 a&r = %f\nchi^2 without a&r = %f",chisqr,chisqr2);
annotation("textbox",[.2 .1 .2 .2],"String",str,"FitBoxToText","on");
%text(0,0,"Intercetta: %f \pm %f",a2(1),err_a2(1));
[yq2,err_yq2] = svdpolyval(a2,xq,covmat2);
plot(xq,yq2,'r-');
plot(xq,yq2+2*err_yq2,'r--',xq,yq2-2*err_yq2,'r--')
legend('Data with a&r','Data without a&r','Location','east');
%{
legend('Data with a&r','Data without a&r','Linear Fit a&r',...
    '95% Prediction Interval a&r','Linear Fit without a&r',...
    '95% Prediction Interval without a&r','Location','east');
%}
str = sprintf("Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f",...
    means_m(2),errors_m(2),a(1),err_a(1),a2(1),err_a2(1));
annotation('textarrow',[0.5,0.2],[0.7,0.8],"String",str);
hold off

fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<m^2>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f\n",...
    means_m(2),errors_m(2),a(1),err_a(1),a2(1),err_a2(1));
fprintf(resultsFile,"Pendenza a&r: %f +- %f\nPendenza without a&r: %f +- %f\n",...
    a(2),err_a(2),a2(2),err_a2(2));
fprintf(resultsFile,"chi^2 a&r = %f\nchi^2 without a&r = %f\n",chisqr,chisqr2);


% <|m|>
figure()
hold on
for i = 1:nMeasures
    errorbar(deltaTau.^2,means1(:,1),errors1(:,1),'ob');
    errorbar(deltaTau.^2,means2(:,1),errors2(:,1),'or');
end
yline(means_m(1),'k','metropolis','Linewidth',1);
yline(means_m(1)+errors_m(1),':k','Linewidth',1);
yline(means_m(1)-errors_m(1),':k','Linewidth',1);
title('\langle|m|\rangle');
xlabel('(\delta\tau)^2');
ylabel('\langle|m|\rangle');

% Fit con svdfit
%fit a&r
inizio = 1;
x = deltaTau(inizio:end).^2;
y1 = means1(inizio:end,1);
err_y1 = errors1(inizio:end,1);
[a,err_a,covmat,chisqr] = svdfit(x,y1,err_y1,1);
%text(0,0,"Intercetta: %f \pm %f",a(1),err_a(1));
xq = linspace(0,deltaTau(inizio)^2);
[yq,err_yq] = svdpolyval(a,xq,covmat);
plot(xq,yq,'b-');
plot(xq,yq+2*err_yq,'b--',xq,yq-2*err_yq,'b--')
%fit without a&r
y2 = means2(inizio:end,1);
err_y2 = errors2(inizio:end,1);
[a2,err_a2,covmat2,chisqr2] = svdfit(x,y2,err_y2,1);
str = sprintf("chi^2 a&r = %f\nchi^2 without a&r = %f",chisqr,chisqr2);
annotation("textbox",[.2 .1 .2 .2],"String",str,"FitBoxToText","on");
%text(0,0,"Intercetta: %f \pm %f",a2(1),err_a2(1));
[yq2,err_yq2] = svdpolyval(a2,xq,covmat2);
plot(xq,yq2,'r-');
plot(xq,yq2+2*err_yq2,'r--',xq,yq2-2*err_yq2,'r--')
legend('Data with a&r','Data without a&r','Location','east');
%{
legend('Data with a&r','Data without a&r','Linear Fit a&r',...
    '95% Prediction Interval a&r','Linear Fit without a&r',...
    '95% Prediction Interval without a&r','Location','east');
%}
str = sprintf("Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f",...
    means_m(1),errors_m(1),a(1),err_a(1),a2(1),err_a2(1));
annotation('textarrow',[0.5,0.2],[0.7,0.8],"String",str);
hold off

fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<|m|>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f\n",...
    means_m(1),errors_m(1),a(1),err_a(1),a2(1),err_a2(1));
fprintf(resultsFile,"Pendenza a&r: %f +- %f\nPendenza without a&r: %f +- %f\n",...
    a(2),err_a(2),a2(2),err_a2(2));
fprintf(resultsFile,"chi^2 a&r = %f\nchi^2 without a&r = %f\n",chisqr,chisqr2);

% suscettibilità
figure()
hold on
for i = 1:nMeasures
    errorbar(deltaTau.^2,means1(:,3),errors1(:,3),'ob');
    errorbar(deltaTau.^2,means2(:,3),errors2(:,3),'or');
end
yline(means_m(3),'k','metropolis','Linewidth',1);
yline(means_m(3)+errors_m(3),':k','Linewidth',1);
yline(means_m(3)-errors_m(3),':k','Linewidth',1);
title('\langle\chi\rangle');
xlabel('(\delta\tau)^2');
ylabel('\langle\chi\rangle');

% Fit con svdfit
%fit a&r
inizio = 1;
x = deltaTau(inizio:end).^2;
y1 = means1(inizio:end,3);
err_y1 = errors1(inizio:end,3);
[a,err_a,covmat,chisqr] = svdfit(x,y1,err_y1,1);
%text(0,0,"Intercetta: %f \pm %f",a(1),err_a(1));
xq = linspace(0,deltaTau(inizio)^2);
[yq,err_yq] = svdpolyval(a,xq,covmat);
plot(xq,yq,'b-');
plot(xq,yq+2*err_yq,'b--',xq,yq-2*err_yq,'b--')
%fit without a&r
y2 = means2(inizio:end,3);
err_y2 = errors2(inizio:end,3);
[a2,err_a2,covmat2,chisqr2] = svdfit(x,y2,err_y2,1);
str = sprintf("chi^2 a&r = %f\nchi^2 without a&r = %f",chisqr,chisqr2);
annotation("textbox",[.2 .1 .2 .2],"String",str,"FitBoxToText","on");
%text(0,0,"Intercetta: %f \pm %f",a2(1),err_a2(1));
[yq2,err_yq2] = svdpolyval(a2,xq,covmat2);
plot(xq,yq2,'r-');
plot(xq,yq2+2*err_yq2,'r--',xq,yq2-2*err_yq2,'r--')
legend('Data with a&r','Data without a&r','Location','east');
%{
legend('Data with a&r','Data without a&r','Linear Fit a&r',...
    '95% Prediction Interval a&r','Linear Fit without a&r',...
    '95% Prediction Interval without a&r','Location','east');
%}
str = sprintf("Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f",...
    means_m(3),errors_m(3),a(1),err_a(1),a2(1),err_a2(1));
annotation('textarrow',[0.5,0.2],[0.62,0.8],"String",str);
hold off

fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<Suscettibilità>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f\n",...
    means_m(3),errors_m(3),a(1),err_a(1),a2(1),err_a2(1));
fprintf(resultsFile,"Pendenza a&r: %f +- %f\nPendenza without a&r: %f +- %f\n",...
    a(2),err_a(2),a2(2),err_a2(2));
fprintf(resultsFile,"chi^2 a&r = %f\nchi^2 without a&r = %f\n",chisqr,chisqr2);

% Binder cumulant
figure()
hold on
for i = 1:nMeasures
    errorbar(deltaTau.^2,means1(:,4),errors1(:,4),'ob');
    errorbar(deltaTau.^2,means2(:,4),errors2(:,4),'or');
end
yline(means_m(4),'k','metropolis','Linewidth',1);
yline(means_m(4)+errors_m(4),':k','Linewidth',1);
yline(means_m(4)-errors_m(4),':k','Linewidth',1);
title('\langleB\rangle');
xlabel('(\delta\tau)^2');
ylabel('\langleB\rangle');

% Fit con svdfit
%fit a&r
inizio = 1;
x = deltaTau(inizio:end).^2;
y1 = means1(inizio:end,4);
err_y1 = errors1(inizio:end,4);
[a,err_a,covmat,chisqr] = svdfit(x,y1,err_y1,1);
%text(0,0,"Intercetta: %f \pm %f",a(1),err_a(1));
xq = linspace(0,deltaTau(inizio)^2);
[yq,err_yq] = svdpolyval(a,xq,covmat);
plot(xq,yq,'b-');
plot(xq,yq+2*err_yq,'b--',xq,yq-2*err_yq,'b--')
%fit without a&r
y2 = means2(inizio:end,4);
err_y2 = errors2(inizio:end,4);
[a2,err_a2,covmat2,chisqr2] = svdfit(x,y2,err_y2,1);
str = sprintf("chi^2 a&r = %f\nchi^2 without a&r = %f",chisqr,chisqr2);
annotation("textbox",[.35 .58 .2 .2],"String",str,"FitBoxToText","on");
%text(0,0,"Intercetta: %f \pm %f",a2(1),err_a2(1));
[yq2,err_yq2] = svdpolyval(a2,xq,covmat2);
plot(xq,yq2,'r-');
plot(xq,yq2+2*err_yq2,'r--',xq,yq2-2*err_yq2,'r--')
legend('Data with a&r','Data without a&r','Location','east');
%{
legend('Data with a&r','Data without a&r','Linear Fit a&r',...
    '95% Prediction Interval a&r','Linear Fit without a&r',...
    '95% Prediction Interval without a&r','Location','east');
%}
str = sprintf("Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f",...
    means_m(4),errors_m(4),a(1),err_a(1),a2(1),err_a2(1));
annotation('textarrow',[0.33,0.2],[0.8,0.3],"String",str);
hold off

fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"<B>\n");
fprintf(resultsFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(resultsFile,"Intercetta metropolis: %f +- %f\nIntercetta a&r: %f +- %f\nIntercetta without a&r: %f +- %f\n",...
    means_m(4),errors_m(4),a(1),err_a(1),a2(1),err_a2(1));
fprintf(resultsFile,"Pendenza a&r: %f +- %f\nPendenza without a&r: %f +- %f\n",...
    a(2),err_a(2),a2(2),err_a2(2));
fprintf(resultsFile,"chi^2 a&r = %f\nchi^2 without a&r = %f\n",chisqr,chisqr2);

fclose(resultsFile);

end

function [observables,errors] = analyze(dati,V)
% Calcolo medie (valori di aspettazione) di |M|, M^2, M^4
N = size(dati,1);
nObs = size(dati,2);
means = mean(dati(:,1:nObs));
% jackknife vectors for primary variables: |M|, M^2, M^4
primary_jacks = means(1:3)-(dati(1:N,1:3)-means(1:3))/(N-1);
primary_errors = sqrt(var(primary_jacks)*((N-1)^2)/N);
% jackknife for derived variables: susceptibility, binder cumulant
susc_jack = susceptibility(primary_jacks(:,1),primary_jacks(:,2),V);
susc_error = sqrt(var(susc_jack)*((N-1)^2)/N);
bc_jack = binderCumulant(primary_jacks(:,2),primary_jacks(:,3));
bc_error = sqrt(var(bc_jack)*((N-1)^2)/N);
% Summary
observables = [means(1:2)/V,susceptibility(means(1),means(2),V),binderCumulant(means(2),means(3))];
errors = [primary_errors(1:2)/V,susc_error,bc_error];
format long g;

function s = susceptibility (absM,M2,V)
s = (M2-absM.^2)/V;
end
function b = binderCumulant (M2,M4)
b = M4./(M2.^2);
end

end

%{
%fit a&r
inizio = 1;
x = deltaTau(inizio:end).^2;
y1 = means1(inizio:end,2);
erry1 = errors1(inizio:end,2);
[yp,S] = wpolyfit(x,y1,erry1,1);
xq = linspace(0,deltaTau(inizio)^2);
[yq,delta] = polyval(yp,xq,S);
plot(xq,yq,'b-');
plot(xq,yq+2*delta,'b--',xq,yq-2*delta,'b--')
%fit without a&r
x = deltaTau(inizio:end).^2;
y2 = means2(inizio:end,2);
erry2 = errors2(inizio:end,2);
[yp,S] = wpolyfit(x,y2,erry2,1);
xq = linspace(0,deltaTau(inizio)^2);
[yq,delta] = polyval(yp,xq,S);
plot(xq,yq,'r-');
plot(xq,yq+2*delta,'r--',xq,yq-2*delta,'r--')
legend('Data with a&r','Data without a&r','Linear Fit a&r',...
    '95% Prediction Interval a&r','Linear Fit without a&r',...
    '95% Prediction Interval without a&r');
hold off
%}
