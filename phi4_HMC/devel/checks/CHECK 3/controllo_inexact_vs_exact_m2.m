data = load('inexact_vs_exact.txt','\t');
data(:,1) = []; %non ci interessa il conteggio
data(:,1) = []; %non ci interessa m
data(:,1) = []; %non ci interessa |m|
data(:,2) = []; %non ci interessa m^4
V = 64;
data = data./V;
N = 1000; % guarda alla prima riga del file inexact_vs_exact.txt
% guarda nella funzione inexact_vs_exact in hmc.c per il range di deltaTau
deltaTau = 1./(5+3.*(0:9))';
nMeasures = length(data)/N/2;

means_ar = zeros(nMeasures,1);
errors_ar = zeros(nMeasures,1);
means_no = zeros(nMeasures,1);
errors_no = zeros(nMeasures,1);

%1
m1_ar = data(1:N);
data(1:N) = [];
m1_no = data(1:N);
data(1:N) = [];
[means_ar(1),errors_ar(1)] = analyze(m1_ar);
[means_no(1),errors_no(1)] = analyze(m1_no);

%2
m2_ar = data(1:N);
data(1:N) = [];
m2_no = data(1:N);
data(1:N) = [];
[means_ar(2),errors_ar(2)] = analyze(m2_ar);
[means_no(2),errors_no(2)] = analyze(m2_no);
%3
m3_ar = data(1:N);
data(1:N) = [];
m3_no = data(1:N);
data(1:N) = [];
[means_ar(3),errors_ar(3)] = analyze(m3_ar);
[means_no(3),errors_no(3)] = analyze(m3_no);
%4
m4_ar = data(1:N);
data(1:N) = [];
m4_no = data(1:N);
data(1:N) = [];
[means_ar(4),errors_ar(4)] = analyze(m4_ar);
[means_no(4),errors_no(4)] = analyze(m4_no);
%5
m5_ar = data(1:N);
data(1:N) = [];
m5_no = data(1:N);
data(1:N) = [];
[means_ar(5),errors_ar(5)] = analyze(m5_ar);
[means_no(5),errors_no(5)] = analyze(m5_no);
%6
m6_ar = data(1:N);
data(1:N) = [];
m6_no = data(1:N);
data(1:N) = [];
[means_ar(6),errors_ar(6)] = analyze(m6_ar);
[means_no(6),errors_no(6)] = analyze(m6_no);
%7
m7_ar = data(1:N);
data(1:N) = [];
m7_no = data(1:N);
data(1:N) = [];
[means_ar(7),errors_ar(7)] = analyze(m7_ar);
[means_no(7),errors_no(7)] = analyze(m7_no);
%8
m8_ar = data(1:N);
data(1:N) = [];
m8_no = data(1:N);
data(1:N) = [];
[means_ar(8),errors_ar(8)] = analyze(m8_ar);
[means_no(8),errors_no(8)] = analyze(m8_no);
%9
m9_ar = data(1:N);
data(1:N) = [];
m9_no = data(1:N);
data(1:N) = [];
[means_ar(9),errors_ar(9)] = analyze(m9_ar);
[means_no(9),errors_no(9)] = analyze(m9_no);
%10
m10_ar = data(1:N);
data(1:N) = [];
m10_no = data(1:N);
data(1:N) = [];
[means_ar(10),errors_ar(10)] = analyze(m10_ar);
[means_no(10),errors_no(10)] = analyze(m10_no);

clear data;

figure
errorbar(deltaTau,means_ar,errors_ar,'bo');
hold on
errorbar(deltaTau,means_no,errors_no,'ro');
legend('with a&r','without a&r');
title('Double check \langlem^2\rangle');
xlabel('\delta\tau');
ylabel('\langlem^2\rangle');
hold off


function [means,errors] = analyze(data)
means = mean(data);
errors = sqrt(cov(data)/length(data));
end