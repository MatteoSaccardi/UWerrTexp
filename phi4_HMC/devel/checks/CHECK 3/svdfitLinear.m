function [a,sa,covmat,chisqr] = svdfitLinear (x,y,sy)
%SVDFITLINEAR computa i parametri a(1), a(2) con i rispettivi errori sa(1),
%   sa(2), la matrice di covarianza e il chi quadro del fit lineare
%   esplicito dei punti (x,y), dove le y sono soggette alle incertezze sy.
%   La procedura è implementata a partire dal libro "Numerical recipes in
%   C", sezione 15.2 "Fitting Data to a Straight Line".
%   L'interpolazione è della forma y = a(1) + a(2)*x

a = zeros(2,1);
sa = zeros(2,1);
covmat = zeros(2,2);

S = sum(1./(sy.^2));
Sx = sum(x./(sy.^2));
Sy = sum(y./(sy.^2));
Sxx = sum(x.^2./(sy.^2));
Sxy = sum(x.*y./(sy.^2));
D = S*Sxx-Sx^2;
a(1) = (Sxx*Sy-Sx*Sxy)/D;
a(2) = (S*Sxy-Sx*Sy)/D;
covab = -Sx/D;
covmat(1,1) = Sxx/D;
covmat(1,2) = covab;
covmat(2,1) = covab;
covmat(2,2) = S/D;
sa = sqrt(diag(covmat));
chisqr = sum(((y-a(1)-a(2)*x)./sy).^2)/(length(x)-2); % chi quadro ridotto

%rab = -Sx/sqrt(S*Sxx);
%Q = gammaq(length(x)/2-1,chisqr/2);
%r = sum((x-mean(x).*(y-mean(y))))/sqrt((sum((x-mean(x)).^2))*(sum((y-mean(y)).^2)));
%NVar = sum(((y-mean(y))./sy).^2);
%chisqr = (1-r^2)*NVar;

end