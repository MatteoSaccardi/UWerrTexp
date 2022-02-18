function [a,sa,covmat,chisqr] = svdfit(x,y,sy,n)
%
%		[a, sa, covmat, chisqr] = svdfit(x,y,z,n)
%
%	returns the parameters of a nth order fit (a in order
%	of const, x, x^2, x^3, etc) for values of n=0 on up, 
%   along with errors (sa), chi squared and the covariance
%   matrix for the data stored in x, and y (N-length column
%   vectors). Routine uses SVD to solve normal equations.
%
%		first check out inputs
%			x.....

M = n+1;

[N L]=size(x);
if L ~= 1, x=x'; [N L]=size(x); end
if N < M, error('Insufficient number of data for fit'); end

%			and then y.....

[Ny L]=size(y);
if L ~= 1, y=y'; [Ny L]=size(y); end
if N ~= Ny, error('x and y must be same length'); end

%
%	 Now construct design matrix
%

A=[ones(N,1)];
if n>0
	for i=1:n
		A = [A x.^i ];
    end
end

%
%   Put the weights
%

w = 1./sy;
w = w(:);
y = y.*w;

for j = 1:n+1
    A(:,j) = w.*A(:,j);
end
%
%		then do fit
%
[U S V]=svd(A,0);		% Then solve by SVD
W=diag(S);			% Get singular values
Wmin=max(W)*1.e-5;		% find minimum s.v. permissable
for i=1:M			% and invert diagonal
	if W(i) > Wmin
		W(i)=1/W(i);
	    else
		W(i)=0;
	end
end
Si=diag(W);			% and remake matrix
a=V*Si*(U'*y);			% compute coefficient vector
covmat=V*Si^2*V';		% compute covariance matrix
%[N,M] = size(A);        % size of design matrix (row by column)
chisqr=sum((A*a-y).^2)/(N-M);	% calculate chisquare
sa=sqrt(diag(covmat));	% errors from chisqr * diag of cov matrix

end