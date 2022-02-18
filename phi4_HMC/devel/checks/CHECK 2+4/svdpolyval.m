function [y,err_y] = svdpolyval(p,x,covmat)
%SVDPOLYVAL Evaluate polynomial.
%   Y = SVDPOLYVAL(P,X) returns the value of a polynomial P evaluated at X.
%   P is a vector of length N+1 whose elements are the coefficients of the
%   polynomial in ascending powers.
%
%       Y = P(1) + P(2)*X + ... + P(N)*X^(N-1) + P(N+1)*X^N
%
%   If X is a matrix or vector, the polynomial is evaluated at all
%   points in X.  See POLYVALM for evaluation in a matrix sense.
%
%   [Y,ERR_Y] = SVDPOLYVAL(P,X,COVMAT) uses the optional output covariance 
%   matrix covmat created by SVDFIT to generate prediction error estimates
%   ERR_Y. ERR_Y is an estimate of the standard deviation of the error in
%   predicting a future observation at X by P(X).

% Check input is a vector
if ~(isvector(p) || isempty(p))
    error(message('MATLAB:svdpolyval:InvalidP'));
end

nc = length(p);
if isscalar(x) && nargin < 3 && nc > 0 && isfinite(x) && all(isfinite(p(:)))
    % Make it scream for scalar x.  Polynomial evaluation can be
    % implemented as a recursive digital filter.
    y = filter(1,[1 -x],p);
    y = y(nc);
    return
end

siz_x = size(x);

y = zeros(siz_x, superiorfloat(x,p));
for i = 0:nc-1
        y = y + p(i+1)*x.^i;
end

if nargout > 1
    if nargin < 3 || isempty(covmat)
        error(message('MATLAB:svdpolyval:RequiresCovmat'));
    end

    % Since the polynomial is linear in the coefficients, its derivatives
    % respect to them will form a vector like A = (1 x x^2 ... x^n)'.
    A = ones(1,max(size(x)));
    for i = 1:nc-1
        A = [A;x.^i];
    end
    % This is the vector we need to multiply the covariance matrix from the
    % left and the right.
    err_y = diag(A'*covmat*A)';
    err_y = sqrt(err_y);
end

