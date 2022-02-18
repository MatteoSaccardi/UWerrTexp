import numpy as np

def svdpolyval(p,x,covmat):
    '''
    SVDPOLYVAL Evaluate polynomial.
    Y = SVDPOLYVAL(P,X) returns the value of a polynomial P evaluated at X.
    P is a vector of length N+1 whose elements are the coefficients of the
    polynomial in ascending powers.

        Y = P(1) + P(2)*X + ... + P(N)*X^(N-1) + P(N+1)*X^N

    If X is a matrix or vector, the polynomial is evaluated at all
    points in X.  See POLYVALM for evaluation in a matrix sense.

    [Y,ERR_Y] = SVDPOLYVAL(P,X,COVMAT) uses the optional output covariance 
    matrix covmat created by SVDFIT to generate prediction error estimates
    ERR_Y. ERR_Y is an estimate of the standard deviation of the error in
    predicting a future observation at X by P(X).
    '''

    nc = len(p)
    siz_x = len(x)
    
    y = np.zeros_like(x)
    for i in range(nc):
        y = y + p[i]*(x**i)

    # Since the polynomial is linear in the coefficients, its derivatives
    # respect to them will form a vector like A = (1 x x^2 ... x^n)'.
    A = np.ones(siz_x)
    for i in range(1,nc):
        A = np.append(A,x**i)
        A = A.reshape((-1,siz_x))
    # This is the vector we need to multiply the covariance matrix from the
    # left and the right.
    err_y = np.diag(A.T@covmat@A).T
    err_y = np.sqrt(err_y)
    
    return y,err_y
