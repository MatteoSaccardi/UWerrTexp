import numpy as np

def svdfit(x,y,sy,n):
    '''
    a, sa, covmat, chisqr = svdfit(x,y,z,n)
    returns the parameters of a nth order fit (a in order
    of const, x, x^2, x^3, etc) for values of n=0 on up, 
    along with errors (sa), chi squared and the covariance
    matrix for the data stored in x, and y (N-length column
    vectors). Routine uses SVD to solve normal equations.
    '''
    
    # first check out inputs
    
    # x.....
    M = n+1
    if type(x) == int or type(x) == float:
        N = len(x)
    else:
        N = x.shape[0]
        if N == 1:
            x = x.T
            N = x.shape[0]
        if N < M:
            raise Exception('Insufficient number of data for fit')
    
    # and then y.....
    if type(y) == int or type(y) == float:
        Ny = len(y)
    else:
        Ny = y.shape[0]
        if Ny == 1:
            y = y.T
            Ny = y.shape[0]
        if N != Ny:
            raise Exception('x and y must be same length')

    #
    # Now construct design matrix
    #
    A = np.ones(N)
    if n > 0:
        for i in range(1,n+1):
            A = np.append(A,x**i)
            A = A.reshape((-1,N))
    A = A.T

    #
    # Put the weights
    #
    w = 1/sy
    y = y*w
    y = np.array([y]).T
    
    for j in range(n+1):
        A[:,j] = w*A[:,j]
    #
    # then do fit
    #
    U,S,V = np.linalg.svd(A,0) # Then solve by SVD
    W = S # Get singular values
    Wmin = np.max(W)*1.e-5 # find minimum s.v. permissable
    for i in range(M): # and invert diagonal
        if W[i] > Wmin:
            W[i] = 1/W[i]
        else:
            W[i] = 0
    Si = np.diag(W) # and remake matrix
    a = (V.T)@Si@(U.T)@y # compute coefficient vector
    covmat = (V.T)@Si**2@V # compute covariance matrix
    # N,M = A.shape # size of design matrix (row by column)
    chisqr = np.sum((A@a-y)**2)/(N-M) # calculate chisquare

    sa = np.sqrt(np.diag(covmat)) # errors from chisqr * diag of cov matrix

    a = a.reshape(-1) # 1d array s.t. coefficients are a[0], a[1], ...
    
    return a,sa,covmat,chisqr