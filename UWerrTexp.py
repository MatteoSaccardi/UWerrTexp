import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
import scipy

def nextpow2(N):
    """ Function for finding the next power of 2 """
    if (type(N) == int) or (type(N) == float):
        n = 1
        exponent = 0
        while n < N:
            n *= 2
            exponent += 1
    else:
        n = np.ones(N.shape)
        exponent = np.zeros(N.shape)
        while (n < N).any():
            n[n < N] *= 2
            exponent[n < N] += 1
    return exponent

def autoCF(data):
    ''' Autocorrelation function of data array. '''
    n = data.shape[0]
    siz = 2**nextpow2(3*n)
    fdata = np.fft.fft(data,siz)
    afdata = fdata*fdata.conjugate()
    g = np.real(np.fft.ifft(afdata))
    g = g[range(n)]
    return g

#----------------------------------------------------------------------------

def computeGamma(delpro, Nr, Tmax):
    gamma = np.zeros(Tmax+1)
    i0 = 0
    for r in range(len(Nr)):
        i1 = i0+Nr[r]
        tmp = autoCF(delpro[range(i0,i1)])
        gamma = gamma + tmp[range(Tmax+1)]
        i0 = i0+Nr[r]
    gamma = gamma[range(Tmax+1)]/(len(delpro) - r*np.arange(Tmax+1).T)
    return gamma

#----------------------------------------------------------------------------

def project(Data,evalQ):
  
  N, Nalpha = Data.shape
  fgrad = np.zeros(Nalpha)
  
  h = np.std(Data, axis = 0)*np.sqrt((N-1)/N)/np.sqrt(N)/2 # spread for num.derivative=sigma_i/sqrt(N)/2
  v = np.mean(Data, axis = 0)
  dv = np.copy(v)
  
  for a in range(Nalpha):
    if h[a] == 0: # Data for this observable do not fluctuate
      fgrad[a] = 0
    else:
      dv[a] = v[a] + h[a]
      fgrad[a] = evalQ(dv)
      dv[a] = v[a] - h[a]
      fgrad[a] = fgrad[a] - evalQ(dv)
      dv[a] = v[a]
      fgrad[a] = fgrad[a]/(2*h[a])
  
  # projected deviations (on total statistics)
  delpro = np.dot(Data,fgrad) - np.dot(v,fgrad)
  return delpro

#--------------------------------------------------------------------------

def findWopt(Gamma, Tmax, S, N):
  
  rint = np.cumsum(Gamma[range(1,Tmax+1)]/Gamma[0])
  tauW = S/(np.log(np.abs((rint+1)/rint)))
  tauW[rint <= 0] = np.finfo(float).eps
  t = np.arange(Tmax)+1
  gW = np.exp(-t/tauW)-tauW/np.sqrt(t*N)
  a = (gW <= 0)
  w = 0
  while (a[w] == False):
    w += 1
    if (w == len(gW)):
      warnings.warn('UWerrTexp:FailedCond \n Windowing condition failed up to W = {}'.format(Tmax))
      break
  return w

#----------------------------------------------------------------------------

def computeDrho(rho, wopt, N):
  ''' changed subleading terms (prev. versions!) here:
  contruct errors acc. to hep-lat/0409106 eq. (E.11)
  pad zeros to simplify summation: '''
  L = len(rho)
  paddedsize = 2*L + wopt
  rho = np.append(rho,np.zeros(paddedsize-L))
  drho = np.zeros(L)
  T = np.arange(L)
  BEG = np.zeros(L)
  for n in range(L):
    BEG[n] = np.max((1,(T-wopt)[n]))
  if type(BEG) == int or type(BEG) == float:
    BEG = np.array([BEG])
  
  for t in T: #0:L-1
    k = np.arange(int(BEG[t]),t+wopt) # summation range
    drho[t] = np.sum((rho[k+t]+rho[abs(k-t)]-2*rho[t]*rho[k])**2)
    
  drho = np.sqrt(drho/N)
  return drho

#----------------------------------------------------------------------------

def computeUpperBound(Rho, DRho, Texp, Nsigma, Wsmall):

  Wmax = len(Rho)
  rhosum = np.cumsum(Rho) - 0.5
  
  t = 0
  for x in (Rho - Nsigma*DRho < 0):
    if x == False:
      t += 1
    else:
      break
                                 
  if t == Wmax:
    warnings.warn('UWerrTexp:WtooSmall \n Could not meet {} sigma criterion. \n \
      Setting Wupper = {} (Wmax)'.forma(Nsigma,Wmax))
    t = Wmax
  else:
    try: # try to get the window closest to Nsigma from zero
      myfun = lambda T,N: np.abs(np.abs(Rho[T]) - N*DRho[T])
      M = np.min(myfun([t,t-1],Nsigma)) # compare t and t-1, select the smallest
      IDX = np.where(myfun([t,t-1],Nsigma) == M)[0][0]
      t = t - (IDX-1)
    except Exception as e:
      print(e)
      t = 1
  
  if t <= 1:
    t = 2 # always get a window of size at least one
    
  if (t <= Wsmall) or (Rho[t] < 0): # small window case
    t1 = t
    t2 = t+1
    tmp1 = rhosum[t1] + np.max((Rho[t1],2*DRho[t1]))*Texp
    tmp2 = rhosum[t2] + np.max((Rho[t2],2*DRho[t2]))*Texp
    
    tauintu = np.min((tmp1,tmp2))
    
    if IDX == 2:
      warnings.warn('UWerrTexp:Changingwopt2 \n Estimating tau_int^u at {},\
        instead of wopt(2) = {}\n'.format(t,t-1))
      t = t+1
  else:
    tauintu = rhosum[t] + Rho[t]*Texp
  
  tau  = tauintu
  w = t-1 # w is in physical notation (0...L-1)
  
  return tau,w

#----------------------------------------------------------------------------

def is_numeric(obj):
  '''The easiest way to check if an object is a number is to do arithmethic
  operations (such as add 0) and see if we can get away with it'''
  try:
      obj + 0
      return True
  except TypeError:
      return False

#----------------------------------------------------------------------------

def UWerrTexp(Data, Parm = np.array([1.5, 0, 3, 5, 1, 0]), Nrep = -1, \
  Name = 'NoName', Quantity = 1):
  '''
  UWERRTEXP Autocorrelation-analysis of MC time-series.
  [value,dvalue] = UWerrTexp(Data) Returns the average of the first column of 
  Data and the error of the average (taking into account autocorrelations)
  following the Gamma-method in ``Monte Carlo errors with less errors''
  by Ulli Wolff, hep-lat/0306017
  and the TauExp bias correction following the method in
  ``Critical slowing down and error analysis in lattice QCD simulations''
  arxiv: 1009.5228
  --------------------------------------------------------------------------
  Based on : Ulli Wolff,   Nov. 2006, UWerr Version V6
  --------------------------------------------------------------------------
  Francesco Virotta, Jan. 2011, Version V4
  --------------------------------------------------------------------------
  
  detailed description of the call:
  [value,dvalue,ddvalue,tauint,dtauint,qval,wopt,gammaFbb,drho] = ...
      UWerrTexp(Data,Parm,Nrep,Name,Quantity,P1,P2,....)
      
  omitted or set to [] input values get default value indicated in [D=..]
  
  Data     -- (N x Nalpha) matrix of measured (equilibrium!) data 
             N = total number of measurements
             Nalpha = number of (primary) observables
             if your data are in a different format, consider 
             the MATLAB commands `reshape' and `permute'
             
  Parm     -- A vector of (up to 6) values:
             Parm[0] : guess for the ratio S of tau/tauint 
             if set = 0, absence of autocorrelations is *assumed*
             Parm[1] : Tauexp used to add a tail and give an upper bound to the
             error. If set = 0, performs only 'standard' UWerr analysis
             Parm[2] : Nsigma, determines where to attach the tail (larger
             values attach the tail at earlier times, the tail *wont* be
             attached past the first negative value of the autocorrelation
             function)
             Parm[3] : Wsmall, check for fast decay of gamma (if wopt < Wsmall
             the amplitude for the tail is given by max(rho(wopt+1),2*drho(wopt+1)), 
             only if this results in a smaller tauexp than the one at wopt).
             This amplitude is always used if rho(wopt)<0
             Parm[4] : MDR, Molecular dynamics times R (factor used for
             rescaling plots)
             Parm[5] : If ~=0 outputs both upper bound and lower bound
             values for error, tauint and wopt
             [D = [1.5 0 3 5 1 0]]
  
  Nrep     -- vector [N1 N2 N3 ...] specifying a breakup of the N rows
             of Data into replica of length N1,N2,.. (N=sum(Nrep)!)  [D=[N]]
             The replica distribution is histogrammed and a Q-value 
             (=probability to find this much or more scatter) is given 
             for R >= 2 replica; if you have one history you may set Nrep to 
             artificially split it into >=2 replica to see their distribution                
  
  Name     -- if string: name of observable for titles of generated plots
             if not string: all plots are supressed [D='NoName']

  Quantity -- either:
             scalar function handle (@functionname) for the derived 
             observable F; it has to operate on a row-vector of length 
             Nalpha as first argument; optional parameters P1,P2,... are 
             passed on to this function as 2nd, 3rd, .... argument
          -- or:
             integer between 1 and Nalpha to analyze primary observable [D=1]
  --------------------------------------------------------------------------
  output::
  value    -- estimate for F
  dvalue   -- its error (incl. autocorrelation effects). If param[5]=1,
              returns a vector of 2 values: [lowerbound upperbound]
  ddvalue  -- statistical error of the error (only lowerbound)
  tauint   -- integrated autocorrelation time: If param[5]=1,
              returns a vector of 2 values: [lowerbound upperbound]
  dtauint  -- statistical error of tauint. If param[5]=1,
              returns a vector of 2 values: [lowerbound upperbound], upperbound
              assumes no error on tauexp.
  qval     -- Q-value of the replica distribution if R >= 2
              (goodness of fit to a constant)
  wopt     -- returns the numerical value of the summation window. If
              param[5]=1, returns a vector of 2 values: [lowerbound upperbound]
  gammaFbb -- Autocorrelation function (only up to 2*wopt)
  drho     -- Error on the normalized autocorrelation function
  '''
  #--------------------------------------------------------------------------
  
  # analyze input arguments:
  N, Nalpha = Data.shape # number of measurements and observables
  Dparm = np.array([1.5, 0, 3, 5, 1, 0]) # default value of Parm
  if Nrep == -1:
    Nrep = N
  L = len(Parm)
  DL = len(Dparm)
  if L > DL:
     raise Exception('UWerrTexp:IncorrectArgin \n Incorrect value of Parm')
  Parm = np.append(Parm,Dparm[range(L,DL)])
  
  Stau  = Parm[0]
  Texp  = Parm[1]
  Nsigma = Parm[2]
  Wsmall = Parm[3]
  MDR = Parm[4]
  BOTH = (Parm[5] != 0)
  
  if (Nrep != np.round(Nrep)).any() or (Nrep < 1) or (np.sum(Nrep) != N):
    raise Exception('UWerrTexp:IncorrectArgin \n inconsistent N,Nrep')
    
  if is_numeric(Quantity):
    if np.round(Quantity) != Quantity or Quantity < 1 or Quantity > Nalpha:
      raise Exception('UWerrTexp:IncorrectArgin \n illegal numeric value of Quantity')
    primary = 1 # logical variable
    primaryindex = Quantity
  else:
   primary = 0
   evalQ = lambda x: Quantity(x)
  if type(Nrep) == int or type(Nrep) == float:
    Nrep = np.array([Nrep])
  Nrep = Nrep[:] # make column
  R = len(Nrep) # number of Replica
    
  #--------------------------------------------------------------------------
  # means of primaries:
  
  abb = np.mean(Data, axis = 0) # total mean, primary obs.
  abr = np.zeros((R,Nalpha)) # replicum-mean, primary obs.
  i0 = 0
  for r in range(R):
    i1 = i0-1+Nrep[r]
    abr[r,:] = np.mean(Data[range(i0,i1+1),:], axis = 0)
    i0 = i0+Nrep[r]
    
  # means of derived (or primary, depending on Quantity):
  if primary:
    Fbb = abb[primaryindex]
    Fbr = abr[:,primaryindex]
  else:
    Fbb = evalQ(abb) # total mean, derived obs.
    Fbr = np.zeros(R) # replica-means, derived obs.
    for i in range(R):
      Fbr[i] = evalQ(abr[i,:])
  Fb = np.sum(Fbr*Nrep)/N # weighed mean of replica means
  
  # --------------------------------------------------------------------------
  # form the gradient of f and project fluctuations:
  
  if primary:
    delpro = Data[:,primaryindex] - abb[primaryindex]
  else:
    delpro = project(Data,evalQ)
  
  # --------------------------------------------------------------------------
  # Computation of Gamma, automatic windowing
  # note: W=1,2,3,.. in Matlab-vectorindex <-> W=0,1,2,.. in the paper!
  # note: W=0,1,2,.. in Matlab-vectorindex <-> W=0,1,2,.. in the paper!
  # sick case:
  Fvar = np.mean(delpro**2)
  if Fvar == 0:
    raise Exception('UWerrTexp:NoFluctuations \n No fluctuations!')
  
  # compute Gamma(t), find the optimal W
  if Stau == 0: # no autocorrelations assumed
    wopt = 0
    tmax = 0
  else:
    wopt = 1
    tmax = int(np.floor(np.min(Nrep)/2)) # do not take t larger than this
  
  gammaFbb = computeGamma(delpro,Nrep,tmax)

  if wopt: # compute wopt
    wopt = findWopt(gammaFbb, tmax, Stau, N)
    tmax = np.min((tmax,2*wopt))
  
  gammaFbb = gammaFbb[range(tmax+1)]
  
  CFbbopt = gammaFbb[0] + 2*np.sum(gammaFbb[range(1,wopt+1)]) # first estimate
  if CFbbopt <= 0:
    raise Exception('UWerrTexp:GammaPathological \n \
      Gamma pathological: estimated error^2 < 0')
  
  gammaFbb = gammaFbb+CFbbopt/N # bias in Gamma corrected
  CFbbopt = gammaFbb[0] + 2*np.sum(gammaFbb[range(1,wopt+1)]) # refined estimate
  sigmaF = np.sqrt(CFbbopt/N) # error of F
  rho = gammaFbb/gammaFbb[0] # normalized autocorr.
  rhosum = np.cumsum(rho)-0.5 # tau_int(t)
  
  # ----------------------------------------------------------------------------
  # bias cancellation for the mean value
  
  if R >= 2:
    bF = (Fb-Fbb)/(R-1)
    Fbbc= Fbb - bF
    if np.abs(bF) > sigmaF/4:
      warnings.warn('UWerrTexp:BiasCancel \n A {} sigma bias of the mean has \
         been cancelled.'.format(bF/sigmaF))
    Fbr = Fbr - bF*N/Nrep
    Fb  = Fb - bF*R
  
  # ----------------------------------------------------------------------------
  # answers to be returned: (besides gammaFbb)
   
  value = Fbb
  dvalue = sigmaF
  ddvalue = dvalue*np.sqrt((wopt+0.5)/N)
  tauint = rhosum[wopt]
  dtauint = tauint*2*np.sqrt((wopt-tauint+0.5)/N)
  drho = computeDrho(rho,wopt,N)
  
  if Texp:
    tauintu, wupper = computeUpperBound(rho,drho,Texp,Nsigma,Wsmall)
    CFbboptu = gammaFbb[0]*2*tauintu # upper estimate
    IDX = 1 + BOTH
    dvalue[IDX]  = np.sqrt(CFbboptu/N)
    taudraw = [tauint, tauintu]
    tauint[IDX] = tauintu
    dtauint[IDX] = np.sqrt(dtauint[0]**2+(drho[wupper+1]*Texp)**2) # assumes error on texp=0
    wdraw = [wopt, wupper]
    wopt[IDX] = wupper
  else:
    taudraw = tauint
    wdraw = wopt
  
  # Q-value for replica distribution if R >= 2:
  if R >= 2:
    chisq = np.sum(((Fbr-Fb)**2)*Nrep)/CFbbopt
    qval = 1-scipy.special.gammainc(chisq/2,(R-1)/2)
  else:
    qval  = []
  
  # ----------------------------------------------------------------------------
  # Plotting:
  doplot = [ Stau != 0, primary !=0, R >= 2]
  # make plots, if demanded:
  if type(Name) != 'str': # no plots
    '''for n in range(3):
      if len(f) > 0 and type(f[n]) == 'function':
        close(f[n])'''
    f = []
    return value,dvalue,ddvalue,tauint,dtauint,qval,wopt,gammaFbb,drho
  elif not('f' in locals()) and not('f' in globals()):
    f = [0, 1, 2]
    for n in range(3):
      if doplot[n]:
        f[n] = plt.figure(n)
      else:
        f[n] = np.nan
  '''elif (type(f) != 'function').any():
    f = [0, 1, 2]
    for n in range(3):
      if doplot[n] and type(f[n]) != 'function':
        f[n] = plt.figure(n)'''
  
  # autocorrelation:
  
  if doplot[0]:
    plt.figure(f[1])
    # plot of GammaF(t)/GammaF(0) = rhoF(t)
    # put the right dimensions to everything
    time = MDR*np.arange(tmax+1)
    w = wdraw*MDR
    texp = Texp*MDR
    rhosum = rhosum*MDR
    drhosum = 2*rhosum*np.sqrt( time/(N*MDR) )
    tau = taudraw*MDR
    
    plt.subplots(1,2,1)
    plt.errorbar(time,rho,drho)
    v = [0, time[-1], np.min(np.min(rho-drho)*1.2, axis = 0), 1]
    matplotlib.axes.Axes.set_xlim(v[0],v[1])
    matplotlib.axes.Axes.set_ylim(v[2],v[3])
    plt.plot([0,time[-1]],[0,0],'--') # zero line
    
    if Texp:
      plt.plot([w[0],w[0]], [v[2],v[3]], 'r--') # Bar at woptl
      plt.plot([w[1],w[1]], [v[2],v[3]], 'r--') # Bar at woptu 
    else:
      plt.plot([w, w], [v[2],v[3]],'r-') # Bar at wopt
    plt.title('normalized autocorrelation of {}'.format(Name))
    plt.ylabel('$\rho$')
    # ----------------------------------------------------------------
    # plot of rhosum
    plt.subplots(1,2,2)
    # tauint vs. W:
    plt.errorbar(time, rhosum, drhosum)
    v = [0, time[-1], np.min(rhosum-drhosum)*0.8, np.max(rhosum+drhosum)*1.2]
    
    if Texp:
      rhosumupper  = rhosum + np.max((rho,2*drho))*texp
      drhosumupper = np.sqrt(drhosum**2+(drho*texp)**2)
      v[3] = np.max(( v[3], np.max(rhosumupper[range(1,len(rhosumupper))] + \
        drhosumupper[range(1,len(drhosumupper))])*1.1))
    matplotlib.axes.Axes.set_xlim(v[0],v[1])
    matplotlib.axes.Axes.set_ylim(v[2],v[3])
    
    if Texp:
      plt.errorbar(time,rhosumupper,drhosumupper)
      plt.plot([0, time[-1]], [tau[0], tau[0]], 'r--') # hor. line at estimate
      plt.plot([0, time[-1]], [tau[1], tau[1]], 'r--') # hor. line at estimate
      plt.plot([w[0], w[0]], [v[2],v[3]], 'r--') # Bar at wopt
      plt.plot([w[1], w[1]], [v[2],v[3]], 'r--') # Bar at wopt
    else:
      plt.plot([0, time[-1]], [tau, tau], 'r--') # hor. line at estimate
      plt.plot([w, w], [v[2],v[3]], 'r-') # Bar at wopt
    plt.title('\tau_{int} with statistical errors of {}'.format(Name))
    plt.ylabel('$\tau_{int}$')
    plt.xlabel('W')
    
  # histogram of data, if primary:
  if doplot[1]:
    plt.figure(f[1])
    plt.hist(Data[:,primaryindex], bins = 20)
    plt.title('Distribution of all data for {}'.format(Name))
  
  # histogram of replica values if R >= 2:
  if doplot[2]:
    plt.figure(f[2])
    p = (Fbr-Fb)/(sigmaF*np.sqrt(N/Nrep-1))
    plt.hist(p)
    plt.title('replica distribution (mean=0,var=1) for {} >> Q = {} <<'\
      .format(Name,str(qval)))
    plt.draw()
    if (type(f) != 'function').any():
      f = []
  
  return value,dvalue,ddvalue,tauint,dtauint,qval,wopt,gammaFbb,drho
