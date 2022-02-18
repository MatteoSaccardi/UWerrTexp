import numpy as np
import scipy
import matplotlib.pyplot as plt
from UWerrTexp import UWerrTexp
from svdfit import svdfit
from svdpolyval import svdpolyval
from analyze import susceptibility, binderCumulant, computeObservablesLK

D = 3
lambdaa = 1.145
kappa = 0.188
L = 8
V = L**D

dataPath = "phi4_HMC/main/DATI/"

n_naccu = 14
naccu = 2**np.arange(n_naccu)
T = 1000
ntraj = T*naccu

plt.figure(1)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\langle|m|\rangle$')
plt.title(r'Binning with jackknife VS $\Gamma$-method for $\langle|m|\rangle$.'+'\n\
     Same statistics plots')

plt.figure(2)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\langle m^2 \rangle$')
plt.title(r'Binning with jackknife VS $\Gamma$-method for $\langle m^2 \rangle$.'+'\n\
     Same statistics plots')

plt.figure(3)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\langle \chi \rangle$')
plt.title(r'Binning with jackknife VS $\Gamma$-method for $\langle \chi \rangle$.'+'\n\
     Same statistics plots')

plt.figure(4)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\langle B \rangle$')
plt.title(r'Binning with jackknife VS $\Gamma$-method for $\langle B \rangle$.'+'\n\
     Same statistics plots')

plt.figure(5)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\sigma_{\langle|m|\rangle}$')
plt.title(r'Error by binning with jackknife VS $\Gamma$-method for $\langle|m|\rangle$.'+'\n\
     Same statistics plots')

plt.figure(6)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\sigma_{\langle m^2 \rangle}$')
plt.title(r'Error by binning with jackknife VS $\Gamma$-method for $\langle m^2 \rangle$.'+'\n\
     Same statistics plots')

plt.figure(7)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\sigma_{\langle \chi \rangle}$')
plt.title(r'Error by binning with jackknife VS $\Gamma$-method for $\langle \chi \rangle$.'+'\n\
     Same statistics plots')

plt.figure(8)
plt.xlabel(r'$log_2(naccu)$')
plt.ylabel(r'$\sigma_{\langle B \rangle}$')
plt.title(r'Error by binning with jackknife VS $\Gamma$-method for $\langle B \rangle$.'+'\n\
     Same statistics plots')

OBS = np.zeros((n_naccu,4))
ERR = np.zeros((n_naccu,4))
OBS_GAMMA = np.zeros((n_naccu,4))
ERR_GAMMA = np.zeros((n_naccu,4))
ERR_ERR_GAMMA = np.zeros((n_naccu,4))

for j in range(n_naccu):
    dataFile = open(dataPath+"ntraj_varying/naccu{}_ntraj{}.txt".format(2**j,1000*2**j),'r')

    data = np.zeros((T,5)) # index for (number of measurements, number of observables)
    counter = 0
    for i,line in enumerate(dataFile):
        if len(line.split()) == 0:
            continue
        if not line.startswith("%"):
            data[[counter]] = np.fromstring(line,sep=' \t ')
            counter += 1
    data = np.delete(data, obj = 0, axis = 1) # non siamo interessati ai conteggi
    data = np.delete(data, obj = 0, axis = 1) # non siamo interessati a <M>

    OBS[j,:], ERR[j,:] = computeObservablesLK(data,V) # jackknife per analisi dati

    OBS_GAMMA[j,0],ERR_GAMMA[j,0],ERR_ERR_GAMMA[j,0], \
        tauint,dtauint,_,wopt,gammaFbb,drho = UWerrTexp(data,Name=0,Quantity = lambda x: x[0]/V)
    print("tauint = {} +- {}, wopt = {}".format(tauint,dtauint,wopt))
    OBS_GAMMA[j,1],ERR_GAMMA[j,1],ERR_ERR_GAMMA[j,1], \
        tauint,dtauint,_,wopt,gammaFbb,drho = UWerrTexp(data,Name=0,Quantity = lambda x: x[1]/V)
    #print("tauint = {} +- {}, wopt = {}".format(tauint,dtauint,wopt))
    OBS_GAMMA[j,2],ERR_GAMMA[j,2],ERR_ERR_GAMMA[j,2], \
        tauint,dtauint,_,wopt,gammaFbb,drho = UWerrTexp(data,Name=0,Quantity = lambda x: (x[1]-x[0]**2)/V)
    #print("tauint = {} +- {}, wopt = {}".format(tauint,dtauint,wopt))
    OBS_GAMMA[j,3],ERR_GAMMA[j,3],ERR_ERR_GAMMA[j,3], \
        tauint,dtauint,_,wopt,gammaFbb,drho = UWerrTexp(data,Name=0,Quantity = lambda x: x[2]/(x[1]**2))
    #print("tauint = {} +- {}, wopt = {}".format(tauint,dtauint,wopt))


for i in range(4):
    plt.figure(i+1)
    plt.errorbar(np.log2(naccu), OBS[:,i], yerr = ERR[:,i], marker='o', \
        linewidth=1, linestyle = '--', label = 'Binning method, different runs')
    plt.errorbar(np.log2(naccu), OBS_GAMMA[:,i], yerr = ERR_GAMMA[:,i], color = 'r', \
        linestyle = '--', label = r'$\Gamma-method$')
    plt.legend()

    plt.figure(i+5)
    plt.plot(np.log2(naccu), ERR[:,i], marker = 'o',linewidth=1, linestyle = '--', \
        label = 'Binning method, different runs')
    plt.errorbar(np.log2(naccu), ERR_GAMMA[:,i], yerr = ERR_ERR_GAMMA[:,i], color = 'r', \
        linestyle = '--', label = r'$\Gamma-method$')
    plt.legend()


plt.show()
