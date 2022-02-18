import numpy as np
import matplotlib.pyplot as plt

def autocorrPlot (nameFile = 'provaHMC_naccu100.txt'):
    
    dataFile = open(nameFile,'r')
    #num_lines = np.sum(1 for line in dataFile)
    nData = 0
    data = np.array([])
    for i,line in enumerate(dataFile):
        if len(line.split()) == 0:
            continue
        if not line.startswith("%"):
            data = np.append(data,np.fromstring(line,sep=' \t '))
            nData += 1
    data = data.reshape((nData,-1))
    
    # Import data
    data = np.delete(data,0,1) # drop out the counters in the first column
    #m1 = data[:,0]
    #absMag = data[:,1]
    #M2 = data[:,2]
    #M4 = data[:,3]
    
    means = np.mean(data,axis = 0)
    #meanMag1 = np.mean(m1) = means[0]
    #meanAbsMag = np.mean(absMag) = means[1]
    #meanM2 = np.mean(M2) = means[2]
    #meanM4 = np.mean(M4) = means[3]

    for i in range(len(means)):
        # Set iniziale grafico
        fig = plt.figure(i+1)
        #ax = fig.add_subplot(111)
        #ax.set_aspect('equal', adjustable='box')
        if i == 0: titolo = "autocorrelation plot for $m$"
        elif i == 1: titolo = "autocorrelation plot for $|m|$"
        elif i == 2: titolo = "autocorrelation plot for $m^2$"
        elif i == 3: titolo = "autocorrelation plot for $m^4$"
        else: print('Neglect other observables')
        plt.title(titolo)
        
        N = nData
        nMax = int(N/10)
        points = np.ones(nMax)
        
        for t in range(nMax):
            gamma = 0
            for j in range(N-t):
                gamma = gamma + (data[j][i]-means[i])*(data[j+t][i]-means[i])
            points[t] = gamma/(N-t-1)

        t = np.arange(1,nMax+1)
        plt.scatter(t,points,c='r',marker='.')
        plt.show()

    

#autocorrPlot('phi4_HMC/main/DATI/ntraj_fixed/naccu1.txt')
