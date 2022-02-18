import numpy as np
import scipy
import matplotlib.pyplot as plt
from UWerrTexp import UWerrTexp
from svdfit import svdfit
from svdpolyval import svdpolyval

def susceptibility (absM,M2,V):
    return (M2-absM**2)/V

def binderCumulant (M2,M4):
    return M4/(M2**2)

def computeObservablesLK(dati,V):
    '''
    Calcolo delle osservabili
    (1/V)<|M|>  (1/V)<M^2>  <suscettibilità>    <binderCumulant>
    a L, kappa fissati.
    '''

    # Calcolo medie (valori di aspettazione) di |M|, M^2, M^4
    N = dati.shape[0]
    means = np.mean(dati[:,range(3)], axis = 0)
    # jackknife vectors for primary variables: |M|, M^2, M^4
    primary_jacks = means-(dati-means)/(N-1)
    primary_errors = np.sqrt(np.var(primary_jacks, axis = 0)*((N-1)**2)/N)
    # jackknife for derived variables: susceptibility, binder cumulant
    susc_jack = susceptibility(primary_jacks[:,0],primary_jacks[:,1],V)
    susc_error = np.sqrt(np.var(susc_jack)*((N-1)**2)/N)
    bc_jack = binderCumulant(primary_jacks[:,1],primary_jacks[:,2])
    bc_error = np.sqrt(np.var(bc_jack)*((N-1)**2)/N)
    # Summary
    observables = np.array([means[0]/V,means[1]/V,susceptibility(means[0],means[1],V),binderCumulant(means[1],means[2])])
    errors = [primary_errors[0]/V,primary_errors[1]/V,susc_error,bc_error]

    return observables,errors

def analyze (nameFolder = "0", wantchecks = 0, wantkexact = 0):
    '''
    ANALYZE
    Analizziamo i dati raccolti con simulazione.c a valori di kappa e L, V
    variabili. In particolare, le simulazioni sono state fatte con
    lambda = 1.145; D = 3;  ntraj = 1000000; ntherm = 10000; naccu = 10000;
    k = {0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23}
    L/a = {4,6,8,10,12,14,16}
    I dati sono presenti nei vari file 'simulazioneL#' dentro le cartelle
    nella presente directory, che possono essere chiamate in modo esplicito, 
    e sono così distribuiti:
   - su cinque colonne: conteggio, <M>, <|M|>, <M^2>, <M^4>;
   - in nove blocchi, ognuno a k fissato e ognuno contente un numero di
       misure (ciascuna è una media ogni naccu) pari a ntraj/naccu;
   - per ogni file, L è fissato (dal valore '#' nel nome).
   
   Esempi di chiamate:
   analyze % senza check, usando i file nella cartella naccu = 1000
   analyze("naccu = 10000/",1) % con check , usando i file nella cartella
   naccu = 10000
   analyze("naccu = 1000/",1) % con check, usando i file nella cartella naccu = 1000,
   equivale a chiamare analyze("0",1)
   
   Se wantkexact è diverso da 0, si utilizzano anche i dati in
   simulazione_k0188, che dovrebbe essere un valore più preciso per il k
   in cui si verifica la transizione di fase.
   Ad esempio:
   analyze("naccu = 1000/",0,1)
   utilizza il file "simulazione_k01888.txt" per l'interpolazione lineare
   del massimo della susettibilità in funzione di L.
   
   OSS: analyze("naccu = 1000/",0,1) è la chiamata finora migliore per
   l'interpolazione lineare finale
   '''
   
    if nameFolder == "0":
       nameFolder = "naccu = 1000/"
    
    # Definiamo le variabili
    D = 3
    lambdaa = 1.145
    nameFiles = [nameFolder+"simulazioneL4.txt",nameFolder+"simulazioneL6.txt", \
        nameFolder+"simulazioneL8.txt", nameFolder+"simulazioneL10.txt", \
        nameFolder+"simulazioneL12.txt", nameFolder+"simulazioneL14.txt", \
        nameFolder+"simulazioneL16.txt"]
    lati = np.array([4,6,8,10,12,14,16])
    volumes = lati**D
    kappa = np.array([0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23])
    T = 100 # ntraj/naccu, numero di misure per ogni kappa, "periodo"

    '''
    Estraiamo i dati dai file
    Vogliamo inserirli in un'unica variabile data(lati,kappa,1:T,1:3) dove il
    primo indice fa riferimento al valore di L in lati, il secondo al valore
    di k in kappa, il terzo al numero della misura da 1 a T e il quarto
    all'osservabile (da 1 a 3, dato che ci interessano solo le osservabili 
    <|M|>, <M^2>, <M^4>). Alla fine, sarà una variabile 7x9x100x3.
    '''
    data = np.zeros((len(lati),len(kappa),T,3))
    for i1 in range(len(lati)):
        dataFile = open(nameFiles[i1],'r')
        data1 = np.zeros((len(kappa)*T,5))
        counter = 0
        for i,line in enumerate(dataFile):
            if len(line.split()) == 0:
                continue
            if not line.startswith("%"):
                data1[[counter]] = np.fromstring(line,sep=' \t ')
                counter += 1
        data1 = np.delete(data1, obj = 0, axis = 1) # non siamo interessati ai conteggi
        data1 = np.delete(data1, obj = 0, axis = 1) # non siamo interessati a <M>
        for i2 in range(len(kappa)):
            data[i1,i2,:,:] = data1[range(i2*T,(i2+1)*T),:]
    
    '''
    Ora, vogliamo analizzare i dati. Creiamo un unico array corrispondente
    alle quattro osservabili che ci interessano. Il primo indice corrisponde
    all'osservabile, nell'ordine di
    (1/V)<|M|>  (1/V)<M^2>  <suscettibilità>    <binderCumulant>
    Il secondo indice corrisponde al valore di L, il terzo a quello di k.
    Inoltre, ad ogni array associamo anche un array analogo contenente gli errori.
    '''
    OBS = np.zeros((4,len(lati),len(kappa)))
    ERR = np.zeros((4,len(lati),len(kappa)))
    for i1 in range(len(lati)):
        for i2 in range(len(kappa)):
            inputData = data[i1,i2,:].reshape((T,3))
            OBS[:,i1,i2],ERR[:,i1,i2] = computeObservablesLK(inputData,volumes[i1])
    
    # Quindi, effettuiamo il plot delle osservabili al variare di kappa (sulle
    # ascisse) e L (varie linee).

    # <|m|>
    plt.figure(1)
    plt.xlabel(r'$\kappa$')
    titolo = r"$\lambda$ = {}, D = {}".format(lambdaa,D)
    plt.title(titolo)
    plt.ylabel(r'$\langle|m|\rangle$')
    for i1 in range(len(lati)):
        lineName = "L/a = {}".format(lati[i1])
        inputOBS = OBS[0,i1,:].reshape(len(kappa))
        inputERR = ERR[0,i1,:].reshape(len(kappa))
        plt.errorbar(kappa, inputOBS, yerr = inputERR, marker='o', \
            linewidth=1, linestyle = '--', label = lineName)
    plt.legend(loc='upper left')
    plt.draw()

    # <m^2>
    plt.figure(2)
    plt.xlabel(r'$\kappa$')
    titolo = r"$\lambda$ = {}, D = {}".format(lambdaa,D)
    plt.title(titolo)
    plt.ylabel(r'$\langle m^2 \rangle / V$')
    # Nel plot, dividiamo per V per avere quantità indipendenti dal volume a k grandi
    for i1 in range(len(lati)):
        lineName = "L/a = {}".format(lati[i1])
        inputOBS = OBS[1,i1,:].reshape(len(kappa))/volumes[i1]
        inputERR = ERR[1,i1,:].reshape(len(kappa))/volumes[i1]
        plt.errorbar(kappa, inputOBS, yerr = inputERR, marker='o', \
            linewidth=1, linestyle = '--', label = lineName)
    plt.legend(loc='upper left')
    plt.draw()

    # <suscettibilità>
    plt.figure(3)
    plt.xlabel(r'$\kappa$')
    titolo = r"$\lambda$ = {}, D = {}".format(lambdaa,D)
    plt.title(titolo)
    plt.ylabel(r'$\langle\chi\rangle$')
    for i1 in range(len(lati)):
        lineName = "L/a = {}".format(lati[i1])
        inputOBS = OBS[2,i1,:].reshape(len(kappa))
        inputERR = ERR[2,i1,:].reshape(len(kappa))
        plt.errorbar(kappa, inputOBS, yerr = inputERR, marker='o', \
            linewidth=1, linestyle = '--', label = lineName)
    plt.legend(loc='upper left')
    plt.draw()

    # <binderCumulant>
    plt.figure(4)
    plt.xlabel(r'$\kappa$')
    titolo = r"$\lambda$ = {}, D = {}".format(lambdaa,D)
    plt.title(titolo)
    plt.ylabel(r'$\langle B \rangle$')
    for i1 in range(len(lati)):
        lineName = "L/a = {}".format(lati[i1])
        inputOBS = OBS[3,i1,:].reshape(len(kappa))
        inputERR = ERR[3,i1,:].reshape(len(kappa))
        plt.errorbar(kappa, inputOBS, yerr = inputERR, marker='o', \
            linewidth=1, linestyle = '--', label = lineName)
    plt.legend(loc='upper left')
    plt.draw()

    #Possiamo effettuare dei controlli sull'autocorrelazione.
    if wantchecks:
        nSigma = 2
        checkFile = open("checkFile.txt","w")
        obsCHECK = np.zeros((4,len(lati),len(kappa)))
        errorsCHECK = np.zeros((4,len(lati),len(kappa)))
        errorsErrorsCHECK = np.zeros((4,len(lati),len(kappa)))
        # tauint si attesta a circa 0.5 se l'autocorrelazione è assente, mentre
        # è superiore in caso contrario
        checkFile.write("Controllo su tauint > 0.5 entro {} sigma.\n".format(nSigma))
        for i1 in range(len(lati)):
            for i2 in range(len(kappa)):
                inputData = data[i1,i2,:,:].reshape((T,3))
                obsCHECK[0,i1,i2],errorsCHECK[0,i1,i2],errorsErrorsCHECK[0,i1,i2], \
                    tauint,dtauint,_,_,_,_ = UWerrTexp(inputData,Name=0,Quantity = lambda x: x[0]/volumes[i1])
                if (tauint-nSigma*dtauint > 0.5):
                    checkFile.write("tauint per L = {}, k = {}, osservabile 1: {} +- {}\n".format(lati[i1],kappa[i2],tauint,dtauint))
                obsCHECK[1,i1,i2],errorsCHECK[1,i1,i2],errorsErrorsCHECK[1,i1,i2], \
                    tauint,dtauint,_,_,_,_ = UWerrTexp(inputData,Name=0,Quantity = lambda x: x[1]/volumes[i1])
                if (tauint-nSigma*dtauint > 0.5):
                    checkFile.write("tauint per L = {}, k = {}, osservabile 2: {} +- {}\n".format(lati[i1],kappa[i2],tauint,dtauint))
                obsCHECK[2,i1,i2],errorsCHECK[2,i1,i2],errorsErrorsCHECK[2,i1,i2], \
                    tauint,dtauint,_,_,_,_ = UWerrTexp(inputData,Name=0,Quantity = lambda x: (x[1]-x[0]**2)/volumes[i1])
                if (tauint-nSigma*dtauint > 0.5):
                    checkFile.write("tauint per L = {}, k = {}, osservabile 3: {} +- {}\n".format(lati[i1],kappa[i2],tauint,dtauint))
                obsCHECK[3,i1,i2],errorsCHECK[3,i1,i2],errorsErrorsCHECK[3,i1,i2], \
                    tauint,dtauint,_,_,_,_ = UWerrTexp(inputData,Name=0,Quantity = lambda x: x[2]/(x[1]**2))
                if (tauint-nSigma*dtauint > 0.5):
                    checkFile.write("tauint per L = {}, k = {}, osservabile 4: {} +- {}\n".format(lati[i1],kappa[i2],tauint,dtauint))
        checkFile.write("Fine controllo su tauint\n\n")
        # Ci aspettiamo che le medie siano uguali nei due casi.
        diffMeans = np.linalg.norm(OBS.reshape(-1)-obsCHECK.reshape(-1))
        checkFile.write("Le medie differiscono per {}\n\n".format(diffMeans))
        # Controlliamo che gli errori non differiscano in modo significativo
        # da quelli calcolati   
        checkFile.write("Controllo su errori compatibili entro {} sigma.\n".format(nSigma))
        LHS = np.abs(ERR-errorsCHECK)-nSigma*errorsErrorsCHECK
        LHS = LHS.reshape(-1)
        RHS = np.abs(ERR-errorsCHECK)+nSigma*errorsErrorsCHECK
        RHS = RHS.reshape(-1)
        CHECKERR = ERR.reshape(-1)
        errorsCHECK = errorsCHECK.reshape(-1)
        errorsErrorsCHECK = errorsErrorsCHECK.reshape(-1)
        # Se 0 è contenuto tra LHS e RHS, allora ERR-errorsCHECK è compatibile
        # con 0 i.e. i due valori sono tra loro compatibili
        bad = np.where((RHS < 0) + (LHS > 0))
        for i in range(len(bad)):
            checkFile.write("[{} {}] con errore {}\n".format(LHS[bad[i]],RHS[bad[i]],errorsErrorsCHECK[bad[i]]))
            checkFile.write("Errore calcolato: {}; errore UWerrTexp: {}; errore sull'errore: {}\n"\
                .format(CHECKERR[bad[i]],errorsCHECK[bad[i]],errorsErrorsCHECK[bad[i]]))
        checkFile.write("Fine controllo su errori\n")
    '''
    Infine, proviamo ad estrarre l'esponente critico di max(suscettibilità)
    in funzione di L. L'andamento atteso è di
        M := max(chi) ~ L^c,
    quindi:
        M = a * L^c
        log(M) = log(a) + c*log(L)
        y = a(1) + a(2)*x, y = log(M), x = log(L)
    Si osserva che qui il massimo pare essere raggiunto a k = 0.19 (in
    realtà, un valore più preciso sarebbe minore; eventualmente, possiamo 
    generare i dati in simulazione_k.c i.e. simulazioni a k fissato, da fare 
    per ogni L, per effettuare questo test).
    Nel file, sarebbero raccolte le misure fatte a k fissato al variare di
    L. Per ogni L sono effettuate T = 100 misure, organizzate come al solito
    in 3 colonne (5, di cui 3 interessanti).
    '''
    if wantkexact:
        # Consideriamo k = 0.188
        datamax = np.zeros((len(lati),T,3))
        dataFile = open("simulazione_k0188.txt",'r')
        datak = np.zeros((len(lati)*T,5))
        counter = 0
        for i,line in enumerate(dataFile):
            if len(line.split()) == 0:
                continue
            if not line.startswith("%"):
                datak[[counter]] = np.fromstring(line,sep=' \t ')
                counter += 1
        datak = np.delete(datak, obj = 0, axis = 1) # non siamo interessati ai conteggi
        datak = np.delete(datak, obj = 0, axis = 1) # non siamo interessati a <M>
        for i1 in range(len(lati)):
            datamax[i1,:,:] = datak[range(i1*T,(i1+1)*T),:]
        T = 100
        OBSk = np.zeros((4,len(lati)))
        ERRk = np.zeros((4,len(lati)))
        for i1 in range(len(lati)):
            inputData = datamax[i1,:,:].reshape((T,3))
            OBSk[:,i1],ERRk[:,i1] = computeObservablesLK(inputData,volumes[i1])
        M = OBSk[2,:] # Suscettibilità
        ErrM = ERRk[2,:]
    else:
        # Approssimiamo il massimo a k = 0.19, che è il quinto valore dei k
        M = OBS[2,:,4]
        ErrM = ERR[2,:,4]
    
    x = np.log(lati)
    y = np.log(M)
    err_y = ErrM/M
    a,err_a,covmat,chisqr = svdfit(x,y,err_y,1)
    # OSSERVAZIONE: INTERPOLANDO IN MODO QUADRATICO, VIENE MOLTO MEGLIO!
    # i.e. log(M) va come log^2(L)
    # [a,err_a,covmat,chisqr] = svdfit(x,y,err_y,2);
    # Potrebbe però essere un effetto dovuto al fatto che a L grande non cresca
    # in modo giusto poiché questo k non è ancora esattamente quello a cui si
    # verifica la transizione di fase.
    xq = np.linspace(x[0],x[-1],100)
    yq,err_yq = svdpolyval(a,xq,covmat)
    
    # Draw
    plt.figure()
    titolo = "Fit critical exponents with $\lambda$ = {}, D = {}".format(lambdaa,D)
    plt.title(titolo)
    plt.xlabel("log L/a")
    plt.ylabel("log max($\chi$)")
    plt.errorbar(x, y, err_y, marker = '.', markersize = 10, label = 'Data')
    plt.plot(xq,yq,'b-',label='Fit')
    plt.plot(xq,yq+2*err_yq,'b--',label='Fit $\pm 2\sigma$')
    plt.plot(xq,yq-2*err_yq,'b--')
    plt.legend(loc = 'right')
    myStr = "$\chi^2$ = {:.3f}".format(chisqr)
    plt.annotate(myStr, xy = (1.4,2.55))
    myStr = "Pendenza = {:.3f} $\pm$ {:.3f}".format(a[1],err_a[1])
    plt.annotate(myStr, xy = (1.4,2.4))

    plt.show()
    
#analyze()
#analyze("0",1)
#analyze("0",0,1)
