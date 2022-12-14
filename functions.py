# load similarities
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.dates as mdates
from matplotlib.pyplot import figure
from scipy.stats import invgamma
from datetime import datetime
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.ticker as ticker


# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 0, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()



def draw_plot(resMH_sumgamma, resMH_sumr, itertot, burnin, simu_01, 
              breakpoints, simu_01_date, simu_01_mean, simu_01_sd, 
              threshold_bp, threshold_func, reconstructiontot, path, 
              title="", showDate=False, save_fig=True):
        
    idx_ = pd.to_datetime(simu_01_date, format='___')
    
    fig_a = plt.figure(constrained_layout=True,figsize=(10, 10), dpi=250)
    gs = fig_a.add_gridspec(2, 2)

    # para ax1 es prob. posterior para breakpoints
    nPlot1 = len(resMH_sumgamma)
    
    # para ax2 es prob. functions
    nPlot2 = len(resMH_sumr)
    
    if showDate == True:

        ### ax1
        f_ax1 = fig_a.add_subplot(gs[0,0])
        # f_ax1.xaxis.set_minor_locator(mdates.MonthLocator())
        f_ax1.xaxis.set_minor_locator(ticker.MultipleLocator(115))
        f_ax1.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
        f_ax1.xaxis.set_major_locator(mdates.YearLocator())
        #f_ax1.xaxis.set_major_locator(ticker.AutoLocator())
        f_ax1.xaxis.set_major_formatter(mdates.DateFormatter('\n%Y'))
        f_ax1.tick_params(axis='x', which='both', labelsize=10)
        f_ax1.tick_params(axis='y', labelsize=10)

        f_ax1.xaxis.set_tick_params(labelsize='x-small')
        f_ax1.plot(idx_, resMH_sumgamma/(itertot-burnin), 
                   marker='o', markeredgecolor='black', markersize=3, 
                   color='white', linestyle='none')
        
        f_ax1.axhline(y=threshold_bp, color='red', linewidth=1)
        f_ax1.set_title('SelecciÃ³n puntos de cambio', size=14) # Breakpoints selection
        f_ax1.set_ylabel('Probabilidad posteriori', fontsize='large') # Posterior probabilities

        ### ax2
        f_ax2 = fig_a.add_subplot(gs[0,1])
        f_ax2.tick_params(axis='x', which='both', labelsize=10)
        f_ax2.tick_params(axis='y', labelsize=10)

        f_ax2.plot(np.arange(nPlot2), resMH_sumr/(itertot-burnin), 
                   marker='o', markersize=3, markeredgecolor='black', 
                   color='white', linestyle='none')

        f_ax2.axhline(y=threshold_func, color='red', linewidth=1)
        f_ax2.set_ylabel('Probabilidad posteriori', fontsize='large') # Posterior probabilities
        f_ax2.set_title('SelecciÃ³n de funciones', size=14) # Functions selection


        ### ax3
        f_ax3 = fig_a.add_subplot(gs[1, :])
        # f_ax3.xaxis.set_minor_locator(mdates.MonthLocator())
        f_ax3.xaxis.set_minor_locator(ticker.MultipleLocator(115))
        f_ax3.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
        f_ax3.xaxis.set_major_locator(mdates.YearLocator())
        f_ax3.xaxis.set_major_formatter(mdates.DateFormatter('\n%Y'))

        simu_01_ = simu_01 # * simu_01_sd + (simu_01_mean)
        reconstruction_ = reconstructiontot #.reshape(-1,1)
        reconstruction_ = reconstruction_ * simu_01_sd + simu_01_mean

        f_ax3.plot(idx_, simu_01_, color='black', ls='solid', lw=2, 
                   label="Serie {}".format(title)) # True serie
        f_ax3.plot(idx_, reconstruction_, color='#2D7ADD', lw=2, 
                   ls=(0, (5, 5)), label='Serie estimada') # Reconstruction    
    
        for i in breakpoints[1:]:
            f_ax3.axvline(x=idx_[i], color='red', linewidth=1)
    
    
        # where some data has already been plotted to ax
        handles, labels = f_ax3.get_legend_handles_labels()
        # Break points
        blue_line = mlines.Line2D([],[], color='red', 
                                  label='Puntos de cambio')

        # handles is a list, so append manual patch
        handles.append(blue_line) 

        f_ax3.legend(handles=handles)
        f_ax3.set_title('Serie {} y modelo estimado'.format(title), size=14) # True serie and reconstruction
        f_ax3.set_ylabel(title, fontsize='large') # Chilean Pesos
        f_ax3.set_xlabel('date', fontsize='large') # day
    
    else:
        ### ax1
        f_ax1 = fig_a.add_subplot(gs[0,0])
        #f_ax1.xaxis.set_tick_params(labelsize='x-small')
        f_ax1.plot(np.arange(nPlot1), resMH_sumgamma/(itertot-burnin), 
                   marker='o', markeredgecolor='black', markersize=3, 
                   color='white', linestyle='none')

        f_ax1.axhline(y=threshold_bp, color='red', linewidth=1)
        f_ax1.set_title('SelecciÃ³n puntos de cambio') # Breakpoints selection
        f_ax1.set_ylabel('Probabilidad posteriori') # Posterior probabilities

        ### ax2
        f_ax2 = fig_a.add_subplot(gs[0,1])
    
        f_ax2.plot(np.arange(nPlot2), resMH_sumr/(itertot-burnin), 
                   marker='o', markersize=3, markeredgecolor='black', 
                   color='white', linestyle='none')
        
        f_ax2.axhline(y=threshold_func, color='red', linewidth=1)
        f_ax2.set_ylabel('Probabilidad posteriori') # Posterior probabilities
        f_ax2.set_title('SelecciÃ³n de funciones') # Functions selection

        ### ax3
        f_ax3 = fig_a.add_subplot(gs[1, :])
        simu_01_ = simu_01 # * simu_01_sd + (simu_01_mean)
        reconstruction_ = reconstructiontot #.reshape(-1,1)
        reconstruction_ = reconstruction_ * simu_01_sd + simu_01_mean

        f_ax3.plot(np.arange(nPlot1), simu_01_, color='black', ls='solid', 
                   lw=2, label='Serie {}'.format(title)) # True serie
        f_ax3.plot(np.arange(nPlot1), reconstruction_, color='#2D7AC0', lw=2, 
                   ls=(0, (5, 5)), label='Serie estimada') # Reconstruction    
    
        for i in breakpoints[1:]:
            f_ax3.axvline(x=i, color='red', linewidth=1)
    
        # where some data has already been plotted to ax
        handles, labels = f_ax3.get_legend_handles_labels()
        #Break points
        blue_line = mlines.Line2D([],[], color='red', label='Puntos de cambio') 

        # handles is a list, so append manual patch
        handles.append(blue_line) 

        f_ax3.legend(handles=handles)
        f_ax3.set_title('Serie {}'.format(title)) # True serie and reconstruction
        f_ax3.set_ylabel('Niveles') # Chilean Pesos
        f_ax3.set_xlabel('DÃ­a') # day

    if save_fig == True:
        fig_a.savefig(path+'/resMH_{}.png'.format(title), dpi=200)
   
##############################################################################
# estimation of the breakpoints and of the functions composing
# the functional part (Metropolis-Hastings algorithm)
##############################################################################

# INPUTS
# serie: the observations
# nbiter: number of MH iterations total (burn-in + post burn-in)
# nburn: number of burn-in iterations
# lec1: chosen value for parameter c1
# lec2: chosen value for parameter c2
# Fmatrix: matrix which gives the values of the functions in the dictionary, 
# at each position.
# nbSegInit: initial number of segments for the segmentation part

# nbToChangegamma: number of gamma components proposed to be changed at each 
# iteration (when we propose to modify gamma), that is the number of inclusion 
# or deletion of breakpoints

# nbFuncInit: initial number of functions from the dictionary for the 
# functional part
# nbToChanger: number of r components proposed to be changed at each 
# iteration (when we propose to modify r), that is the number of inclusion or 
# deletion of functions from the dictionary
# Pi: vector with prior probabilities for the breakpoints: for position l 
# it is the prior probability that we observe a difference between the 
# segmentation part at time l and the segmentation part at time (l-1). 
# By convention Pi[1] = 1.
# eta: vector with prior probabilities for the functions in the dictionary: 
# eta[j] gives the prior proba that the function j from the dictionary will 
# be included in the functional part. By convention eta[1] = 1.
# printiter: if TRUE, the number of the actual iteration is plotted

# OUTPUTS
# sumgamma: for the breakpoints, for component l: number of iterations during 
# which a difference between the segmentation part at time l and the 
# segmentation part at time (l-1) was non nul (during post-burn-in).
# sumr: for functions from the dictionary, for component l: number of 
# iterations during which the function j from the dictionary was included in 
# the functional part (during post-burn-in).
# nbactugamma: number of iterations during which a MH proposal for the 
# breakpoints has been accepted (during which gamma has been modified) 
# (among iterations during burn-in and post-burn-in).
# nbactur: number of iterations during which a MH proposal for the functions 
# from the dictionary has been accepted (during which r has been modified) 
# (among iterations during burn-in and post-burn-in).
# gammamatrix: matrix to store all the post-burn-in simulations for gamma: 
# one line corresponds to one simulation (iteration)
# rmatrix: matrix to store all the post-burn-in simulations for gamma r: one 
# line corresponds to one simulation (iteration)

# serie tiene que ser un arreglo de nx1
def segmentation_bias_MH(serie, nbiter, nburn, lec1, lec2, Fmatrix, nbSegInit,
                         nbToChangegamma, nbFuncInit, nbToChanger, Pi, eta,
                         printiter=True):
    
    # de matrix to vector
    y = serie.reshape(-1)
    n = len(y)
    X = np.tri(n, n, 0, dtype=int)
    J = Fmatrix.shape[1]

    # to store results    
    gammamatrix = np.empty((nbiter-nburn, n))
    gammamatrix.fill(np.nan)
    rmatrix = np.empty((nbiter-nburn,J))
    rmatrix.fill(np.nan)
    
    sumgamma = np.zeros(n, int) 
    nbactugamma = 0
    sumr = np.zeros(J, int)
    nbactur = 0
        
    # initialization gamma 
    indgamma10 = np.random.choice(np.arange(1, n), size=nbSegInit-1, 
                                  replace=False)    
    gamma0 = np.zeros(n, int)
    gamma0[0] = 1
    
    # asigna desde 0 a (nbSegInit-1)-1
    for i in np.arange(nbSegInit-1):
        # printProgressBar(i, (nbSegInit-1), prefix = 'Progress:', suffix = 'Complete', length = 50)
        gamma0[indgamma10[i]] = 1

    indgamma1 = np.concatenate((np.array([0]), indgamma10))
    gamma = gamma0.copy()
    nbSeg = nbSegInit
    Xgamma = X[:,indgamma1] # fancy indexing
    invUgamma = np.diag(np.ones(n)) - lec1/(1+lec1) * (
        Xgamma @ np.linalg.inv(Xgamma.T @ Xgamma) @ Xgamma.T)
    
        
    # initialization r
    indr10 = np.random.choice(np.arange(1, J), size=nbFuncInit-1, 
                                  replace=False)
    
    r0 = np.zeros(J, int)
    r0[0] = 1
    for i in np.arange(nbFuncInit-1):
        # printProgressBar(i, (nbFuncInit-1), prefix = 'Progress:', suffix = 'Complete', length = 50)
        r0[indr10[i]] = 1

    indr1 = np.concatenate((np.array([0]), indr10))
    r = r0
    nbFunc = nbFuncInit
    Fmatrixr = Fmatrix[:,indr1]
    
    temp1 = np.linalg.det(
        np.linalg.inv(
            Fmatrixr.T @ (invUgamma + np.identity(n)/lec2) @ Fmatrixr
            )
        )
    
    temp2 = y.T @ (
        invUgamma - invUgamma @ Fmatrixr @ np.linalg.inv(
                Fmatrixr.T @  (invUgamma + np.identity(n)/lec2) @ Fmatrixr
            ) @ Fmatrixr.T @ invUgamma) @ y

    temp3 = np.linalg.det(np.dot(Fmatrixr.T, Fmatrixr))
    

    # iterations MH
    for iter in np.arange(nbiter):
        printProgressBar(iter, len(np.arange(nbiter)), prefix = 'Progress:', suffix = 'Complete MH', length = 30)

        if printiter == True:
            print("iter " + str(iter+1))
        # uniform choice between several movements
        choix = np.random.choice([1,2], size=1)
        
        # movement proposing to change only gamma
        if choix == 1:
            gammaprop = gamma.copy()
            indgamma1prop = indgamma1.copy()
            nbSegprop = nbSeg
            indToChange = np.random.choice(np.arange(1,n), 
                                           size=nbToChangegamma,
                                           replace=False)
            
            for i in np.arange(nbToChangegamma):
                if gamma[indToChange[i]] == 0:
                    gammaprop[indToChange[i]] = 1
                    indgamma1prop = np.concatenate([indgamma1prop, 
                                                   [indToChange[i]]])                   
                    nbSegprop = nbSegprop + 1
                else:
                    gammaprop[indToChange[i]] = 0                    
                    indremove = np.nonzero(indgamma1prop==indToChange[i])[0]
                    indgamma1prop = np.delete(indgamma1prop, indremove[0])          
                    nbSegprop = nbSegprop - 1
            
            Xgammaprop = X[:,indgamma1prop]
 
            invUgammaprop = np.identity(n) - (lec1/(1+lec1) * (
                Xgammaprop @ np.linalg.inv(Xgammaprop.T @ Xgammaprop) @ 
                Xgammaprop.T))

            # new nnarria
            tmpsolve_ = np.linalg.inv(
                Fmatrixr.T @ (invUgammaprop + np.identity(n)/lec2) @ Fmatrixr)
           
            temp1prop = np.linalg.det(tmpsolve_)
            temp2prop = y.T @ (invUgammaprop - invUgammaprop @ Fmatrixr @ (
                    tmpsolve_
                ) @ Fmatrixr.T @ invUgammaprop
            ) @ y
    
            A = np.power(1+lec1, (nbSeg-nbSegprop)/2) * np.prod(
                    np.power(Pi[1:]/(1-Pi[1:]),(gammaprop-gamma)[1:])
                ) * np.power(temp1prop/temp1, 1/2) * (
                np.float_power(temp2/temp2prop, n/2))
            
            probaccept1 = min(1, A)
            seuil = np.random.uniform(0,1,1)[0] # runif(1)
            
            if seuil < probaccept1:
                gamma = gammaprop.copy()                
                indgamma1 = indgamma1prop.copy()
                nbSeg = nbSegprop
                Xgamma = Xgammaprop.copy()
                invUgamma = invUgammaprop.copy()
                temp1 = temp1prop
                temp2 = temp2prop.copy()
                nbactugamma = nbactugamma + 1
            
            
        # movement proposing to change only r
        if choix == 2:
            rprop = r.copy()
            indr1prop = indr1.copy()
            nbFuncprop = nbFunc
            indToChange = np.random.choice(np.arange(1,J), 
                                           size=nbToChanger,
                                           replace=False)

            for i in np.arange(nbToChanger):
                if r[indToChange[i]] == 0:
                    rprop[indToChange[i]] = 1
                    indr1prop = np.concatenate([indr1prop, [indToChange[i]]])
                    nbFuncprop = nbFuncprop + 1
                else:
                    rprop[indToChange[i]] = 0
                    indremove = np.nonzero(indr1prop == indToChange[i])[0]
                    indr1prop = np.delete(indr1prop, indremove[0])
                    nbFuncprop = nbFuncprop - 1

            Fmatrixrprop = Fmatrix[:, indr1prop]
            
            # new nnarria
            tmpsolve_ = np.linalg.inv(Fmatrixrprop.T @ (
                    invUgamma + np.identity(n)/lec2) @ Fmatrixrprop)
            
            temp1prop = np.linalg.det(tmpsolve_)
            temp2prop = y.T @ (invUgamma - invUgamma @ Fmatrixrprop @  
                               tmpsolve_ @ Fmatrixrprop.T @ invUgamma) @ y
            temp3prop = np.linalg.det(Fmatrixrprop.T @ Fmatrixrprop)  ###  warning in algorithn ###


            A = np.power(lec2, (nbFunc-nbFuncprop)/2) * np.prod(
                    np.power(eta[1:]/(1-eta[1:]),(rprop-r)[1:])
                ) * np.power(temp1prop/temp1, 1/2) * (
                    np.float_power(temp2/temp2prop, n/2)) * ( 
                        np.power(temp3prop/temp3, 1/2))
                         
                        
            probaccept1 = min(1, A)
            seuil = np.random.uniform(0,1,1)[0] # runif(1)   
            
            if seuil < probaccept1:
                r = rprop.copy()
                indr1 = indr1prop.copy()
                nbFunc = nbFuncprop
                Fmatrixr = Fmatrixrprop.copy()
                temp1 = temp1prop
                temp2 = temp2prop.copy()
                temp3 = temp3prop
                nbactur = nbactur + 1
        
        # store results when we are in post-burn-in
        if iter >= nburn:
            sumgamma = sumgamma + gamma
            sumr = sumr + r
            gammamatrix[(iter-nburn),:] = gamma.copy()
            rmatrix[(iter-nburn),:] = r.copy()
    
    full = len(np.arange(nbiter))
    printProgressBar(full, full, prefix = 'Progress:', suffix = 'Complete', length = 30)
    # devuelve un diccionario
    return dict(sumgamma=sumgamma, sumr=sumr, nbactugamma=nbactugamma, 
            nbactur=nbactur, gammamatrix=gammamatrix, rmatrix=rmatrix)

##############################################################################
# After the segmentation: 
# Estimation of betagamma, 
# lambdar and sigma2 (Gibbs sampler)   
##############################################################################

# INPUTS
# serie: the observations
# nbiter: total number of iterations for the Gibbs sampler (burn-in + 
# post burn-in)
# nburn: number of burn-in iterations
# lec1: chosen value for the c1 parameter
# lec2: chosen value for the c2 parameter
# Fmatrix: matrix which gives the values of the functions in the dictionary, 
# at each position.
# gammahat: estimated gamma vector (using the preceding MH algo) 
# rhat: estimated r vector (using the preceding MH algo) 
# priorminsigma2: prior minimum value for sigma2
# priormaxsigma2: prior maximum value for sigma2
# printiter: if TRUE, the number of the actual iteration is plotted

# OUTPUTS
# resbetagammahat: matrix to store all the post-burn-in simulations for 
# betagamma: one line corresponds to one simulation (iteration)
# reslambdarhat: matrix to store all the post-burn-in simulations for 
# lambdar: one line corresponds to one simulation (iteration)
# ressigma2hat: matrix to store all the post-burn-in simulations for sigma2: 
# one line corresponds to one simulation (iteration)
# estbetagamma: estimation of betagamma (mean of all the simulated betagamma)
# estlambdar: estimation of lambdar (mean of all the simulated lambdar)
# estsigma2: estimation of sigma2 (mean of all the simulated sigma2)

def estimation_moy_biais (serie, nbiter, nburn, lec1, lec2, Fmatrix, gammahat,
                          rhat, priorminsigma2, priormaxsigma2, 
                          printiter=True):    
    y = serie[0:,0]
    n = len(y)
    X = np.tri(n, n, 0, dtype=int)    
   
    dgammahat = np.sum(gammahat)
    drhat = np.sum(rhat)
   
    if drhat != 0:
        
        ### 
        directions_x3 = np.array([[ 0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0],
                       [2, 0, 1], [2, 1, 0]])
        matrix_dir_x3 = np.random.choice(np.arange(directions_x3.shape[0]), 
                                                   size=nbiter)
        ###
        # useful quantities for calculations
        Xgamma = X[:,np.nonzero(gammahat==1)[0]]
        Fmatrixr = Fmatrix[:,np.nonzero(rhat==1)[0]]
        temp1 = np.linalg.inv((1+lec1)/lec1 * Xgamma.T @ Xgamma)
        temp2 = np.linalg.inv((1+lec2)/lec2 * Fmatrixr.T @ Fmatrixr)
        alpha = (n+dgammahat+drhat)/2
        temp3 = Xgamma.T @ Xgamma/lec1
        temp4 = Fmatrixr.T @ Fmatrixr/lec2
       
        # to store results     
        resbetagammahat = np.zeros(dgammahat*(nbiter-nburn)).reshape(
            dgammahat, -1)       
        reslambdarhat = np.zeros(int(drhat*(nbiter-nburn))).reshape(
            int(drhat), -1)
        ressigma2hat = np.zeros(nbiter-nburn)

       
        # initializations
        sigma2hat = np.exp(np.random.uniform(np.log(priorminsigma2), 
                                          np.log(priormaxsigma2), 1))[0]
         
        mean_ = np.zeros(dgammahat)
        cov_ = lec1*sigma2hat*np.linalg.inv(Xgamma.T @ Xgamma)
        betagammahat = np.random.multivariate_normal(mean_, cov_, 1)[0]
      
        mean_ = np.zeros(int(drhat))
        cov_ = lec2*sigma2hat*np.linalg.inv(Fmatrixr.T @ Fmatrixr)
        lambdarhat = np.random.multivariate_normal(mean_, cov_, 1)[0]
        
        # iterations
        for iter in np.arange(nbiter):
            # printProgressBar(iter, nbiter, prefix = 'Progress:', suffix = 'Complete', length = 30)
            if printiter == True:
                print('iter ' + str(iter))
            
            ordre = directions_x3[matrix_dir_x3[iter]]
            
            for j in np.arange(3):
                if ordre[j] == 0:
                    mean_ = temp1 @ Xgamma.T @ (y-Fmatrixr @ lambdarhat)
                    cov_ = sigma2hat*temp1

                    # if exist any nan or inf or -inf
                    if not(np.any(np.isnan(mean_))):
                        betagammahat = np.random.multivariate_normal(mean_, cov_, 1)[0]
                    else:
                        np.nan
                
                if ordre[j] == 1:
                    mean_ = temp2 @ Fmatrixr.T @ (y-Xgamma @ betagammahat)
                    cov_ = sigma2hat*temp2

                    # if exist any nan or inf or -inf
                    if not(np.any(np.isnan(mean_))):
                        lambdarhat = np.random.multivariate_normal(mean_, cov_, 1)[0]
                    else:
                        np.nan
                    
                if ordre[j] == 2:
                    scale_ = 1/2 * (
                        (y-Xgamma @ betagammahat-Fmatrixr @ lambdarhat).T @ (
                            y-Xgamma @ betagammahat-Fmatrixr @ lambdarhat
                            ) + betagammahat.T @ temp3 @ betagammahat +
                            lambdarhat.T @ temp4 @ lambdarhat)

                    # if exist any nan or inf or -inf
                    if not(np.any(np.isnan(scale_))):
                        sigma2hat = invgamma.rvs(alpha, scale=scale_, size=1)
                    else:
                        np.nan
            
            if iter >= nburn:
                resbetagammahat[:, iter-nburn] = betagammahat
                reslambdarhat[:, iter-nburn] = lambdarhat
                ressigma2hat[iter-nburn] = sigma2hat
        
        estbetagamma = np.sum(resbetagammahat, axis=1)/(nbiter-nburn)
        estlambdar = np.sum(reslambdarhat, axis=1)/(nbiter-nburn)
        estsigma2 = np.sum(ressigma2hat)/(nbiter-nburn)
            
    if drhat == 0:
        
        ###
        directions_x2 = np.array([[ 0, 1],[1, 0]])
        matrix_dir_x2 = np.random.choice(np.arange(directions_x2.shape[0]), 
                                                   size=nbiter)
        ###
        
        # useful quantities for calculations
        Xgamma = X[:,np.nonzero(gammahat==1)[0]]
        temp1 = np.linalg.inv((1+lec1)/lec1 * Xgamma.T @ Xgamma)
        alpha = (n+dgammahat)/2
        temp3 = Xgamma.T @ Xgamma/lec1

        # to store results
        resbetagammahat = np.zeros(dgammahat*(nbiter-nburn)).reshape(
            dgammahat, -1)
        reslambdarhat = None
        ressigma2hat = np.zeros(nbiter-nburn)
        
        # initializations
        sigma2hat = np.exp(np.random.uniform(np.log(priorminsigma2), 
                                          np.log(priormaxsigma2), 1))[0]
        mean_ = np.zeros(dgammahat)
        cov_ = lec1*sigma2hat*np.linalg.inv(Xgamma.T @ Xgamma)
        betagammahat = np.random.multivariate_normal(mean_, cov_, 1)[0]
     
        # iterations
        for iter in np.arange(nbiter):
            # printProgressBar(iter, nbiter, prefix = 'Progress:', suffix = 'Complete', length = 30)
            if printiter == True:
                print("iter " + str(iter))
            
            ordre = directions_x2[matrix_dir_x2[iter]]
            
            for j in np.arange(2):
                if ordre[j] == 0:
                    mean_ = temp1 @ Xgamma.T @ y
                    cov_ = sigma2hat*temp1
                    
                    # if exist any nan or inf or -inf
                    if not(np.any(np.isnan(mean_))):
                        betagammahat = np.random.multivariate_normal(mean_, cov_, 1)[0]
                    else:
                        np.nan

                if ordre[j] == 1:                    
                    scale_ = 1/2 * ((y-Xgamma @ betagammahat).T @ (
                            y-Xgamma @ betagammahat) + betagammahat.T @ (
                                temp3 @ betagammahat))

                    # if exist any nan or inf or -inf
                    if not(np.any(np.isnan(scale_))):
                        sigma2hat = invgamma.rvs(alpha, scale=scale_, size=1)
                    else:
                        np.nan
                    
            
            if iter >= nburn:
                resbetagammahat[:, iter-nburn] = betagammahat
                ressigma2hat[iter-nburn] = sigma2hat
        
        estbetagamma = np.sum(resbetagammahat, axis=1)/(nbiter-nburn)
        estlambdar = None
        estsigma2 = np.sum(ressigma2hat)/(nbiter-nburn)     

    return list([resbetagammahat, reslambdarhat, ressigma2hat, 
                           estbetagamma, estlambdar, estsigma2])

##############################################################################
# AThe two steps of M-H and Gibbs sampler are integrated for the 
# estimation of points of change with functional effect 
##############################################################################

# INPUTS
# Fmatrix: Matrix with the dictionary data
# data_serie:data series to which you want to detect change points
# itertot: Total number of iterations to be used in both M-H and Gibbs Sampler
# burnin: Number of iterations that will not be considered
# lec1: chosen value for the c1 parameter
# lec2: chosen value for the c2 parameter
# nbseginit: initial number of segments for the segmentation part
# nfuncinit: initial number of functions from the dictionary for the 
# functional part
# nbtochangegamma: number of gamma components proposed to be changed at each 
# iteration (when we propose to modify gamma), that is the number of inclusion 
# or deletion of breakpoints
# nbtochanger: number of r components proposed to be changed at each 
# iteration (when we propose to modify r), that is the number of inclusion or 
# deletion of functions from the dictionary
# Pi: vector with prior probabilities for the breakpoints: for position l 
# it is the prior probability that we observe a difference between the 
# segmentation part at time l and the segmentation part at time (l-1). 
# By convention Pi[1] = 1.
# eta: vector with prior probabilities for the functions in the dictionary: 
# eta[j] gives the prior proba that the function j from the dictionary will 
# be included in the functional part. By convention eta[1] = 1.
# threshold_bp: probability threshold from which the change points 
# are selected
# threshold_fnc: probability threshold from which the functions are selected
# printiter: if TRUE, the number of the actual iteration is plotted

    
# OUTPUTS
# list with;
# [0]: list with
# sumgamma: for the breakpoints, for component l: number of iterations 
# during which a difference between the segmentation part at time l and the 
# segmentation part at time (l-1) was non nul (during post-burn-in).
# sumr: for functions from the dictionary, for component l: number of 
# iterations during which the function j from the dictionary was included 
# in the functional part (during post-burn-in).
# nbactugamma: number of iterations during which a MH proposal for the 
# breakpoints has been accepted (during which gamma has been modified) 
# (among iterations during burn-in and post-burn-in).
# nbactur: number of iterations during which a MH proposal for the functions 
# from the dictionary has been accepted (during which r has been modified) 
# (among iterations during burn-in and post-burn-in).
# gammamatrix: matrix to store all the post-burn-in simulations for gamma: 
# one line corresponds to one simulation (iteration).
# rmatrix: matrix to store all the post-burn-in simulations for gamma r: one 
# line corresponds to one simulation (iteration)
# [1] resbetagammahat: matrix to store all the post-burn-in simulations for.
# [2] reslambdarhat: matrix to store all the post-burn-in simulations for. 
# [3] ressigma2hat: matrix to store all the post-burn-in simulations for.
# [4] estbetagamma: estimation of betagamma (mean of all the simulated 
# betagamma)
# [5] estlambdar: estimation of lambdar (mean of all the simulated lambdar)
# [6] estsigma2: estimation of sigma2 (mean of all the simulated sigma2)
# [7] reconstructiontot
# [8] basefunctions

def dbp_with_function_effect (
        Fmatrix, data_serie, itertot, burnin, 
        lec1, lec2, nbseginit, nbfuncinit, nbtochangegamma, nbtochanger,
        Pi, eta, threshold_bp, threshold_fnc, printiter=False):
    
    n = len(data_serie)
    printProgressBar(0, n, prefix = 'Progress:', suffix = '', length = 30)

    # result Metropolis Hastings
    resMH = segmentation_bias_MH (
        data_serie, itertot, burnin, lec1, 
        lec2, Fmatrix, nbseginit, 
        nbtochangegamma, nbfuncinit, 
        nbtochanger, Pi, eta, 
        printiter=printiter)
    
    print("resMH calculado")
    
    # Selection of the points of change
    breakpoints = np.asarray(
    np.nonzero(resMH['sumgamma']/(itertot-burnin) > threshold_bp))[0]
    print(f"breakpoints: {str(breakpoints)}")

    gammahat = np.zeros(n, int)
    for i in (breakpoints):
        gammahat[i] = 1
    printProgressBar(10, n, prefix = 'BPoints:', suffix = 'Complete', length = 30)

    # Selection of the functions
    basefunctions = np.asarray(
        np.nonzero(resMH['sumr']/(itertot-burnin) > threshold_fnc))[0]
    print(f"basefunctions: {str(basefunctions)}")
    
    rhat = np.zeros(Fmatrix.shape[1])
    for i in basefunctions[1:]:
        rhat[i] = 1
    printProgressBar(20, n, prefix = 'BPoints:', suffix = 'Complete', length = 30)
    # Estimation of betagamma, lambdar and sigma2
    # To calculate the hope of the series, the segmentation and 
    # functional part is estimated and then added.        
    priorminsigma2 = 0.001
    priormaxsigma2 = 5
        
    estim = estimation_moy_biais(
        data_serie, itertot, burnin, lec1, lec2, 
        Fmatrix, gammahat, rhat, priorminsigma2, 
        priormaxsigma2, printiter=printiter)
    
    muest = np.zeros(breakpoints.shape[0])
    muest[0] = estim[3][0]
    reconstructionmu = np.zeros(Fmatrix.shape[0])
    breakpoints_extend = np.concatenate([breakpoints, [Fmatrix.shape[0]]])
        
    if breakpoints.shape[0] > 1:
        muest = np.cumsum(estim[3])
                
        for i in np.arange(len(breakpoints_extend)-1):
            printProgressBar(i, len(breakpoints_extend), prefix = 'Progress:', suffix = 'MU Complete', length = 30)
            reconstructionmu[breakpoints_extend[i]:breakpoints_extend[i+1]] = muest[i]
        
    reconstructionf = np.zeros(Fmatrix.shape[0])
    if basefunctions.shape[0] > 1:
        for i in np.arange(1, basefunctions.shape[0]):
            printProgressBar(i, len(basefunctions), prefix = 'Progress:', suffix = 'Fun Complete', length = 30)
            reconstructionf = reconstructionf + estim[4][i-1] * (Fmatrix[:,basefunctions[i]])
        
    reconstructiontot = reconstructionmu + reconstructionf
    # print(f"Reconstruction ok")

    
    printProgressBar(n, n, prefix = 'Progress:', suffix = 'DBP Complete', length = 30)

    return list([resMH, estim[0], estim[1], estim[2], estim[3], estim[4], 
                 estim[5], reconstructiontot, basefunctions, Fmatrix])

##############################################################################
#
##############################################################################
# 

def dbp_function_with_resMH(resMH,
        Fmatrix, data_serie, itertot, burnin, 
        lec1, lec2, threshold_bp, threshold_fnc, printiter=False):
    
    n = len(data_serie)    
    printProgressBar(0, n, prefix = 'dpb_function:', suffix = '', length = 30)
    
    # Selection of the points of change
    breakpoints = np.asarray(
    np.nonzero(resMH['sumgamma']/(itertot-burnin) > threshold_bp))[0]
    print(f"breakpoints: {str(breakpoints)}")

    gammahat = np.zeros(n, int)
    for i in (breakpoints):
        gammahat[i] = 1
    printProgressBar(10, n, prefix = 'gammahat:', suffix = 'Complete', length = 30)

    # Selection of the functions
    basefunctions = np.asarray(
        np.nonzero(resMH['sumr']/(itertot-burnin) > threshold_fnc))[0]
    print(f"basefunctions: {str(basefunctions)}")
    
    rhat = np.zeros(Fmatrix.shape[1])
    for i in basefunctions[1:]:
        rhat[i] = 1
    printProgressBar(20, n, prefix = 'rhat:', suffix = 'Complete', length = 30)
    # Estimation of betagamma, lambdar and sigma2
    # To calculate the hope of the series, the segmentation and 
    # functional part is estimated and then added.        
    priorminsigma2 = 0.001
    priormaxsigma2 = 5
        
    estim = estimation_moy_biais(
        data_serie, itertot, burnin, lec1, lec2, 
        Fmatrix, gammahat, rhat, priorminsigma2, 
        priormaxsigma2, printiter=printiter)

    printProgressBar(30, n, prefix = 'moy_biasis:', suffix = 'Complete', length = 30)
    
    muest = np.zeros(breakpoints.shape[0])
    muest[0] = estim[3][0]
    reconstructionmu = np.zeros(Fmatrix.shape[0])
    breakpoints_extend = np.concatenate([breakpoints, [Fmatrix.shape[0]]])
        
    if breakpoints.shape[0] > 1:
        muest = np.cumsum(estim[3])
                
        for i in np.arange(len(breakpoints_extend)-1):
            printProgressBar(i, len(breakpoints_extend), prefix = 'breakpoints:', suffix = 'MU Complete', length = 30)
            reconstructionmu[breakpoints_extend[i]:breakpoints_extend[i+1]] = muest[i]
        
    reconstructionf = np.zeros(Fmatrix.shape[0])
    if basefunctions.shape[0] > 1:
        for i in np.arange(1, basefunctions.shape[0]):
            printProgressBar(i, len(basefunctions), prefix = 'functions:', suffix = 'Fun Complete', length = 30)
            reconstructionf = reconstructionf + estim[4][i-1] * (Fmatrix[:,basefunctions[i]])
        
    reconstructiontot = reconstructionmu + reconstructionf
    
    printProgressBar(n, n, prefix = 'Estims: ', suffix = 'end', length = 30)

    return list([resMH, estim[0], estim[1], estim[2], estim[3], estim[4], 
                 estim[5], reconstructiontot, basefunctions, breakpoints_extend, Fmatrix])


##############################################################################
#
##############################################################################
# 

def dbp_extern_reconstruction(Fmatrix=None, estim=[], breakpoints=[], basefunctions=[]):
    """
    """
    reconstructionmu = np.zeros(Fmatrix.shape[0])
    breakpoints_extend = np.concatenate([breakpoints, [Fmatrix.shape[0]]])
    
    print(">>>>>>>>> ", breakpoints_extend.shape)
        
    if breakpoints.shape[0] > 1:
        muest = np.cumsum(estim[0])
        print(len(breakpoints_extend), len(muest))

        for i in np.arange(len(breakpoints_extend)-1):
            printProgressBar(i, len(breakpoints_extend), prefix = 'breakpoints:', suffix = 'MU Complete', length = 30)
            reconstructionmu[breakpoints_extend[i]:breakpoints_extend[i+1]] = muest[i]
        
    reconstructionf = np.zeros(Fmatrix.shape[0])
    if basefunctions.shape[0] > 1:
        for i in np.arange(1, basefunctions.shape[0]):
            printProgressBar(i, len(basefunctions), prefix = 'functions:', suffix = 'Fun Complete', length = 30)
            reconstructionf = reconstructionf + estim[1][i-1] * (Fmatrix[:,basefunctions[i]])
        
    reconstructiontot = reconstructionmu + reconstructionf
    printProgressBar(0, 30, prefix = 'Reconstruction: ', suffix = 'end', length = 30)
    return reconstructiontot


##############################################################################
#
##############################################################################
# 

def dbp_extern_reconstruction2(Fmatrix=None, estim=[], basefunctions=[]):
    """
    """
    reconstructionf = np.zeros(Fmatrix.shape[0])
    if basefunctions.shape[0] > 1:
        for i in np.arange(1, basefunctions.shape[0]):
            # print("i: ", i, estim[i-1])
            reconstructionf = reconstructionf + estim[i-1] * (Fmatrix[:,basefunctions[i]])
    
    return reconstructionf



##########################
### Simulation of series
##########################
# function to simulate M series, each one having their own breakpoints, but 
# all sharing the same number of segments, and the same functional part which 
# is specified in the article.

# INPUTS
# M: the number of series (if M>1, we are in case of multplie series)
# n: the number of points of the series
# K: the number of segments for each serie
# mu.choix: vector which gives the possible values for each segment of the 
# segmentation part of the series. The values of each segment are unifomly 
# sampled from this vector, and the value of a segment should not be equal to 
# the value of the preceding segment.
# varianceError: vector of variances for the error term
# pHaar: the positions of the Haar functions composing the functional part
# p1: A breakpoint should be positioned at a distance from the Haar functions 
# of at least p1
# p2: Each segment is at least of length p2

# OUTPUTS
# K: the number of segments for each serie
# muMat: matrix (M x n) which gives the value of the segmentation part of 
# each serie at each position
# tauMat: matrix (M x k.max) which gives the index of the breakpoints for each 
# series (K points for each series)
# errorMat: matrice (M x n) which gives the error term at every time point, 
# for each series
# biais: evaluation of the functional part at every time point
# Profile: data.frame which summarizes all the preceding informations in 
# columns: number (id) of the series, position in the grid of time points, 
# value of the segmentation part at this position, presence or not of a 
# breakpoint at this position, value of the functional part at this position, 
# error term at this position (one column for each specified variance of the 
# error term)

def simuSeries(M, n, K, muChoix, varianceError, pHaar, p1, p2):
    
    #### only 1 series with 4 segments, with a functional part, 100 points. 
    #We specify here only one variance for the error term only.
    #M = 1
    #n = 100
    #K = 4
    #muChoix = np.array([0, 1, 2, 3, 4, 5])
    #standard_deviation = np.array([0.1])
    #varianceError = np.power(standard_deviation, 2)
    #pHaar = np.array([10, 50, 60])
    #p1 = 3
    #p2 = 5
    #nbsimu = 1
    
    # Construction of the functional part f
    t = np.arange(n)
    A = 0.3
    Part1 = A*np.sin(2*np.pi*t/20)
    Part2 = np.zeros(n)
    t1 = pHaar[0]
    t2 = pHaar[1]
    t3 = pHaar[2]
    Part2[t1] = 1.5
    Part2[t2] = -2
    Part2[t3] = 3
    biais = Part1 + Part2
    
    # construction of the M series
    erreurs = list()
    muMat = np.empty((M, n), int)
    muMat[:] = -1
    tauMat = np.empty((M, K), int)

    # errors: a matrix (M x n) for each possible value of variance, 
    # giving the errors of each series at each position
    for i in np.arange(len(varianceError)):
        errorMat = np.empty((M, n))
        errorMat[:] = np.NaN
        for m in np.arange(M):
            errorMat[m,] = np.random.normal(0, np.sqrt(varianceError[i]), n)   
        erreurs.append(errorMat)
  
    for m in np.arange(M):
        # positions of breakpoints pour series m
        cond = 0
        while (cond == 0):
            # generar un numero aleatoriao uniforme entre
            # los n posibles puntos excluyendo el ultimo
            tauTmp = np.ceil((n-1)*np.random.uniform(0,1,K-1))
            tauTmp = np.sort(tauTmp)
            tauTmp = np.concatenate([tauTmp,[n]]).astype(int).copy()
            cond2 = 1
            
            for i in np.arange(1,K):
                cond2 = cond2 * ((tauTmp[i]-tauTmp[i-1])>=p2)
                
            for i in np.arange(K):
                for j in np.arange(len(pHaar)):
                    cond2 = cond2 * np.abs((tauTmp[i]-pHaar[j]))>=p1
            cond = cond2
        tauMat[m, np.arange(K)] = tauTmp
        
        # values of segments (segmentation part) for series m
        mutemp = np.random.choice(muChoix, size=1)
        muMat[m, np.arange(tauMat[m,0])] = np.repeat(mutemp, tauMat[m,0])
        muChoixTemp = muChoix
        
        if K > 1:
            for k in np.arange(1, K):
                toremove = np.nonzero(muChoixTemp == mutemp)
                muChoixTemp = np.delete(muChoixTemp, toremove).copy()
                mutemp = np.random.choice(muChoixTemp, size=1)

                muMat[m,np.arange(tauMat[m,k-1], tauMat[m,k])] = (
                      np.repeat(mutemp, (tauMat[m, k]-tauMat[m, k-1]))
                      )

    series = np.array([], int)
    for i in np.arange(M):
        series = np.concatenate([series, np.repeat(i, n)])
    
    position = np.repeat(np.arange(n), M)
    mu = (muMat.T).copy().reshape((M*n,1))
    errors = erreurs[0].T.copy().reshape((M*n, 1))
    nomcolonnes = "erreur1"
    
    if len(varianceError) > 1:
        for i in np.arange(1,len(varianceError)):
            errors = pd.concat(errors, erreurs[i].T.copy().reshape(M*n, 1))
            nomcolonnes = np.concatenate([[nomcolonnes], "erreur"+str(i)])
        
    tau = np.array([], int)
    for m in np.arange(M):
        tauPos = np.repeat(0, n)
        tauPos[tauMat[m,]-1] = 1
        tau = np.concatenate([tau,tauPos])

    # Create a zipped list of tuples from above lists
    zippedList =  list(zip(series, position, mu.reshape(n), tau, biais, 
                       errors.reshape(n)))

    Profile = pd.DataFrame(zippedList, columns = ['series' , 'position', 'mu', 
                                                  'tau', 'biais', nomcolonnes])
    
    # output
    return(list([[K],muMat,tauMat,errorMat,biais,Profile]))

###
# A recursive funtions to found 
# offset in serie dataset
# :param: position to found
# :param: serie dataset
#Â 
def translation(a_serie=[]):    
    arr = np.array(range(len(a_serie)))
    arr = arr[[bool(val) for val in  ~np.isnan(a_serie)]]
    return arr

###
# 
# :param: position to found
# :param: serie dataset
# 
# def draw_serie(serie, breakpoints, date, 
#                 serie_sd, serie_mean, result,
#                 save_fig=True, path=".", title="",
#                 style="ggplot"):
    # idx_ = pd.to_datetime(date, format='___')
    # fig = figure(figsize=(20, 5), dpi=500)
    # plt.style.use(style)

    # plt.plot(idx_, serie)

    # # plt.plot(idx_, serie, color='black', ls='solid', lw=2, 
    # #             label="Serie {}".format(title)) # True serie
    
    # plt.plot(idx_, result, color='#2D7AC0', lw=3, 
    #             ls=(0, (5, 1)), label='Serie estimada') # Reconstruction    

    # # for i in breakpoints[1:]:
    # #     plt.axvline(x=idx_[i], color='red', linewidth=2)

    # plt.title("{}".format(title))

    # if save_fig == True:
    #     plt.savefig(path+'/fig_{}_processed.png'.format(title))
    #     plt.close(fig)

def draw_serie(resMH_sumgamma, resMH_sumr, itertot, burnin, simu_01, 
              breakpoints, simu_01_date, simu_01_mean, simu_01_sd, 
              threshold_bp, threshold_func, reconstructiontot, path, 
              style="grayscale", reconstructionalt=None, title="", 
              showDate=False, save_fig=True):
        
    idx_ = pd.to_datetime(simu_01_date, format='___')
    
    fig_a = plt.figure(constrained_layout=True, figsize=(10, 10), dpi=250)    
    gs = fig_a.add_gridspec(2, 2)

    # para ax1 es prob. posterior para breakpoints
    nPlot1 = len(resMH_sumgamma)
    
    # para ax2 es prob. functions
    nPlot2 =  len(resMH_sumr)

   

    plt.style.use(style)
    
    if showDate == True:

        ### ax1
        f_ax1 = fig_a.add_subplot(gs[0,0])
        # f_ax1.xaxis.set_minor_locator(mdates.MonthLocator())
        f_ax1.xaxis.set_minor_locator(ticker.MultipleLocator(115))
        f_ax1.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
        f_ax1.xaxis.set_major_locator(mdates.YearLocator())
        #f_ax1.xaxis.set_major_locator(ticker.AutoLocator())
        f_ax1.xaxis.set_major_formatter(mdates.DateFormatter('\n%Y'))
        f_ax1.tick_params(axis='x', which='both', labelsize=10)
        f_ax1.tick_params(axis='y', labelsize=10)

        f_ax1.xaxis.set_tick_params(labelsize='x-small')
        f_ax1.plot(idx_, resMH_sumgamma/(itertot-burnin), 
                   marker='o', markeredgecolor='black', markersize=3, 
                   color='white', linestyle='none')
        
        f_ax1.axhline(y=threshold_bp, color='red', linewidth=1)
        f_ax1.set_title('Breakpoints selection', size=14) # Breakpoints selection
        f_ax1.set_ylabel('Posterior Probability', fontsize='large') # Posterior probabilities

        ### ax2
        f_ax2 = fig_a.add_subplot(gs[0,1])
        f_ax2.tick_params(axis='x', which='both', labelsize=10)
        f_ax2.tick_params(axis='y', labelsize=10)

        f_ax2.plot(np.arange(nPlot2), resMH_sumr/(itertot-burnin), 
                   marker='o', markersize=3, markeredgecolor='black', 
                   color='white', linestyle='none')

        f_ax2.axhline(y=threshold_func, color='red', linewidth=1)
        f_ax2.set_ylabel('Posterior Probability', fontsize='large') # Posterior probabilities
        f_ax2.set_title('Breakpoints selection', size=14) # Functions selection


        ### ax3
        f_ax3 = fig_a.add_subplot(gs[1, :])
        # f_ax3.xaxis.set_minor_locator(mdates.MonthLocator())
        f_ax3.xaxis.set_minor_locator(ticker.MultipleLocator(115))
        f_ax3.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
        f_ax3.xaxis.set_major_locator(mdates.YearLocator())
        f_ax3.xaxis.set_major_formatter(mdates.DateFormatter('\n%Y'))

        simu_01_ = simu_01 # * simu_01_sd + (simu_01_mean)
        reconstruction_ = reconstructiontot #.reshape(-1,1)
        reconstruction_ = reconstruction_ * simu_01_sd + simu_01_mean

        f_ax3.plot(idx_, simu_01_, ls='solid', lw=2, 
                   label="Serie {}".format(title)) # True serie
        f_ax3.plot(idx_, reconstruction_, color='#0f3cf1', lw=2, 
                   ls=(0, (5, 5)), label='Reconstruction') # Reconstruction    
        
        if not(reconstructionalt is None):
            f_ax3.plot(idx_, reconstructionalt, color='#AF601A', lw=2, 
                    ls=(0, (5, 5)), label='Serie Alt') # Reconstruction  
    
        for i in breakpoints[1:]:
            f_ax3.axvline(x=idx_[i], color='#D35400', linewidth=1.5)
            
        # where some data has already been plotted to ax
        handles, labels = f_ax3.get_legend_handles_labels()
        # Break points
        blue_line = mlines.Line2D([],[], color='#873600', 
                                  label='Breakpoints')

        # handles is a list, so append manual patch
        handles.append(blue_line) 

        f_ax3.legend(handles=handles)
        f_ax3.set_title('True serie {} and reconstruction'.format(title), size=14) # True serie and reconstruction
        f_ax3.set_ylabel(title, fontsize='large') # Chilean Pesos
        f_ax3.set_xlabel('date', fontsize='large') # day
    
    else:
        ### ax1
        f_ax1 = fig_a.add_subplot(gs[0,0])
        #f_ax1.xaxis.set_tick_params(labelsize='x-small')
        f_ax1.plot(np.arange(nPlot1), resMH_sumgamma/(itertot-burnin), 
                   marker='o', markeredgecolor='black', markersize=3, 
                   color='white', linestyle='none')

        f_ax1.axhline(y=threshold_bp, color='#F8C471', linewidth=1)
        f_ax1.set_title('Breakpoints selection') # Breakpoints selection
        f_ax1.set_ylabel('Posterior probabilities') # Posterior probabilities

        ### ax2
        f_ax2 = fig_a.add_subplot(gs[0,1])
    
        f_ax2.plot(np.arange(nPlot2), resMH_sumr/(itertot-burnin), 
                   marker='o', markersize=3, markeredgecolor='black', 
                   color='white', linestyle='none')
        
        f_ax2.axhline(y=threshold_func, color='red', linewidth=1)
        f_ax2.set_ylabel('Posterior probabilities') # Posterior probabilities
        f_ax2.set_title('Functions selection') # Functions selection

        ### ax3
        f_ax3 = fig_a.add_subplot(gs[1, :])
        simu_01_ = simu_01  # * simu_01_sd + (simu_01_mean)
        reconstruction_ = reconstructiontot # .reshape(-1,1)
        # reconstruction_ = reconstruction_ * simu_01_sd + simu_01_mean

        f_ax3.plot(np.arange(nPlot1), simu_01_, ls='solid', 
                   lw=2, label='True serie {}'.format(title)) # True serie
        f_ax3.plot(np.arange(nPlot1), reconstruction_, color='#0f3cf1', lw=2, 
                   ls=(0, (5, 5)), label='Reconstruction') # Reconstruction    
    
        if not(reconstructionalt is None):
            f_ax3.plot(np.arange(nPlot1)[0:len(reconstructionalt)], reconstructionalt, color='#AF601A', lw=2, 
                    ls=(0, (5, 5)), label='Reconstruction') # Reconstruction  

        for i in breakpoints[1:]:
            f_ax3.axvline(x=i, color='#D35400', linewidth=1.5)
    
        # where some data has already been plotted to ax
        handles, labels = f_ax3.get_legend_handles_labels()
        #Break points
        blue_line = mlines.Line2D([],[], color='#873600', label='Breakpoints') 

        # handles is a list, so append manual patch
        handles.append(blue_line) 

        f_ax3.legend(handles=handles)
        f_ax3.set_title('Serie {}'.format(title)) # True serie and reconstruction
        f_ax3.set_ylabel('Levels') # Chilean Pesos
        f_ax3.set_xlabel('Time') # day

    if save_fig == True:
        fig_a.savefig(path+'/res_P0i_{}.png'.format(title), dpi=400)
        
    