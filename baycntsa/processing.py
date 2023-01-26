####################################################################################
# main script
####################################################################################

import sys
params = sys.argv

if len(params) < 4:
    print("At least two argument must be supplied (input file), (bases), (output folder) \n")
    sys.exit(100)

INPUT_FILE = params[1]
BASE_FILE = params[2]
OUTPUT_FILE = params[3]

import numpy as np
from functions import * 
from settings import * 

import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, title

## carga de la serie de datos
if DEBUG: print("\u001b[36m\n** Loading data \u001b[0m")
data_fnct_dp = pd.read_csv(INPUT_FILE, decimal='.', delimiter=',')
## cambio de formato a las fechas
data_fnct_dp['date'] = pd.to_datetime(data_fnct_dp['date'], format='%Y-%m-%d')
## cambios en el dataset
# se ordena la serie y se resetea el index
data_fnct_dp = data_fnct_dp.sort_values(by='date')
data_fnct_dp = data_fnct_dp.reset_index(drop=True)   
### listado de series
data_serie_list = []
for col in range(0, 5):   ### cantidad de series
    data = data_fnct_dp.iloc[:, col:col + 1]
    data_serie = np.array(data)
    data_serie_list.append(data_serie)    

if DEBUG: print("\u001b[36m** Draw raw series \u001b[0m")
### pintar los datos
for index, data in enumerate(data_serie_list): 
    # plt.ioff()
    fig = figure(figsize=(20, 5), dpi=200)
    plt.style.use(plt_style)
    plt.plot(data_fnct_dp['date'], data, color='black', linestyle='solid', linewidth=2, label='')
    plt.savefig(f"{OUTPUT_FILE}/fig_Pi0_{data_fnct_dp.columns[index]}")
    plt.close(fig)

### se escalan los datos con media=0 y sd=1
data_serie_means = []  ## medias de cada serie
data_serie_sd = []   ## sd de cada serie
pX = []
#pix = pd.read_csv('./inputs/Pi_0.csv', decimal='.', delimiter=',')
#P0_i_list = []
#for col in range(0, 5):   ### cantidad de series
#    P0_i = pix.iloc[:, col:col + 1]
#    data_P0_i = np.array(P0_i)
#    P0_i_list.append(P0_i)   
for index, a_serie in enumerate(data_serie_list):
    serie_without_nan = a_serie[~np.isnan(a_serie)]   ## remove NAs
    data_serie_means.append(np.mean(serie_without_nan)) ### media
    data_serie_sd.append(np.std(serie_without_nan))   ## sd
    data_serie = (serie_without_nan - np.mean(serie_without_nan)) / np.std(serie_without_nan)
    data_serie_list[index] = data_serie.reshape(len(data_serie), 1) ## save data without NAs
    
    # data_serie_list[index] = a_serie.reshape(len(a_serie), 1) ## save data with NAs
    pI = np.repeat(0.01, len(a_serie))
    for idx, val in enumerate(a_serie):
        if (idx+1) < len(a_serie) and not(np.isnan(val)) and np.isnan(a_serie[idx+1]):
            #print("nan-----> ", idx)
            pI[idx] = 0.001
        elif (idx+1) < len(a_serie) and np.isnan(val) and not(np.isnan(a_serie[idx+1])):
            #print("nan> ", idx+1)
            pI[idx+1] = 0.001
            pI[idx] = np.nan
        #else:
            #print("<<<<< ", idx)
        #    pI[idx] = 0.01
    pI[0] = 1.0
    pI = pI[~np.isnan(pI)]
    #print(pI)
    pX.append(pI)

result_list = []
for index, a_serie in enumerate(data_serie_list):
    
    n = len(a_serie) 
    # Probabilidad de ser punto de cambio para cada punto de la serie (0.01)
    # Cambias las posiciones 0.001
    # Pi = np.concatenate([np.array([1.00]), np.repeat(0.01, n-1)], axis=0)
    Pi = pX[index]
    data_fnct = pd.read_csv(BASE_FILE, decimal='.', delimiter=',') 

    data_fnct = data_fnct.iloc[:, 1:]
    Fmatrix = np.array(data_fnct)
    Fmatrix.shape

    # Probabilidad que una funcion sea seleccionada
    eta = np.concatenate([[1], np.repeat(0.01, Fmatrix.shape[1]-1)])

    start=datetime.now()
    if DEBUG: print(f"started at = {str(start)}")
    result_list.append(
        dbp_with_function_effect( Fmatrix, a_serie, itertot, burnin, lec1, lec2,
                        nbseginit, nbfuncinit, nbtochangegamma, nbtochanger, Pi, eta,
                        threshold_bp, threshold_fnc, printiter=False)
    )
    crono = datetime.now() - start
    if DEBUG: print(f"ended at = { crono/60 }\n")   

##
# Plot Result
##
mpl.rcParams.update(mpl.rcParamsDefault)

if DEBUG: print("\u001b[36m\n** Ploting results \u001b[0m")
##
# plot result
##
## plt.style.use("classic")
for index, result_ in enumerate(result_list):
    resMH = result_[0]
    # Puntos de cambio segn umbral
    breakpoints_bp = np.where(resMH["sumgamma"]/(itertot-burnin) > threshold_bp)[0]

    a_serie = data_fnct_dp.iloc[:, index:index + 1] 

    print("breakpoints_bp ", breakpoints_bp)
    translations = translation(np.array(a_serie))
    breakpoints_bp = [translations[val] for val in breakpoints_bp]
    if DEBUG: print("breakpoints_bp translated", breakpoints_bp)


    title = data_fnct_dp.columns[index]
    
    pickle.dump( result_, open( f'{OUTPUT_FILE}/result_Pi02_{title}.pkl', "wb" ) )
