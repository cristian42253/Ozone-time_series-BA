# -*- coding: utf-8 -*-

# Load libraries
# %%
import numpy as np
from functions import * 
from config import * 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot
import pickle

from sklearn.cluster import KMeans

# %%

import os
# cwd = os.getcwd()
# print(cwd)
results_pkl=[]
all_station = ["O3_BA" ,"O3_Comp","O3_ERA","O3_Flora","O3_Pance"]
for a_station in all_station:
    file = open(f"./objects/result_Pi0_{a_station}.pkl",'rb')
    results_pkl.append(pickle.load(file))

# %%
## carga de la serie de datos
if DEBUG: print("loading data")
data_fnct_dp = pd.read_csv('./inputs/processed_o3_data.csv', decimal='.', delimiter=',')
## cambio de formato a las fechas
data_fnct_dp['date'] = pd.to_datetime(data_fnct_dp['date'], format='%Y-%m-%d')

# Cambios en el Dataset
# # se ordena la serie y se resetea el index
data_fnct_dp = data_fnct_dp.sort_values(by='date')
data_fnct_dp = data_fnct_dp.reset_index(drop=True)   
# Todas las series
data_serie_list = []
for station in all_station:   ### cantidad de series
    #data = data_fnct_dp.iloc[:, col:col + 1]
    data = data_fnct_dp.loc[:, station]
    data_serie = np.array(data)
    data_serie_list.append(data_serie)    

data_serie_means = []
data_serie_sd = []   ## sd de cada serie
for index, a_serie in enumerate(data_serie_list):
    serie_without_nan = a_serie[~np.isnan(a_serie)]   ## remove NAs
    data_serie_means.append(np.mean(serie_without_nan)) ## means by serie
    data_serie_sd.append(np.std(serie_without_nan)) # sd by serie

    data_serie = (serie_without_nan - np.mean(serie_without_nan)) / np.std(serie_without_nan)
    data_serie_list[index] = data_serie.reshape(len(data_serie), 1) ## save data without NAs
    # data_serie_list[index] = a_serie.reshape(len(a_serie), 1) ## save data with NAs


# %%
## Load config 
from config import *

# umbrales de corte para cada serie
threshold_bp_list = [0.45, 0.4, 0.5, 0.45, 0.55]
threshold_fnc_list = [0.3, 0.2, 0.4, 0.25, 0.4]

# TODO: Improve the conditions with assertion, is more polite
if len(threshold_bp_list) != len(threshold_fnc_list):
    f"Should consider an equal number of thredsholds"

if len(all_station) != len(threshold_bp_list):
    f"Should consider an equal number of thredsholds and stations"

for index, result_ in enumerate(results_pkl):

    print(">>>>>>>>>", f"Station_{all_station[index]}", 
            f"threshold_bp_{threshold_bp_list[index]}",
             f"threshold_fnc_{threshold_fnc_list[index]}"
    )

    resMH = result_[0]
    breakpoints_bp = np.where(resMH["sumgamma"]/(itertot-burnin) > threshold_bp_list[index] )[0]
    print("breakpoints_bp ", breakpoints_bp)
    #a_serie = data_fnct_dp.iloc[:, index:index + 1] 
    a_serie = data_fnct_dp.loc[:, all_station[index]] 

    d_fnct = pd.read_csv(
        './inputs/Fmatrix.csv', decimal='.', delimiter=','
        ) 
    d_fnct = d_fnct.iloc[:, 1:]
    TFmatrix = np.array(d_fnct)
    TFmatrix.shape

    fmatrix_nan = pd.read_csv(
        f"./inputs/Fmatrix{all_station[index]}.csv", decimal='.', delimiter=','
    ) 
    fmatrix_nan = fmatrix_nan.iloc[:, 1:]

    ## type conversion
    IFmatrix= np.array(fmatrix_nan)
    IFmatrix.shape

    serie_without_nan = np.asarray( a_serie.dropna() )
    n = len(serie_without_nan)
    serie_w_nan = serie_without_nan.reshape(-1,1)

    new_estims = dbp_function_with_resMH(
        resMH, Fmatrix=IFmatrix, data_serie=serie_w_nan, 
        itertot=itertot, burnin=burnin, lec1=lec1, lec2=lec2,
        threshold_bp=threshold_bp_list[index], 
        threshold_fnc=threshold_fnc_list[index], 
        printiter=False
    ) 
    ## Draw incomplete
    draw_serie(resMH["sumgamma"], resMH["sumr"], itertot, burnin, 
        serie_without_nan, breakpoints_bp, data_fnct_dp['date'][:len(serie_without_nan)], 
        data_serie_means[index], data_serie_sd[index], 
        threshold_bp_list[index], threshold_fnc_list[index],
        new_estims[7], showDate=False,
        path='./outputs', title=f"noNAs_Pi0_{all_station[index]}"
    )

    basefunctions = new_estims[8]
    
    print(">>>>>>>>>", f"Station_{all_station[index]}")

    print("basefunctions ", basefunctions)
    basefunctions_names = [fmatrix_nan.columns[idx] for idx in basefunctions]
    print("basefunctions_names ", basefunctions_names)
    basefunctions =[ 
        idx for idx,name in
        enumerate(d_fnct) if name in basefunctions_names
    ]
    basefunctions = np.asarray(basefunctions)
    print("basefunctions ", basefunctions)

    ## Complete serie

    translations = translation(np.array(a_serie))
    breakpoints_bp = [translations[val] for val in breakpoints_bp]
    if DEBUG: print("breakpoints_bp translated", breakpoints_bp)

    reconstruction1 = reconstruction = dbp_extern_reconstruction(
        Fmatrix= TFmatrix,
        estim= [new_estims[4], new_estims[5]],
        breakpoints= np.asarray(breakpoints_bp),
        basefunctions= basefunctions
    ) 

    zeros = [0]*(len(a_serie)- len(result_[7]))
    resMH["sumgamma"] = list(resMH["sumgamma"]) + zeros
    resMH["sumgamma"] = np.array(resMH["sumgamma"])
    resMH["sumr"] = list(resMH["sumr"]) 
    resMH["sumr"] = np.array(resMH["sumr"])

    ## Draw complete
    draw_serie(resMH["sumgamma"], resMH["sumr"], itertot, burnin, 
        np.array(a_serie), breakpoints_bp, data_fnct_dp['date'], 
        data_serie_means[index], data_serie_sd[index],
        threshold_bp_list[index], threshold_fnc_list[index],
        reconstruction1, showDate=False, path='./outputs', 
        title=f"FULL_P0i_{all_station[index]}"
    )


#%%
# Change Plots and identifiying change point after missing values posterior probabilities.
# Clustering
for index, result_ in enumerate(results_pkl):
    
    resMH = result_[0]
    resMHg = resMH["sumgamma"]/(itertot-burnin)
    resMHs = resMH["sumr"]/(itertot-burnin)

    resMHg = resMHg.reshape(-1, 1)
    resMHs = resMHs.reshape(-1, 1)

    fig = figure(figsize=(20, 20), dpi=400)
    plt.xlabel('Time', fontsize=30)
    plt.ylabel('Probability', fontsize=30)
    plt.rc('xtick', labelsize=30)
    plt.rc('ytick', labelsize=30)
    plt.scatter(x=range(0, len(resMH["sumgamma"])), y=resMHg, 
         s=300, marker='o',color='black')
    plt.axhline(y=threshold_bp_list[index],color="red")
    x_rc = breakpoints_bp
    #y_rc = (resMH["sumgamma"][(resMH["sumgamma"]/(itertot-burnin))>0.6] )/(itertot-burnin)
    #y_rc = [
    #    resMHg[pos]
    #    for pos in breakpoints_bp
    #    ]
    plt.scatter(x=x_rc,y=y_rc,s=600 ,marker="^")
    plt.ylim(-0.05, 1.05)
    plt.xlim(-50, 1200)
    plt.savefig(f'outputs/sccater_sumgamma_{all_station[index]}.png')
    plt.close(fig)

    # Extrac values for step by step
    ((resMH["sumgamma"][(resMH["sumgamma"]/(itertot-burnin))>0.6] )/(itertot-burnin))[breakpoints_bp.index(361)]
#%%
    #
    
    fig = figure(figsize=(20, 20), dpi=400)
    plt.xlabel('Time', fontsize=30)
    plt.ylabel('Probability', fontsize=30)
    plt.rc('xtick', labelsize=30)
    plt.rc('ytick', labelsize=30)
    plt.scatter(x=x_rc, y=y_rc, s=400 ,cmap='viridis')
    plt.axhline(y=threshold_bp_list[index],color="red")
    plt.savefig(f'outputs/sccater_rc_{all_station[index]}.png')
    plt.close(fig)


    fig = figure(figsize=(40, 20), dpi=400)
    #plt.style.use(plt_style)
    plt.plot(data_fnct_dp['date'], data, color='black', linestyle='solid', linewidth=2, label='')
    plt.savefig("./outputs/fig_{}".format( data_fnct_dp.columns[index] ))
    plt.close(fig)


