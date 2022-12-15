# Bayes_BPD_NaN [Bayesian approach for detecting change points in time series with missing data]


  In this repository, you will find the implementation in Python and R for detecting change points in the Ozone time series with functional parts and the presence of missing data. The files reproduce the results obtained in the manuscript entitled "Ozone time series segmentation in presence of missing data using Bayesian approach"

## Data

Air quality monitoring stations are where ozone is measured; for this case, a calendar year between 2019 and 2021 was taken. 

## Goal 

We aim to detect change points with functional parts using a Bayesian approach in time series with missing data. For this reason, we create function dictionaries according to each time series. Each dictionary has fewer Haar functions. 

## Getting started

To replicate the results and run the files it is necessary to install the dependencies found in the ```requirements.txt``` file

We begin treatment as follows:

1. To run the models we must start with the creation of the dictionaries, for which we use the file ```bases_gen.r```

2. Change the addresses of the file containing the time series and function dictionaries, in ```processing.py```

3. The ```config.py``` file contains the parameters on which the algorithm is set to perform the runs. You can consider increasing or decreasing them depending on the problem to study; this will imply program execution time.
```python
############## initial parameters

itertot = 160000
burnin = 40000
lec1 =  50
lec2 =  50
nbseginit = 10
nbfuncinit = 12
nbtochangegamma = 1
nbtochanger = 1

plt_style  = "classic" #"ggplot" #"seaborn"

DEBUG = True
```
4. When obtaining the estimates, the program creates a file with extension ```name.pkl``` ; it contains the estimates and the runs used to obtain the reconstruction and the locations of the change points.

5. To plot and reconstruct the estimates to the complete data series, including the missing data, we use the file ```drawler_pkl.py```. In this program, we can change the thresholds for the choice of change points and functions that appear significant according to the posterior probabilities. An example is the following:

```python
threshold_bp_list = [0.45, 0.4, 0.5, 0.45, 0.55]

threshold_fnc_list = [0.3, 0.2, 0.4, 0.25, 0.4]
```

For each series, a threshold must be established for the change points and the choice of functions; in our case, we consider five Ozone series.
