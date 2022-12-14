# Bayes_BPD_NaN [Bayesian approach for detecting change points in time series with missing data]


  In this repository, you will find the implementation in Python and R for detecting change points in the Ozone time series with functional parts and the presence of missing data. The files reproduce the results obtained in the manuscript entitled "Ozone time series segmentation in presence of missing data using Bayesian approach"

## Data

Air quality monitoring stations are where ozone is measured; for this case, a calendar year between 2019 and 2021 was taken. 

## Goal 

We aim to detect change points with functional parts using a Bayesian approach in time series with missing data. For this reason, we create function dictionaries according to each time series. Each dictionary has fewer Haar functions. 

## Getting started
We begin treatment as follows:

1. To run the models we must start with the creation of the dictionaries, for which we use
   > bases_gen.R
