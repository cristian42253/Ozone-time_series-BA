BaycnTSA 
==========

A Bayesian approach for detecting change points in time series with missing data

References
==========

1.Bayesian approach for ozone time series segmentation in presence of missing data.

## Installation

Create bases

```sh
Rscript bases_generator.R data/samples/processed_o3_data.csv data/outputs/
```

Create environtment
```sh
python3 -m venv env
```
Activate environtment
```sh
source env/bin/activate
```
Install the dependencies.
```sh
pip install -r requirements.txt
```