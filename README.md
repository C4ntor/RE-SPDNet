# ReSPDNet
A Geometric Deep Learning Approach to
 Forecast the Time Series of Covariance Matrices.


Note that the copyright of the manopt toolbox is reserved by https://www.manopt.org/  
All rights of SPDNet code Version 1.0 is reserved by Zhiwu Huang and Luc Van Gool. A Riemannian Network for SPD Matrix Learning, In Proc. AAAI 2017. 

## Installation
Clone the repository 
```bash
    git clone https://github.com/C4ntor/RE-SPDNet
    cd RE-SPDNet
```

### Usage

1. Launch  
```bash
    Simulation_pub.R
```
to generate data and perform predictions of alternative models

2. Save locally the simulated data (from R workspace) in the RE-SPDNET/data folder. Ensure the file is .csv of vectorized covariance matrix, where columns correspond to distinct elements of the matrix, and rows to time observations. First row contains headers

3.  Launch  
```bash
    spdnet_afew.m
```
Ensure the project path, and data filename have been changed accordingly.

