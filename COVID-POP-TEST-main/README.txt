Code for: Optimization of population-level testing and contact tracing for the next phase of the COVID-19 pandemic: a modelling study based on data from China

## Citation

 Optimization of population-level testing and contact tracing for the next phase of the COVID-19 pandemic: a modelling study based on data from China

Zengmiao Wang, Bingying Li, Jose Lourenco, Ruixue Wang, Yidan Li, Yonghong Liu, Ottar N. Bjornstad, Huaiyu Tian.

## Abstract

The COVID-19 pandemic entered a new phase in 2022 because of the widely distributed vaccines and the emerged variants. However, the decreased vaccines efficiency and asymptomatic transmission add great uncertainty to the future of the pandemic. Control measures are still necessary to contain the flare-ups of SARS-CoV-2. Here, we investigate the effectiveness of population-level testing (PLT) strategies (including a combination of testing interval, testing frequency, break interval for testing crew, and response time) and contact tracing (CT) to curb COVID-19 resurgence using a data set that includes both symptomatic and asymptomatic infections detected from multiple rounds of population-level PCR tests and CT in China during the flare-ups in 2021. We show that resurgent in a city could be controlled by 3–7 rounds of PLT, supplemented by CT as long as testing is frequent. However, data and models suggest that the window for control narrows for resurgence triggered by more contagious variants such as lineage B.1.1.7 (Alpha), P.1 (Gamma), B.1.617.2 (Delta) and B.1.1.529 (Omicron). Four rounds of PLT would still be needed to curb the resurgence induced by Gamma and Delta even with more than 50% of vaccine coverage. The needed rounds would be 7 for Omicron. Higher vaccination coverage would reduce significantly the number of rounds necessary. Our analysis potentially provides examples to control future flare-ups of SARS-CoV-2 via PLT and CT, even when immunity coverage has been established through natural infection or vaccines.

## Notes on the code

To run, you need a Matlab toolbox called "DRAM": DRAM is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. This page explains the basic ideas behind DRAM and provides examples and Matlab code for the computations.(see http://helios.fmi.fi/~lainema/dram/)

About code folder:

About code folder: We used three code folders to estimate the parameters for Tonghua, Changchun and Xingtai, respectively. For each folder: main file 1 (fitting the model): f1.m; function dependencies:f2.m,f3.m,f4.m,ContactTracing.m; The setting of the fixed parameters for each city can be found in the main text and supplemented material.

The code to simulate the scenarios can be found in the Tonghua folder: main file 2: f_simu2.m; function dependencies:fsimu_deterministic.m,fsimu_P_deterministic.m,fsimu_R_deterministic.m,fsimu_MiddleSize_deterministic_v2.m,fsimu_runFromMiddle_deterministic_v2.m,fsimu_Vaccine_deterministic_asym_v2.m,fsimu_Vaccine_deterministic_v2.m; Note that do not run this file directly since running time is quite long.




## Data

### Epidemiological data

We collected the daily official case reports from the health commission of 2 provincial-level administrative units and 3 city-level units, the website’s links are provided. The information was collected by Bingying Li.

### Population-level testing strategy

We prepared the population-testing strategy files (for example: population_testing_strategy_break0.xlsx. "break0" represents the break interval is 0 day.). Each column represents one popualtion-level testing strategy.



