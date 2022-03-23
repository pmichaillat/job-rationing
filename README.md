# Do Matching Frictions Explain Unemployment? Not in Bad Times: Code and Data

This repository contains the code and data accompanying the paper "Do Matching Frictions Explain Unemployment? Not in Bad Times", written by [Pascal Michaillat](https://pascalmichaillat.org), and published in the [American Economic Review](https://doi.org/10.1257/aer.102.4.1721) in June 2012.

## Paper webpage

The paper and its online appendix are available at https://pascalmichaillat.org/1/.

## Data

The folder `data` contains text files with data from the Bureau of Labor Statistics (BLS) and other sources. The data describe key macroeconomic variables and are used in the simulations:

+ `CPI-URBAN.txt` - Urban consumption price index from BLS
+ `CPS-UL.txt` - Unemployment level from the BLS Current Population Survey (CPS)
+ `CPS-UR.txt` - Unemployment rate from CPS
+ `MSPC-EMP.txt` - Employment level from the BLS Major Sector Productivity and
Costs program (MSPC)
+ `MSPC-OUTPUT.txt` - Output level from MSPC
+ `HELPWANT.txt` - Help-wanted advertising index from the Conference
Board
+ `JOLTS-JOLNF.txt` - Job-opening level in the nonfarm business sector from
the BLS Job Opening and Labor Turnover Survey (JOLTS)
+ `CES-HWAGEPROD.txt` - Average hourly earnings of production and  nonsupervisory workers from the BLS Current Employment Survey (CES)
+ `ECI.txt` - Various measures of worker compensation from the BLS Employment Cost Index (ECI)
+ `FERNALD.txt` -  Utilization-adjusted total factor productivity from Fernald (2009)

The readme files `CES_README.txt`, `CPI_README.txt`, `CPS_README.txt`, `HELPWANT_README.txt`, `MSPC_README.txt`, `JOLTS_README.txt`, and `ECI_README.txt` provide additional details on the data.

## Code

The simulations are conducted with MATLAB.

### Helper scripts and functions

The simulations rely on a number of helper scripts and functions:

+ `setup_equilibrium.m`, `setup_elasticity.m`, `setup.m`, `setup_g.m`,  `setup_robusthigh.m`, `setup_robustlow.m` - Calibrate parameters used in simulations. Calibrated values are summarized in Table 1 and Table A6.
+ `QTOW.m`, `QUARTER.m`, `hpfilter.m`, `TECHNO_1600.m`, `TECHNO.m`, `TFP_1964_2009.m`, `W2QUARTER.m`, `data_1964_2009.m`, `bayesdata.m` - Fetch, prepare, and transform data
+ `FINDTH.m` - Perform useful calculations
+ `reduform.m`, `shftrght.m`, `vech.m`, `STEADYLL.m`, `aim_eig.m`, `aimerr.m`, `build_a.m`, `copy_w.m`, `eigsys.m`, `ex_shift.m`, `numshift.m`, `obstruct.m`, `penta2.m`, `LIN_DSGE.m`, `LIN_DSGE_c.m`, `LIN_DSGE_g.m`, `setupsimul.m`, `setupsimul_1600.m`, `setupsimul_c.m`, `setupsimul_g.m`, `setupsimul_robusthigh.m`, `setupsimul_robustlow.m` - Solve the loglinear DSGE model with job rationing and matching frictions
+ `SUMSTAT.m`, `ACF.m`, `AR.m`, `AUTOCORREL.m` - Perform statistical analysis
+ `EXPECTEDMC.m`, `MAKEMC.m`, `MCSOLVE.m`, `OBJEULER.m`, `SHOOTING.m`, `SIMULFT.m`, `SOLVESYS.m`, `STEADYGE.m` - Solve the nonlinear DSGE model using the Fair-Taylor algorithm


### Scripts for main-text results

The results of the numerical simulations in the main text are obtained by running the following scripts:

+ `diagram_equilibrium.m` - Draw the equilibrium diagram for various search-and-matching models. Technology and recruiting cost take different values. Results presented in Figure 1.
+ `diagram_elasticity.m` - Compute the elasticity of labor market tightness with respect to recruiting cost in the model with wage rigidity and in the model with job rationing. Results are presented in Figure 2.
+ `moments_US.m` - Compute the moments of key labor market variables using US data for the 1964–2009 period. US data are detrended using a Hoddrick-Prescott (HP) filter with parameter of 10000. Results are presented in Table 2.
+ `moments.m` - Compute the moments of key labor market variables simulated using the loglinear DSGE model with matching frictions and job rationing. Repeat the simulation of the model many times to obtain more precise estimates of the simulated moments. Report estimated second moments and standard deviation of these estimates. Results are presented in Table 3.
+ `irf_rationing.m` - Compute the impulse response functions of the loglinear DSGE model with job rationing and matching frictions. Results are presented in Figure 3.
+ `script_FT.m` - Solve the nonlinear DSGE model with job rationing and matching frictions, when it is subject to the actual technology shock measured in US data. Compare actual unemployment to simulated unemployment. Decompose the simulated unemployment series into a frictional and a rationing component. Results are presented in Figures 4 and 5.

### Scripts for online-appendix results

The results of the numerical simulations in the online appendix are obtained by running the following scripts:

+ `moments_robusthigh.m`, `moments_robustlow.m`, `moments_c.m`, `moments_g.m` - Compute the moments of key labor market variables simulated using the loglinear DSGE model with job rationing and matching frictions for different calibrations of the recruiting cost, a different specification of recruiting expenses, and gradual wage adjustment. Results are reported in Tables A1, A2, and A3.
+ `decomposition_shocks.m` - Decompose actual unemployment into frictional and rationing components using a variety of specifications for the type of shocks in the economy. Results are reported in Figures A1, A2, and A3.
+ `moments_US_1600.m` - Compute the moments of key labor market variables using US data for the 1964—2009 period; US data are detrended using a HP filter with conventional parameter of 1600. Results are presented in Table A4.
+ `moments_1600.m` - Compute the moments of key labor market variables simulated using the loglinear DSGE model with matching frictions and job rationing for a  HP-filter parameter of 1600 instead of 10000. Results are reported in Table A5.
+ `irf_robusthigh.m`, `irf_robustlow.m` , `irf_c.m`, `irf_g.m`  - Compute the impulse response functions of the loglinear DSGE model with job rationing and matching frictions for different calibrations of the recruiting cost, a different specification of recruiting expenses, and gradual wage adjustment. Results are reported in Figures A4, A7, and A8.
+ `script_FT_robusthigh.m`, `script_FT_robustlow.m`, `script_FT_TFP.m` - Solve the nonlinear DSGE model with job rationing and matching frictions when it is subject to the actual technology shock measured in US data. These scripts consider three variants from the analysis in the article: high recruiting cost, low recruiting cost, and capacity-adjusted technology series instead of measured technology series. The scripts compare actual unemployment to simulated unemployment and decompose the simulated unemployment series into a frictional and a rationing component. Results are reported in Figures A5, A6, A9, and A10.
+ `script_FT.m` - The end of the script compares the time series for unemployment and labor market tightness generated by the model with two different numerical solution methods: (i) a series of equilibria in static environments that abstract from aggregate shocks to technology and dynamics of unemployment; and (ii) the exact solution to the nonlinear model, which accounts fully for the dynamics of unemployment and rational expectations of stochastic process of technology and labor market variables. Results are presented in Figure A11.

## Software

The results were obtained using MATLAB R2010a on macOS Snow Leopard.

## License

This repository is licensed under the [MIT License](LICENSE.md).