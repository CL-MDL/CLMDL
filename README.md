Implementation of CLMDL
================

This GitHub repository contains R functions to implement the CLMDL algorithm proposed in the paper arXiv:1904.06340. The CLMDL algorithm is designed for conducting multiple change-point estimation in a spatio-temporal process based on a parametric model (e.g. a linear regression based mean function + a Matern covariance function).

The repository consists of three folders, *code_demo*, *code_manuscript* and *output*. The *code_demo* folder contains R codes for a demonstrative example of CLMDL, the *code_manuscript* folder contains R codes for simulation 1 and 2 and real data analysis in the paper arXiv:1904.06340, and the *output* folder contains the simulated data and the corresponding RData files generated by the R codes. We give a detailed description of the code_demo, code_manuscript and output folders, and we then give the workflow for implementing the code.


### The code folder

The code folder has three R files, `DGP.R`, `Main. R` and `Util_Functions.R`.
1. The `DGP.R` file contains the R code to simulate a spatio-temporal process based on the four-parameter autoregressive spatial model (see equation (22) in the manuscript), which is the parametric model employed by CLMDL for the numerical studies in the paper. In particular, on each stationary segment, we have $$y_t-\mu = \phi (y_t-\mu) +\epsilon_t,$$ where $y_t=(y_{t,\textbf{s}}, \textbf{s}\in \mathcal S)\in \mathbb R^S$ collects the observations on $S$ spatial locations at time $t$, and $\mu$ is the mean parameter, $\phi$ is the temporal dependence parameter, and $\epsilon_t=(\varepsilon_{t,\textbf{s}}, \textbf{s}\in \mathcal S)$ is a Gaussian process with exponential covariance function ${\rm Cov}(\varepsilon_{t,\textbf{s}}, \varepsilon_{t,\textbf{s}'}) = \sigma^2\exp(-\left\lVert \textbf{s}-\textbf{s}'\right\lVert_2/\rho)$ with $\sigma^2$ being the variance parameter and $\rho$ being the spatial dependence parameter. This is a special case of the Matern covariance function with $\nu=0.5$. In total, we have four parameters $(\phi, \rho, \sigma^2, \mu)$ on each stationary segment.

   The `DGP.R` file generates the above four-parameter autoregressive spatial model on a regular two-dimensional rectangle grid $\mathcal S=(1,2,\cdots,s)^2$ of spatial dimension $S=s^2$ and time dimension $T$. This is the setting used in the numerical studies of the paper. It returns a data matrix of size $T \times S$ that stores the spatio-temporal process and a distance matrix of size $S\times S$ that stores the distance between the $S$ spatial locations.

2. The `Main.R` file contains the R code to implement the CLMDL algorithm based on the four-parameter autoregressive spatial model. In particular, it requires the input of (a) a data matrix of size $T \times S$ that stores the spatio-temporal process and (b) a distance matrix of size $S\times S$ that stores the distance between the $S$ spatial locations. It then conducts multiple change-point estimation in the spatio-temporal process and further constructs confidence intervals for each detected change-point (if any). The key tuning parameters are `s.lag` and `t.lag`, which specifies which pairs of observations are used in the construction of the composite likelihood in CLMDL. The `s.lag` is the spatial distance parameter $d$ in the manuscript and the `t.lag` is the time lag parameter $k$ in the manuscript. Another tuning parameter is minimum spacing between two change-points, we set it as $0.1 T$.

   We give a description of several important functions used in `Main.R`.
   
   2.1. Based on the distance matrix, `s.lag` and `t.lag`, the function `D.cal()` computes the average number of times that an observation is used in the composite likelihood (see equation (7) in the manuscript) and the number of times a marginal likelihood is needed for correcting the edge effect.

   2.2. The `pl()` function computes the proposed composite likelihood (including both the pairwise likelihood and the marginal likelihood) based on the four-parameter space autoregressive model for a given segment of the spatio-temporal process (see equation (5) in the manuscript).

   2.3. The key function is `pelt()`, which implements the PELT algorithm to minimize the proposed objective function CLMDL (see equation (8) in the manuscript) and returns its minimizer, which is the change-point estimator. The PELT algorithm is a modified dynamic programming algorithm with an additional pruning step. For every time point $t=1,2,\cdots, T$, the `pelt()` function outputs the current change-point estimator along with the (pruned) candidate for the dynamic programming at the time point $t$. To further improve computational efficiency, we further include a `res` argument in the `pelt()` function, which stands for *resolution*. By default, `res=1`, which means that the `pelt()` function considers all $t=1,2,\cdots T$ as potential change-points, this is the algorithm implemented in the manuscript. For computational efficiency, we can also set `res=2` when $T$ is large, in such a case, the `pelt()` function will only consider even time points as potential change-points. We implement both `res=1` and `res=2` in the `Main.R` file for comparison.

   2.4. The `Q.grid()` function is used to construct the confidence interval of the estimated change-point based on a parametric bootstrap procedure (see equation (18) in the manuscript).

4. The `Util_Functions.R` contains auxiliary functions used in the `DGP.R` and `Main.R` files.

### The output folder

The output folder has three RDS files, `y0.RDS`, `y1.RDS` and `y2.RDS`, which are generated by the `DGP.R` file and contains a spatio-temporal process with 0, 1, and 2 change-points, respectively. The folder also has three RData files, `Result_y0.RData`, `Result_y1.RData` and `Result_y2.RData`, which are generated by the `Main.R` file and contains the estimated change-points and constructed confidence intervals based on CLMDL for the dataset `y0.RDS`, `y1.RDS` and `y2.RDS`.

### The overall workflow

1. Download all three R files `DGP.R`, `Main. R` and `Util_Functions.R` in the code folder and put them into one single folder.
2. Run `DGP.R`, which will generate three simulated spatio-temporal processes `y0.RDS`, `y1.RDS` and `y2.RDS` with 0,1,2 change-points, respectively
3. Run `Main.R`, which will run the proposed CLMDL algorithm to conduct multiple change-point estimation on the three simulated spatio-temporal processes and save the analysis results into `Result_y0.RData`, `Result_y1.RData` and `Result_y2.RData` files, respectively. It takes around 90 mins on a desktop to analyze each simulated data using the `Main.R` file.

