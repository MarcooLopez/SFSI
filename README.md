# SFSI: Sparse Family and Selection Indices

[![CRAN status](https://www.r-pkg.org/badges/version/SFSI?color=green)](https://CRAN.R-project.org/package=SFSI)
[![CRAN checks](https://badges.cranchecks.info/worst/SFSI.svg)](https://cran.r-project.org/web/checks/check_results_SFSI.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SFSI)](http://www.r-pkg.org/pkg/SFSI)
[![Downloads/month](http://cranlogs.r-pkg.org/badges/SFSI?color=blue)](http://www.r-pkg.org/pkg/SFSI)

The SFSI R-package solves penalized regression problems offering tools for the solutions to penalized selection indices. In this repository we maintain the latest (developing) version.

*Last update: Nov 17, 2023*

## Package installation

Installation of SFSI package requires a R-version greater than 3.5.0

From CRAN (stable version)
```r
 install.packages('SFSI',repos='https://cran.r-project.org/')
```

From GitHub (developing version)
```r
 install.packages('remotes',repos='https://cran.r-project.org/')  # 1. install remotes
 library(remotes)                                                 # 2. load the library
 install_github('MarcooLopez/SFSI')                               # 3. install SFSI from GitHub
```

## Selection Indices

Prediction of **breeding values** ($u_i$) for a target trait ($y_i$) is usually done using a **Selection Index (SI)**.
In the selection index all the available information contribute to the prediction of the $i^{th}$ candidate of selection as:

$$
\color{NavyBlue}{\mathcal{I}_ i= x_{i1}\beta_{i1} + x_{i2}\beta_{i2} + \cdots + x_{ip}\beta_{ip}}
$$

or (in matrix notation)

$$
\color{blue}{\mathcal{I}_ i = \boldsymbol{x}_{i}'\boldsymbol{\beta}_i}
$$

where the predictors $\boldsymbol{x}_ i= (x_{i1},...,x_{ip})'$ can be indirect information from either:

- Correlated traits measured in the same candidates
- Measurements on the same trait of interest collected on related individuals

### Standard Selection Index

The weights $\boldsymbol{\beta}_ i = (\beta_{i1},...,\beta_{ip})'$ are regression coefficients derived by minimizing the optimization problem:

$$
\color{blue}{\hat{\boldsymbol{\beta}}_ i = \text{arg min}\{\frac{1}{2}\mathbb{E}(u_ i - \boldsymbol{x}_{i}'\boldsymbol{\beta}_i)\}}
$$

This problem is equivalent to:

$$
\color{blue}{\hat{\boldsymbol{\beta}}_ i = \text{arg min}\[\frac{1}{2}\boldsymbol{\beta}'_ i\textbf{P}_ x\boldsymbol{\beta}_ i - \textbf{G}'_ {xy}\boldsymbol{\beta}_i\]}
$$

where $\textbf{P}_ x$ is the phenotypic variance-covariance matrix of predictors and $\textbf{G}_{xy}$ is a vector with the genetic covariances between predictors and response.

Under standard assumptions, the solution to the above problem is

$$
\color{blue}{\hat{\boldsymbol{\beta}}_ i = \textbf{P}^{-1}_ x\textbf{G}_{xy}}
$$

### Sparse Selection Index
The weights can be derived by impossing a sparsity-inducing penalization in the above optimization function as

$$
\color{blue}{\hat{\boldsymbol{\beta}}_ i = \text{arg min}\[\frac{1}{2}\boldsymbol{\beta}'_ i\textbf{P}_ x\boldsymbol{\beta}_ i - \textbf{G}'_ {xy}\boldsymbol{\beta}_i + \lambda f(\boldsymbol{\beta}_i)\]}
$$

where $\lambda$ is a penalty parameter and $f(\boldsymbol{\beta}_i)$ is a penalty function on the weights. A value of $\lambda = 0$ yields the coefficients for the standard (un-penalized) selection index. Commonly used penalty functions are based on the L1- (i.e., **LASSO**) and L2-norms (i.e., **Ridge Regression**). **Elastic-Net** considers a combined penalization of both norms,

$$
\color{blue}{f(\boldsymbol{\beta}_ i) = \alpha\sum^p_{j=1}|\beta_{ij}| + (1-\alpha)\frac{1}{2}\sum^p_{j=1}\beta^2_{ij}}
$$

where $\alpha$ is a number between 0 and 1. The LASSO and Ridge Regression appear as special cases of the Elastic-Net when $\alpha = 1$ and $\alpha = 0$, respectively.

Functions `LARS()` and `solveEN()` can be used to obtain solutions to the above penalized optimization problem taking $\textbf{P}_ x$ and $\textbf{G}_{xy}$ as inputs. The former function provides LASSO solutions for the entire $\lambda$ path using Least Angle Regression (Efron et al., 2004), and the later finds solutions for the Elastic-Net problem for given values of $\alpha$ and $\lambda$ via the Coordinate Descent algorithm (Friedman, 2007). 

## Documentation (two applications)
* **Application with high-throughput phenotypes:**
Lopez-Cruz *et al.* (2020). [[Manuscript](https://www.nature.com/articles/s41598-020-65011-2)]. [[Documentation](http://htmlpreview.github.io/?https://github.com/MarcooLopez/SFSI/blob/master/inst/doc/PSI-documentation.html)].

* **Application to Genomic Prediction:**
Lopez-Cruz and de los Campos (2021). [[Manuscript](https://doi.org/10.1093/genetics/iyab030)]. [[Documentation](http://htmlpreview.github.io/?https://github.com/MarcooLopez/SFSI/blob/master/inst/doc/SSI-documentation.html)].

## How to cite SFSI R-package
* Lopez-Cruz M, Olson E, Rovere G, Crossa J, Dreisigacker S, Mondal S, Singh R & de los Campos G **(2020)**. Regularized selection indices for breeding value prediction using hyper-spectral image data. *Scientific Reports*, 10, 8195.

A BibTeX entry for LaTeX is
```
  @Article{,
    title = {Regularized selection indices for breeding value prediction using hyper-spectral image data},
    author = {Marco Lopez-Cruz and Eric Olson and Gabriel Rovere and Jose Crossa and Susanne Dreisigacker
              and Suchismita Mondal and Ravi Singh and Gustavo {de los Campos}},
    journal = {Scientific Reports},
    year = {2020},
    volume = {10},
    pages = {8195},
  }
```

## Dataset
The SFSI R-package contains a reduced version of the full data used in Lopez-Cruz *et al.* (2020) for the development of penalized selection indices. This full data can be found in this [repository](https://github.com/MarcooLopez/Data_for_Lopez-Cruz_et_al_2020).

## References
* Efron B, Hastie T, Johnstone I & Tibshirani R **(2004)**. Least angle regression. *The Annals of Statistics*, 32(2), 407–499.
* Friedman J, Hastie T, Höfling H & Tibshirani R **(2007)**. Pathwise coordinate optimization. *The Annals of Applied Statistics*, 1(2), 302–332.
