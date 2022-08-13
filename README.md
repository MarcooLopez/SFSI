# SFSI: Sparse Family and Selection Indices

[![CRAN status](https://www.r-pkg.org/badges/version/SFSI?color=green)](https://CRAN.R-project.org/package=SFSI)
[![CRAN checks](https://cranchecks.info/badges/worst/SFSI)](https://cran.r-project.org/web/checks/check_results_SFSI.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SFSI)](http://www.r-pkg.org/pkg/SFSI)
[![Downloads/month](http://cranlogs.r-pkg.org/badges/SFSI?color=blue)](http://www.r-pkg.org/pkg/SFSI)

The SFSI R-package solves penalized regression problems offering tools for the solutions to penalized selection indices. In this repository we maintain the latest (developing) version.

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

## Selection Indices (SI)

Prediction of **breeding values** (u<sub><i>i</i></sub>) for a target trait (y<sub><i>i</i></sub>) is usually done using a **Selection Index (SI)**.
In the selection index all the available information contribute to the prediction of the *i*<sup>th</sup> candidate of selection as:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img1.png" height="26"/>
</p>

or (in matrix notation)

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img2.png" height="27"/>
</p>

where the predictors <b>x</b><sub><i>i</i></sub> can be indirect information from either:

- Correlated traits measured in the same candidates
- Measurements on the same trait of interest collected on related individuals

In the first case, the borrowing of information is provided by the **genetic covariance between the target trait and measured traits**. The second case is a kinship-based prediction approach since the borrowing of information is taken from **genetic relateness among individuals**. These relationships can be provided through either a pedigree- or marker-based relationship matrix (**genomic prediction**).

### Standard Selection Index

The weights <b>&beta;</b><sub><i>i</i></sub> = (&beta;<sub><i>i1</i></sub>,...,&beta;<sub><i>ip</i></sub>)'
are derived by minimizing the optimization problem:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img3.png" height="42"/>
</p>

This problem is equivalent to:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img4.png" height="33"/>
</p>

Under standard assumptions, the solution to the above problem is

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img5.png" height="28"/>
</p>

where <b>P</b><sub>x</sub> is the phenotypic variance-covariance matrix among predictors, <b>x</b><sub><i>i</i></sub>, and <b>G</b><sub>xy</sub> is a vector with the genetic covariances between predictors <b>x</b><sub><i>i</i></sub> and response y<sub><i>i</i></sub>.

### Penalized Selection Index
The regression coefficients can be derived by impossing a penalization in the above optimization function as

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img6.png" height="35"/>
</p>

where &lambda; is a penalty parameter and <i>F</i>(<b>&beta;</b><sub><i>i</i></sub>)
is a penalty function on the regression coefficients. A value of &lambda;=0 yields the coefficients for the standard (un-penalized) selection index. Commonly used penalty functions are based on the L1- and L2- norms.

* **L1-penalized Selection Index:** is obtained using the L1-norm:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img7.png" height="28"/>
</p>

This problem does not have a closed form solution; however, a solution can be obtained using Least Angle Regression (LARS) (Efron, 2004) or Coordinate Descent algorithms (Friedman, 2007). These algorithms are implemented in the SFSI R-package using <b>P</b><sub>x</sub> and <b>G</b><sub>xy</sub> as inputs.

* **L2-penalized Selection Index:** is obtained using the L2-norm:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img8.png" height="30"/>
</p>

In this case, the solution has the following closed form:

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img9.png" height="27"/>
</p>

where <b>I</b> is an identity matrix.

* **Elastic-Net-penalized SI:** considers a weighted penalization of both norms,

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/main/vignettes/Img10.png" height="30"/>
</p>

where &alpha; is a weighting parameter.

The L1-penalized and L2-penalized SI appear as special cases of the Elastic-Net-penalized index when &alpha;=1 and &alpha;=0, respectively.

**NOTE:** As for the L1-penalized SI, if &alpha;>0, no closed form solution exists. Solutions can be obtained using only the Coordinate Descent algorithm (Friedman, 2007). This method is implemented in the SFSI R-package using <b>P</b><sub>x</sub>, <b>G</b><sub>xy</sub> and &alpha; as inputs.

### Sparse Selection Index.
A Penalized Selection Index involving the L1-norm is refered to as **Sparse Selection Index** as the norm induces variable selection (sparsity in the coefficients).


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
