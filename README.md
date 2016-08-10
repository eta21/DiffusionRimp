
# DiffusionRimp
Data-imputation and density approximations for diffusion processes.

## What is DiffusionRimp?
__DiffusionRimp__ is a package for performing data imputation on discretely observed diffusion processes as well as calculating numerical approximations to transition and first passage time densities.  

## Why use DiffusionRimp?
__DiffusionRimp__ provides routines for exploring and analysing diffusion processes with highly non-linear dynamics. Although the analysis of non-linear diffusion processes is often quite challenging, the __DiffusionRimp__ package makes it possible to analyse models with complicated dynamics with basic coding skills and experience in numerical analysis.  

## Get DiffusionRjgqd?
Check out [DiffusionRimp](https://github.com/eta21/DiffusionRimp) for the package source files, vignettes and other downloadable content.


# Installation Notes
Mac users may have to carry out some additional installation procedures in order for __DiffusionRimp__ to operate optimally. 

## Mac users:
To install the latest version of __Rcpp__, the latest version of R is needed.
To install __RcppArmadillo__, the __Fortran__ version used by R needs to be updated.
To install __rgl__, the computer needs to have X11 installed.
Update R to the latest version.
Run the following code:

```
install.packages("Rcpp", type = "source", dep = TRUE) 
```

#### Open a Terminal and run the following code:
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2 
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C / 

#### Back in R, run the following code:
```
install.packages("RcppArmadillo", dep = TRUE) 
```

#### Make sure you have X11 installed. 
Go to Applications/Utilities and see if X11 is there. If not, you’ll need to install X11 or XQuartz. These are available from http://xquartz.macosforge.org/landing/

#### Back in R, run the following code:
```
install.packages(“rgl", dep = TRUE) 
```

#### Download the DiffusionRjgqd package and run the code:
```
install.packages("~/DiffusionRimp_0.1.0.tar.gz", repos = NULL, type = "source”)
```

#### Run the following code in R to see if the package works:

```
library(DiffpackRimp) 
example(MOL.density)
```