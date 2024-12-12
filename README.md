# Gaia - Globular Clusters Analysis Software
 ![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pietroferraiuolo/ggcas/python-package.yml)

This package is a software which uses GAIA data releases to retrieve and analyze Globular Cluster (mainly) astrometric data.

## Table of Contents

- [Installation](#installation)
- [Package](#package-content-description)

## Installation

As of now, that there is not yet a PyPi version available, intall the package from the github repository:

```bash
git clone https://github.com/pietroferraiuolo/G-GCAS.git
cd G-GCAS
pip install .
```

## Package Content Description

### functions.py
This module provides functions to compute various astronomical and physical quantities, such as angular
separation, line-of-sight distance, radial distances in 2D and 3D, total velocity, and effective 
gravitational potential. The classes which define those quantities provide both the analythical formula
and the methods for computing it, provided the data (and errors and computations if available).

### plots.py
This module contains various plotting functions to visualize the useful data of globular clusters. 
It includes functions for creating scatter plots with histograms, color-magnitude diagrams, proper 
motion distributions, spatial distributions, histograms with optional kernel density estimation, 
scatter plots with error bars, and more. Each function is designed to be flexible and customizable 
through various parameters.

### statistics.py
This module provides a series of functions for the statistical analysis of astronomical data. It 
includes functions for extreme deconvolution estimation, Gaussian mixture modeling, and regression 
model estimation. The module leverages R scripts for statistical computations and integrates them into 
Python using the `rpy2` package. R scripts and other supporting scripts are in the
`ggcas.analyzers._Rcode` module.

### _query.py
This module contains the `GaiaQuery` class, which facilitates querying the Gaia database using ADQL
(Astronomical Data Query Language). The class provides methods for performing custom queries as well
as pre-defined queries for astrometric, photometric, and radial velocity data. It also includes functionality
for saving query results and their associated metadata.

### _cluster.py
This module contains the `Cluster` class, which encapsulates all the information of a specified cluster. 
The class is designed to load and manage parameters of globular clusters from the Harris Catalogue 2010 Edition.

### analyzers
This module provides functions for numerical computation and error analysis in the context of globular
cluster studies. It includes methods for computing numerical values of functions, error propagation, and
integration using the Gauss-Legendre method. Additionally, it offers an interface to a Fortran90 code for
integrating the Single-Mass King model and an interface for the `McLuster` C/F code 
(Kuepper, A. H. W. ; Maschberger, Th. ; Kroupa, P. ; Baumgardt, H.).

#### dynamics.py
This module provides functions for analyzing the dynamics of globular clusters. It includes methods 
for converting proper motion to velocity, computing 2D and 3D radial distances, total velocity,
and density profiles.

#### calculus.py
This module contains various functions to compute numerically analythical equations (written with
sympy), as well as error formulae derivation and computation. It also has the Gauss-Legendre Quadrature
numerical integrator, supported by a sub-module, `_glpoints.py`, which contains predefined sets of points
and weights for 20, 40, 80, and 96-point Gauss-Legendre quadrature.

#### mcluster.py
This module is a bridge from python to the C/Fortran code `McLuster` and `McLusterSEE`. This code generates
syntetic data of globular clusters, while it's other, more complicated, version `SEE` includes stellar 
evolution when generaing the cluster. This code is often used to generate imput for numerical simulations
tools, such as Nbody6.

#### king.py
This module is a bridge between python and the fortran code `king.f90`. This Fortran program calculates
the equilibrium configurations of a gravitational gas using the King distribution function as described
in the 1966 article by I.P. King (King, I.P. (1966). AJ 71, 64.). The program is designed to analyze
isotropic King models and can be modified to handle anisotropic Boltzmann functions.

#### _Rcode
This sub-module contains R source code for statistical analysis, which are:
- `regression.R`: a code which performs regresison on data of many equations, like 'gaussian', 'linear', 'poisson', and so on...
- `gaussian_mixture.R`: a code which trains a Gaussian Mixture Model using ``Mclust``.
In support of bridging python to R, there are two modules, `r_check.py` and `r2py_models.py`, which handle
the connection and translation between the two languages.

### _utility
This module provides utility functions for file and data management within the Gaia - Globular Clusters
Analysis Software as well as defining supportive classes of better data/sample management.

#### folder_paths.py
This module defines various folder paths used throughout the Gaia - Globular Clusters Analysis Software.
It includes paths for base data, catalog files, query data, models, and untracked data. Additionally, it
provides functions to generate specific paths for cluster data and models based on cluster names.

#### base_formula.py
This module defines the base class inherited by all `functions.py` astronomical formulae.

#### sample.py
This module provides the `Sample` class, which encapsulates data samples retrieved from the Gaia database.
It includes methods for data manipulation, filtering, and basic statistical analysis, making it easier
to work with subsets of data within the G-GCAS framework.

#### osutils.py
This module provides utility functions for operating system interactions within the Gaia - Globular 
Clusters Analysis Software. It includes functions for handling file operations, directory management,
and environment variable manipulation. These utilities help streamline the workflow by providing a
consistent interface for common OS-related tasks.
