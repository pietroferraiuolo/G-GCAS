
# Gaia - Globular Clusters Analysis Software

A brief description of what your project does.

## Table of Contents

- [Installation](#installation)
- [Package](#package-content-description)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

Instructions on how to install and set up your project.

```bash
# Example command
git clone https://github.com/yourusername/yourproject.git
cd yourproject
npm install
```

# Package Content Description

### king.f90
This Fortran program calculates the equilibrium configurations of a gravitational gas using the King distribution function as described in the 1966 article by I.P. King (King, I.P. (1966). AJ 71, 64.). The program is designed to analyze isotropic King models and can be modified to handle anisotropic Boltzmann functions.

### cluster.py
This module contains the `Cluster` class, which encapsulates all the information of a specified cluster. The class is designed to load and manage parameters of globular clusters from the Harris Catalogue 2010 Edition.

### plots.py
This module contains various plotting functions to visualize the kinetic data of globular clusters. It includes functions for creating scatter plots with histograms, color-magnitude diagrams, proper motion distributions, spatial distributions, histograms with optional kernel density estimation, scatter plots with error bars, and more. Each function is designed to be flexible and customizable through various parameters.

#### Functions

- `scatter_2hist(x, y, kde=False, **kwargs)`: Creates a 2D scatter plot with histograms projected on each axis.
- `colorMagnitude(g, b_r, teff_gspphot, **kwargs)`: Generates a color-magnitude diagram using BP and RP photometry and temperature information.
- `properMotion(pmra, pmdec, **kwargs)`: Plots the proper motion distribution of the sample.
- `spatial(ra, dec, **kwargs)`: Creates a scatter plot in the RA-Dec plane.
- `histogram(data, kde=False, **kwargs)`: Plots a histogram of the data with optional kernel density estimation.
- `scat_xhist(x, y, xerr=None, **kwargs)`: Makes a scatter plot with a projected histogram relative to a quantity 'y'.
- `errorbar(data, dataerr, x=None, xerr=None, **kwargs)`: Creates a scatter plot with error bars.

Each function provides additional parameters for customization, allowing users to tailor the plots to their specific needs.

### query.py
This module contains the `GaiaQuery` class, which facilitates querying the Gaia database using ADQL (Astronomical Data Query Language). The class provides methods for performing custom queries as well as pre-defined queries for astrometric, photometric, and radial velocity data. It also includes functionality for saving query results and their associated metadata.

#### Methods

- `free_query(radius, gc=None, save=False, **kwargs)`: Executes a custom ADQL query with specified parameters and conditions.
- `get_astrometry(radius, gc=None, save=False, **kwargs)`: Retrieves principal astrometric parameters with optional query conditions.
- `get_photometry(radius, gc=None, save=False, **kwargs)`: Retrieves principal photometric parameters with optional query conditions.
- `get_rv(radius, gc=None, save=False, **kwargs)`: Retrieves radial velocity parameters with optional query conditions.

Each method allows for saving the query results and their metadata, providing a comprehensive tool for Gaia data retrieval and analysis.

### functions.py

This module provides functions to compute various astronomical and physical quantities, such as angular separation, line-of-sight distance, radial distances in 2D and 3D, total velocity, and effective gravitational potential.

#### Functions

- `angular_separation(ra0=None, dec0=None)`: Calculates the angular separation between two points in the sky.
- `los_distance()`: Calculates the line-of-sight distance based on parallax.
- `radial_distance_2d(analytical_w=False, **params)`: Returns the formula for the 2D-projection of the radial distance of a source from the cluster center or given RA/DEC coordinates.
- `radial_distance_3d(gc_distance=None, analytical_r2d=False, analytical_w=False)`: Calculates the 3D radial distance of a source from the cluster center.
- `total_velocity()`: Calculates the total velocity based on x and y components.
- `effective_potential(shell=False)`: Calculates the effective gravitational potential.

### analyzers.calculus.py
This module provides functions for numerical computation and error analysis in the context of globular cluster studies. It includes methods for computing numerical values of functions, error propagation, and integration using the Gauss-Legendre method. Additionally, it offers an interface to a Fortran90 code for integrating the Single-Mass King model.

#### Functions

- `compute_numerical_function(func, variables, var_data)`: Computes the numerical value of a sympy function given its variables and associated data.
- `compute_error(func, variables, var_data, var_errors, corr=False, corr_values=None)`: Computes the numerical error for a given function using error propagation formulas.
- `error_propagation(func, variables, correlation=False)`: Computes the analytical error of a function using standard error propagation techniques.
- `gaus_legendre_integrator(fnc, a, b, points)`: Integrates a function using the Gauss-Legendre method.
- `king_integrator(w0, output='profile')`: Calls the king.f90 code to integrate the Single-Mass King model and returns the output data.

Each function is designed to handle specific computational tasks, providing robust tools for the analysis of globular clusters.

### analyzers.dynamics.py
This module provides functions for analyzing the dynamics of globular clusters. It includes methods for converting proper motion to velocity, computing 2D and 3D radial distances, total velocity, and density profiles.

#### Functions

- `velocity_conversion(mu, gc_distance, mu_error=0, gc_distance_error=0)`: Converts proper motion into velocities in km/s, with optional error propagation.
- `compute_r2d()`: Computes the 2D radial distance of a source from the cluster center.
- `compute_r3d()`: Computes the 3D radial distance of a source from the cluster center.
- `compute_Vtot()`: Calculates the total velocity of a source.
- `density_profile(data)`: Computes the density profile of the cluster based on the provided data.

Each function is designed to facilitate the analysis of globular cluster dynamics, providing essential tools for astrophysical research.

### analyzers.effective_potential.py

### analyzers._glpoints.py
This module provides the Gauss-Legendre quadrature points and weights for numerical integration. It includes predefined sets of points and weights for 20, 40, 80, and 96-point Gauss-Legendre quadrature.

#### Variables

- `x20`, `w20`: Points and weights for 20-point Gauss-Legendre quadrature.
- `x40`, `w40`: Points and weights for 40-point Gauss-Legendre quadrature.
- `x80`, `w80`: Points and weights for 80-point Gauss-Legendre quadrature.
- `x96`, `w96`: Points and weights for 96-point Gauss-Legendre quadrature.

These variables can be used to perform high-precision numerical integration over a specified interval.

### utility.osutils.py
This module provides utility functions for file and data management within the Gaia - Globular Clusters Analysis Software. It includes functions for loading data, retrieving file lists, generating timestamps, and searching for tracking numbers.

#### Functions

- `load_query(file)`: Loads data from a file into an Astropy QTable.
- `get_kwargs(names, default, kwargs)`: Retrieves a value from keyword arguments, returning a default if not found.
- `get_file_list(tn=None, fold=None, key=None)`: Returns a sorted list of files from a specified directory, with optional filtering by a key.
- `tnlist(gc_name)`: Lists tracking numbers for a given globular cluster.
- `_timestamp()`: Generates a timestamp in the format 'yyyymmdd_HHMMSS'.
- `_findTracknum(tn, complete_path=False)`: Searches for a tracking number within data path subfolders, returning paths or folder names.

Each function is designed to streamline data handling and file operations, enhancing the overall functionality of the software.

### utility.folder_paths.py
This module defines various folder paths used throughout the Gaia - Globular Clusters Analysis Software. It includes paths for base data, catalog files, query data, models, and untracked data. Additionally, it provides functions to generate specific paths for cluster data and models based on cluster names.

#### Variables

- `BASE_PATH`: The base path for the data, retrieved from the `GCASDATA` environment variable.
- `BASE_DATA_PATH`: The main data directory.
- `CATALOG_FILE`: Path to the catalog file.
- `QUERY_DATA_FOLDER`: Directory for query data.
- `KING_MODELS_FOLDER`: Directory for King models.
- `KING_INTEGRATOR_FOLDER`: Directory for King integrator files.
- `UNTRACKED_DATA_FOLDER`: Directory for untracked data.

#### Functions

- `CLUSTER_DATA_FOLDER(name: str)`: Returns the path for a specific cluster's data folder.
- `CLUSTER_MODEL_FOLDER(name: str)`: Returns the path for a specific cluster's model folder.

## Usage

Instructions and examples on how to use your project.

```python
# Example usage
import yourmodule

yourmodule.do_something()
```

## Contributing

Guidelines for contributing to your project.

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes.
4. Commit your changes (`git commit -am 'Add new feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a new Pull Request.

## License
MIT License