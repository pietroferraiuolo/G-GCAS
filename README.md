# GRASP - Globular clusteR Astrometry and Photometry Software
 ![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pietroferraiuolo/GRASP/python-test.yaml)

The GRASP package is a tool for astrophysical data analysis, mainly thought for Globular CLusters 
astrometric and photometric data retrievement and GCs dynamical evolution analysis.

## Table of Contents

- [Installation](#installation)
    - [Build R](#build-r)
    - [ANTLR4](#install-latex-parser-dependencies)
    - [GRASP](#installing-grasp)
- [Examples](#retrieving-data)

## Installation
Since this package uses a combination of Python (mainly), R, C and Fortran, some additional steps, as
to avoid errors in the code execution, must be done before installing the package.

<ins>NOTE</ins>: The use of a conda environment with `python >= 3.10` is highly recommended for the
correct functioning of the package.

<ins>NOTE2</ins>: All the extra steps necessay will be accounted for in the setup of the package 
itself, sooner or later...

### Build R
Firstly, one must install R, as the `statistics` module of the package uses R code to perform regression routines:

```bash
sudo apt update
sudo apt install r-base r-base-dev
```

This will install the latest R distribution, as well as all the base packages that coome with it. The
`dev` install makes sure the appropriate compilers are installed, so that R packages can be correctly
built upon installing them. If compiler errors arise (mentioning `GLIBCXX_3.4.XX not found for librosa`), try installing directly the `gcc` compiler:

```bash
conda install -c conda-forge gcc
```

### Install latex parser dependencies: ANTLR4
Since the `formulary` module, which handles formulas definitions and computations, use both the sympy
and the latex syntax interpreter, the latter needs an additional package to run: the ANTLR4 python 
routine. As of version `0.2.0` of GRASP, the precise version to install is:

```bash
pip install antlr4-python3-runtime==4.11
```

Equivalently:

```bash
conda install -c conda-forge antlr4-python3-runtime==4.11
```

### Installing GRASP
Finally, when everything is set up, install the package through the repository, as no PyPI release
exists yet:

```bash
pip install git+'https://github.com/pietroferraiuolo/GRASP.git'
```

<details>
<summary>Utilization Examples</summary>

### Retrieving data
Right now, the only implemented archive available for data retrievement is the GAIA archive.
The GAIA archive is comprehensive of various data table, with the main table for data release `X`
being `gaiadrX.main_table`. To list all the available data tables:

```python
In [2]: import grasp

In [3]: grasp.available_tables() # or equivalentely grasp.gaia.query.available_tables()
    ```
INFO: Retrieving tables... [astroquery.utils.tap.core]
INFO: Parsing tables... [astroquery.utils.tap.core]
INFO: Done. [astroquery.utils.tap.core]
external.apassdr9
external.catwise2020
external.gaiadr2_astrophysical_parameters
.
. 
. # continuing with all available data tables
    ```
```

### Data visualization

### Computing formulas

</details>