# Ensemble Kalman Square Root Filter

Implementation of the Ensemble Kalman Square Root Filter using BLAS. The algorithm is after [Whitaker and Hamill (2002)](https://journals.ametsoc.org/view/journals/mwre/130/7/1520-0493_2002_130_1913_edawpo_2.0.co_2.xml) and [Peters et al. (2005)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2005JD006157). That said, the package is experimental.

## Installation

### python api 

clone the repository and cd into the main directory
Build with
```sh
pip install . 
```
or use
``` 
python setup.py build
```
in the latter case, make sure that the resulting shared object file is located in your python path.

#### Deprecation of distutils

The build currently relies on the distutils python packages which is deprecated, and support will [drop somewhere in 2023](https://numpy.org/devdocs/reference/distutils_status_migration.html). Either meson or cmake is the solution, but not currently implemented.

### Fortran

Currently, no automatic installation of either a static library nor shared library exist. To include it in a project, copy the source. 

## Usage

### python


