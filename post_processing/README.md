# Pylbo
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pylo tests <master>](https://github.com/n-claes/legolas/workflows/pylbo-core/badge.svg?branch=master)](https://github.com/n-claes/legolas/actions?query=workflow%3Apylbo-core)

Pylbo (**Py**thon for **L**egolas **B**inary **O**utput) is the post-processing framework that we developed dedicated to the Legolas code.
This allows users to easily read datfiles into python and do analysis. See the **[website]** for more information.

## Requirements
Below is a list of requirements (including packages) needed in order to run Pylbo:
- Python 3.6+
- [Numpy](https://numpy.org), for obvious reasons.
- [Matplotlib](https://matplotlib.org), for plotting.
- [f90nml](https://f90nml.readthedocs.io/en/latest/), to handle reading and writing namelists.
- [tqdm](https://tqdm.github.io), used for progress bars.
- [psutil](https://psutil.readthedocs.io/en/latest/), for management of multiprocessing resources during parallel runs.

## Installation
Pylbo comes with a `setup.py` installation script, which automatically checks and includes the above listed packages. To install, simply do
```bash
cd post-processing
python setup.py develop
```
where the `develop` argument is included so you can easily pull new updates without having to reinstall again.

## Formatting

Pylbo uses [Black](https://github.com/psf/black) and [flake8](https://flake8.pycqa.org/en/latest/) for code formatting and style, which are both explicitly checked in the CI tests.
