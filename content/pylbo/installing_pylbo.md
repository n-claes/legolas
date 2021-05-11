---
title: Installing Pylbo
layout: single
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-09-03
---
Pylbo comes as a ready-to-install Python package, contained in the `post-processing` folder of the legolas repository.
It depends on a few standard packages, listed below:

**Requirements**
- [Numpy](https://numpy.org), for obvious reasons.
- [Matplotlib](https://matplotlib.org), for plotting.
- [f90nml](https://f90nml.readthedocs.io/en/latest/), to handle reading and writing namelists.
- [tqdm](https://tqdm.github.io), used for progress bars.
- [psutil](https://psutil.readthedocs.io/en/latest/), for management of multiprocessing resources during parallel runs.

## Installing as a package (recommended)
This is the recommended installation option. To do so, navigate to the `post-processing` folder and execute
the `setup.py` script. This will automatically install the above listed packages as well if they are not already there.
Activate any conda/virtualenv environment beforehand in order to install it there.
```bash
cd post-processing
python setup.py develop
```
All of the dependencies listed above will be automatically handled if you run this.
The `develop` argument means that you can pull updates to the package from the upstream repository straight away,
without having to install it again.

## Sourcing the folder
Another possibility is adding the folder to your PYTHONPATH. If the above option does not work, 
i.e. if you can not manage your python environment for some reason, this is a valid alternative.


**macOS/Linux**  
Edit `.zshrc`/`.bashrc` or similar, we assume that you have already set the `$LEGOLASDIR` environment variable 
(see [installation](../../getting-started/installation#environment-variables)).
```bash
export PYTHONPATH=$LEGOLASDIR/post-processing:$PYTHONPATH
```
**Windows**  
Go to `Start` > `Control panel` > `Edit the system environment variables` > tab `"advanced"` > `environment variables`.
Create or edit the PYTHONPATH environment variable and add the _full path_ to the `post-processing` folder.
Multiple paths are separated by a semicolon `";"`.
