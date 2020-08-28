---
title: Running Legolas
layout: single
sidebar:
  nav: "leftcontents"
toc: true
---

Once you have compiled Legolas it will have placed an executable in your current working directory (or in the
source directory if you chose an [in-repository build](../installation/#1-in-repository-build)).
The executable can be run in two ways, either through direct calling or through the Pylbo interface.

# 1. Running Legolas directly
Passing parameters to Legolas is done through the use of a parameter file.
If you've set up your current running directory using `setuplegolas.py`, the script will have asked to copy over the
parfile. Edit this file to your needs, take a look at the different [possibilities](../../general/parameter_file).
Afterwards you simply call the executable from the command line, passing the parfile as an option:
```bash
./legolas -i parfile.par
```
The output is one single file, placed by default in a subdirectory `output` which is created during compilation.

# 2. Interfacing with Pylbo
A second way to run Legolas is through the Pylbo-interface, which is is a bit more cumbersome than running directly
but extremely handy when you are doing parametric studies or want to run multiple setups at once. This obviously
requires you to have Pylbo installed, explained [here](../../pylbo/installing_pylbo).
To start, fire up Python (or create a script) with the following lines
```python
import pylbo

parfiles = ["path_to_file1", "path_to_file2"]
pylbo.run_legolas(parfiles, remove_parfiles=False, nb_cpus=8)
```
the variable `parfiles` contains the paths to your parfiles as Strings or 
[pathlike](https://docs.python.org/3/library/pathlib.html) objects. The optional `remove_parfiles` keyword argument 
(`False` by default) removes the parfiles after the runs are completed, 
in case you want to automatically clean up afterwards. 
The `nb_cpus` kwargs is also optional, with a default value of 1. If you supply multiple parfiles and set this
number larger than one Pylbo will use Python's [multiprocessing](https://docs.python.org/3/library/multiprocessing.html)
module to parallelise the runs. A progressbar will be printed so you can keep track of the progress.

**Advice**: Please note that Legolas loads the _full_ matrices into memory, meaning that if you run multiple large
setups in parallel you will probably run out of memory and your system will start using swap memory to keep up, 
_greatly_ slowing down the calculations. 
A rule of thumb is running a maximum (combined) value of 3000 _base_ gridpoints at once, 
meaning 30 runs of 100 gridpoints each (`nb_cpus=30`), or 15 runs of 200 gridpoints each (`nb_cpus=15`), etc.
This depends on your available computational resources ofcourse, but if you start noticing slowdown you're
probably using too much memory.
{: .notice--warning}
