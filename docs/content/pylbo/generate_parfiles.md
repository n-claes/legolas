---
title: Generating parfiles
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-09-03
---

Pylbo comes with a method to generate parfiles for you, something that is extremely handy when doing multiruns
and you want to vary one or more parameters. All that is needed is to specify the namelist items in a dictionary,
and pass them to the [`generate_parfiles`](../../src-docs/sphinx/autoapi/pylbo/index.html#pylbo.generate_parfiles) method.
You don't have to specify which variable goes where, Pylbo will automatically take care of that and places them in their
corresponding namelist. Additionally, you can ask that filepaths and output folders are automatically resolved and placed
in the parfile as well, see below for an example.

**Tip:** when specifying paths using Python, the standard library package [`pathlib`](https://docs.python.org/3/library/pathlib.html)
is really recommended. This will resolve filepaths on any file system, and lets you easily work with relative paths as well.
That way you can specify the relative path to a file, copy the script over to another machine, run it, and the path strings
get automatically updated.
{: .notice--success}

# Generating a single parfile
To generate a single parfile, specify the setup as a dictionary. An example for Suydam cluster modes is given below:
```python
import pylbo

config = {
    "geometry": "cylindrical",
    "x_start": 0,
    "x_end": 1,
    "gridpoints": 101,
    "parameters": {
        "k2": 1.0,
        "k3": -1.2,
        "cte_rho0": 1,
        "cte_v02": 0,
        "cte_v03": 0.14,
        "cte_p0": 0.05,
        "p1": 0.1,
        "alpha": 2.0
    },
    "equilibrium_type": "suydam_cluster",
    "logging_level": 0,
    "show_results": False,
    "write_eigenfunctions": False,
    "write_matrices": False,
    "basename_datfile": "suydam_modes",
    "output_folder": "output"
}
parfile = pylbo.generate_parfiles(parfile_dict=config, 
                                  basename_parfile="suydam_parfile", 
                                  output_dir="output_parfile")
```
When the key `parameters` (this should be a dictionary itself) is found, Pylbo automatically sets `use_defaults` to 
`False`. The above code sample will create a parfile named `suydam_parfile.par` in the directory `output_parfile`
(relative to the directory where the script is run from). If the latter directory is not found it is automatically created.
The variable `parfile` here can be immediately used to run Legolas
(see [interfacing with Pylbo](../../getting-started/running/#2-interfacing-with-pylbo)), like this:
```python 
pylbo.run_legolas(parfile)
```
The actual parfile looks like this:
```fortran 
&equilibriumlist
    equilibrium_type = 'suydam_cluster'
    use_defaults = .false.
/

&gridlist
    geometry = 'cylindrical'
    gridpoints = 101
    x_end = 1
    x_start = 0
/

&paramlist
    alpha = 2.0
    cte_p0 = 0.05
    cte_rho0 = 1
    cte_v02 = 0
    cte_v03 = 0.14
    k2 = 1.0
    k3 = -1.2
    p1 = 0.1
/

&savelist
    basename_datfile = 'suydam_modes'
    logging_level = 0
    output_folder = 'output'
    show_results = .false.
    write_eigenfunctions = .false.
    write_matrices = .false.
/
```
which is a perfectly formatted Legolas namelist.

# Generating multiple parfiles
You'll probably be able to generate a single parfile faster by hand instead of using Pylbo. However, the method
described above excels in varying parameters, and modifying the parfiles accordingly. Say, you want to vary the parameter
`alpha` in the above example between 1 and 4, and do 50 runs. This means you have to create and edit 50 parfiles by hand, 
which is a long and painful chore. Instead, simply edit the above `config` dictionary with the following:
```python 
import numpy as np
...
"number_of_runs": 50,
"parameters": {
    "k2": 1.0,
    "k3": -1.2,
    "cte_rho0": 1,
    "cte_v02": 0,
    "cte_v03": 0.14,
    "cte_p0": 0.05,
    "p1": 0.1,
    "alpha": np.linspace(1, 4, 50)
},
...
```
where you have to add the `number_of_runs` entry. This is purely for a sanity check, and to make sure that everything
is nice and consistent. Pylbo will throw an appropriate error if there is an inconsistency between this number and the
number of runs you want.

Now you call `generate_parfiles` in exactly the same way, which will create 50 parfiles where a number `xxxx` is prepended
in front of the filename, the same holds true for the datfiles. 
You can even combine multiple variations. The example below illustrates this:

```python
import numpy as np
import pylbo

config = {
    "geometry": "cylindrical",
    "x_start": 0,
    "x_end": 1,
    "gridpoints": np.ones(50) * 100,
    "parameters": {
        "k2": 1.0,
        "k3": -1.2,
        "cte_rho0": 1,
        "cte_v02": 0,
        "cte_v03": np.linspace(0.05, 0.25, 100),
        "cte_p0": 0.05,
        "p1": 0.1,
        "alpha": np.linspace(1, 4, 100)
    },
    "number_of_runs": 100,
    "equilibrium_type": "suydam_cluster",
    "logging_level": 0,
    "show_results": False,
    "write_eigenfunctions": False,
    "write_matrices": False,
    "basename_datfile": "suydam_modes",
    "output_folder": "output"
}
config["gridpoints"][0:10] = np.ones(10) * 250
parfiles = pylbo.generate_parfiles(parfile_dict=config)
```
This generates 50 parfiles, in which the parameter `cte_v02` is varied together with the parameter `alpha`.
Additionally, all runs are done using 100 gridpoints, except that we redefined the first 10 runs with 250 gridpoints.

Generating parfiles using Pylbo is fully customisable and quite convenient for a large number of runs.
