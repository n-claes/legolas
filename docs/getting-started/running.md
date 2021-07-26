---
title: Running Legolas
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: false
last_modified_at: 2021-07-26
---

Once you cloned the repository and installed the Pylbo framework you are set to run your first problem.
As an example here we will run one of the pre-implemented equilibria already present in the source code.

## Running your first problem
Navigate to a directory of your choosing where you want to place the configuration files and output, and call the setup script:
```bash
cd my_favourite_directory
setuplegolas.py
```
Follow the instructions, this script will do and ask you a couple of things:
1. It copies over a `CMakeLists.txt` file if not already present
2. Ask if you want to copy over `pylbo_wrapper.py`, you need this if you want to plot your results immediately after running. **Answer yes**
3. Ask if you want to copy over a parfile if none was found. **Answer no**
4. Ask if you want to copy over a default submodule template `smod_user_defined.f08` if none was found. **Answer no**

Since we are running a pre-implemented equilibrium we do not need the `smod_user_defined.f08` file, this one is only needed if you want to
set up and run your own problem.

Next call the build script:
```bash
buildlegolas.sh
```
If everything went well you should have 5 items in your current working directory: a `CMakeLists.txt` file, a `build` folder, an `output` folder,
the file `pylbo_wrapper.py` and the `legolas` executable.

<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
It is NOT recommended to create custom build or run directories _inside_ the main legolas folder.
Keeping the repository clean ensures that you don't run into merge conflicts when updating the code.
Run the setup scripts in a separate directory and let them take care of linking.
{: .notice--warning}

## Creating the configuration file
Configuring Legolas is done through use of a "parfile", where various options are passed in the form of Fortran namelists. For an overview
of different possibilities you can take a look [here](../../general/parameter_file).
During setup we answered "no" when the script asked to copy over the parfile.
 For now, create a new file called `my_setup.par` in which you copy the following lines:
```fortran
&gridlist
  gridpoints = 50
/

&equilibriumlist
  equilibrium_type = "rayleigh_taylor"
/

&savelist
  write_eigenfunctions = .true.
  show_results = .true.
/
```

## Running the code
Finally call `legolas` and supply the parfile you just created as argument:
```bash
./legolas -i my_setup.par
```
This will run a setup which defines Rayleigh-Taylor instabilities in a Cartesian geometry,
the spectrum itself should correspond to case **a** in Section 13.2, p. 487 of [MHD of Laboratory and Astrophysical Plasmas](http://doi.org/10.1017/9781316403679).

When the run is completed (this should take a few seconds) the code will immediately fire up the post-processing framework since we set `show_results = .true.` in the parfile.
A few interactive figures will pop up, including the spectrum. Click on one of the eigenvalues in the figure and press `Enter`, this will draw the corresponding eigenfunctions.
You can cycle through multiple eigenfunctions using the arrow keys, right-clicking a point unselects it and pressing the `d` key clears all selected points.

The legend shows the various continuum regions, which can also be interactively toggled by clicking on their legend entries.
More information on the interface can be found [here](../../pylbo/using_pylbo/#interactive-continua--eigenfunctions).

This covers the basics of running Legolas on a single equilibrium. You can take a look at the various [pre-implemented equilibria](../../general/equilibria) and
modify the `equilibriumlist` or `paramlist` accordingly to play around with different options.
