---
layout: single
sidebar:
  nav: "leftcontents"
toc: true
toc_label: "Installation Guide"
toc_icon: "cogs"
last_modified_at: 2020-11-27
---

Due to the heavy use of object-oriented features (Fortran 2003) and submodules (Fortran 2008), Legolas
requires relatively recent Fortran compilers. This page gives a detailed overview of required
and optional dependencies to successfully build and run both Legolas and the post-processing framework
[Pylbo](../../pylbo/installing_pylbo/). Note that the compilers noted below are only recommended,
and it is quite possible that Legolas also builds with lower versions (we haven't tested that).

# Dependencies
## Compilation
- Fortran compiler
    - gfortran v8.x+ (recommended and tested)
    - Intel compilers v18.0+ (not tested)
- CMake v3.12+
- Python v3.6+
- make

**Note:** Python is needed for Pylbo, not for Legolas itself. You can still run Legolas
if the Python requirements are not satisfied, however you will not be able to immediately
see the results after the run finishes (so set `show_results=.false.` in the parfile).
{: .notice--info}

## Linear algebra
### BLAS and LAPACK
The [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/)
linear algebra packages are **required** dependencies, and you will
not be able to compile without them. We recommend version 3.5 or higher.
CMake is configured in such a way that both libraries should be found and linked automatically if
they are installed. In general, it is sufficient if you have them installed through
```bash
sudo apt-get install libblas-dev
sudo apt-get install liblapack-dev
```
if you're on Linux, or through
```bash
brew install openblas
brew install lapack
```
using [HomeBrew](https://brew.sh) on macOS. Note that macOS ships with a default BLAS/LAPACK installation
as part of the [vecLib](https://developer.apple.com/documentation/accelerate/veclib) framework.
If you did a manual compilation of the BLAS and LAPACK libraries (if you don't have sudo rights, for example),
CMake may not find the libraries by default. In that case it will throw a warning, and you may have to set the
`$BLAS_LIBRARIES` and `$LAPACK_LIBRARIES` variables which link to the compiled libraries.


### ARPACK
The [ARPACK](https://www.caam.rice.edu/software/ARPACK/) library (and its parallel counterpart
PARPACK is an **optional** dependency, so
Legolas will compile and run just fine if you don't have this installed (related modules are
conditionally compiled). Also here CMake will try to automatically find and link the libraries if installed.
We recommend using maintained [arpack-ng](https://github.com/opencollab/arpack-ng) repository to
install ARPACK. If CMake fails to find the library and you have it installed, you can set the
`$ARPACK_ROOT` environment variable, pointing to the top-level of your local `arpack-ng` repository.


### SCALAPACK
The [SCALAPACK](http://www.netlib.org/scalapack/) package is an **optional** depencency.
CMake has been configured to find and link this library, but this feature is currently
disabled since Legolas does not have SCALAPACK routines implemented (yet).


### MUMPS
The [MUMPS](http://mumps.enseeiht.fr) library is an **optional** dependency.
CMake has been configured to find and link this library, but this feature is currently
disabled since Legolas does not have MUMPS routines implemented (yet).


# Building Legolas
## Obtaining the code
Legolas can be obtained by cloning the online repository:
```bash
git clone https://github.com/n-claes/legolas.git
```
which will put it in the local repository `legolas`.
CMake has been configured for an out-of-source build, meaning that the repository stays clean.
This allows for easy updates through `git pull`.

## Environment variables
For an easy setup we recommend setting the environment variable `$LEGOLASDIR` and adding
the `setup_tools` folder to your PATH. This can be done by editing your `.bashrc` (or `.zshrc` on macOS)
as follows
```bash
export LEGOLASDIR='path_to_the_legolas_directory'
PATH="$PATH:$LEGOLASDIR/setup_tools"
```
The last line allows for easy access to the `buildlegolas.sh` and `setuplegolas.py` scripts in the `setup_tools`
folder, such that they can be called from any directory.

## Compiling the code
Compiling the code is quite straightforward, and for your convenience we have provided a simple shell script
to do everything at once. Compiling the code can either be done from inside the repository or from a dedicated folder somewhere.
The latter option will be particularly useful when you're setting up your own problems.

**Note:** whichever of the two options you choose, `Legolas` is _always_ compiled in a directory called `build`
_inside_ the main repository, which is ignored by git. Because we make heavy use of submodules (which prevent
compilation cascades when changes are made), you don't have to recompile the code every time you set up a new problem.
Only the modified user-defined submodule is recompiled, and the "new" executable is placed in the same directory.
{: .notice--info}

### 1. In-repository build
The first option is building inside the repository. To do so, navigate to the legolas source directory and do
```bash
sh setup_tools/buildlegolas.sh
```
This will create a directory `build` inside the repository (which is ignored by git), and places the `legolas`
executable in the topmost directory. You can also do it manually:
```bash
mkdir build
cd build
cmake ..
make
```
An in-repository build is fine if you want to run pre-implemented problems, but note that when you start modifying
equilibria (or adding your own) you are essentially modifying the source files. This means that you will probably run
into merge conflicts soon when you're updating the code (and we all hate those!) so we don't recommend doing it like this.

### 2. Out-of-repository build (recommended)
This is the recommended way to build Legolas. This requires you to have set the environment variable and PATH
modification described above. To setup a folder for a Legolas build, call `setuplegolas.py` from the command line and
follow the instructions. The script will do (and ask you) a couple of things:
1. Copy over `CMakeLists.txt`.
2. Copy over `pylbo_wrapper.py`, which you will need if you want to plot the results afterwards.
3. Copy over a default parfile if none was found.
4. Copy over the default user submodule template if none was found.
5. Modify the build script to include the custom submodule if found.

Next you simply build as described before, either through `buildlegolas.sh` which does everything automatically, or
manually through
```bash
mkdir build
cd build
cmake ..
make
```
In both cases the executable is placed in the same directory as the parfile and user submodule.

# Doing a clean build
Since we use CMake for compilation there is no exact `make clean` equivalent as there is with GNU Make.
Instead, you can do one of the following things:
1. Navigate to the `build` directory (inside the repository or the local one) and do `make clean`.
   This removes the compiled object files and lets you do a fresh `make` compilation.
2. Remove the (local) `build` directory, re-create it and compile again. This is perhaps best way to 
   do it, since that means that the `CMakeCache` is removed as well.
   
You can also do this automatically from any directory by supplying an additional argument to `buildlegolas.sh`.
Say you just compiled in a local directory and you want to do a fresh compilation, simply call
```bash
buildlegolas.sh clean
```
This will remove the `build` folder in the main legolas directory, the local `build` directory and the local executable (if one is found).
