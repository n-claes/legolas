---
title: The Legolas testing framework
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "laptop-code"
last_modified_at: 2021-07-27
---

This page gives a detailed overview on the various Legolas testing suites, how to run them and what they are testing.

As the Legolas codebase becomes larger when new features are implemented, it is important to test everything
thoroughly to make sure that new or modified code does not break backwards compatibility or introduces bugs
in pieces of code that previously worked.

Legolas has two levels of testing suites: unit tests and regression tests, all of which run on both the `master` and
`develop` branches for every commit or pull request. The Pylbo framework has its own dedicated testing suite.
All of these test suites are decouples from each other and test distinctly different things: the unit tests are low-level
and test individual subroutines and functions in the Legolas source code using [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).
The regression tests on the other hand are high-level and compare code output from a recent commit to previously known results, while
the Pylbo tests check data loading and management, interfacing with Legolas and visualisations.

The entire testing framework is supplemented with code coverage indications.

## Our testing policy
When developing new features or modifying existing code we have a few general guidelines that we follow:

{% capture guidelines %}
1. **All** functions/subroutines in the Legolas or Pylbo source codes should be accompanied by one or more unit tests.
2. Adding a new equilibrium means adding **at least one** regression test.
3. Code coverage should not decrease (although in some cases a minor decrease is unavoidable).
{% endcapture %}
<div class="notice--success">
  {{ guidelines | markdownify }}
</div>

## Style checks
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

All Python files in the legolas repository (so Pylbo, setup files and test files) adhere to [Black](https://github.com/psf/black)
formatting standards, which automatically take care of trailing whitespaces, improper indents etc.
We also run [flake8](https://flake8.pycqa.org/en/latest/) as style enforcement.

We explicitly check for both Black and flake8 formatting and style during the automated tests and strictly enforce this.
Many IDE's support Black/flake8 plugins which can automatically format the file on save, so it may be useful to set this up
if you are making edits to the source code.
## The Legolas unit tests
### Running the tests
To run the unit tests you will need to have pFUnit installed, which has some
[prerequisites](https://github.com/Goddard-Fortran-Ecosystem/pFUnit#prerequisites). Make sure Legolas is compiled, then navigate to
the `tests` folder and then down one level to the unit tests folder. In there, compile the tests and run:
```bash
mkdir build
cd build
cmake ..
make
cd ..
./test_legolas
```
This will execute the entire unit test suite.

### Creating a new unit test
All similar tests are grouped together in dedicated `.pf` files. Either add a new test to an already existing file, or create
a new file. The module `mod_suite_utils.f90` has a few convenient routines which can be used to set up or tear down new tests.
As a general rule of thumb all arrays should be deallocated and all global variables reset when a test finishes (on either failure or success),
to ensure that every test can run independently of the others.
If you added a new `.pf` file you will have to add it to the list in the `CMakeLists.txt` file in the unit tests folder.

Note that unit tests should test only _a single_ (part of a) subroutine/function and should finish quickly. If you want to test multiple things,
just split them in multiple tests as much as possible.

## The Legolas regression tests
The regression tests run the various pre-implemented equilibria and compare the results with previously stored datfiles. Contrary to the Legolas unit tests
this suite is Python-based and uses [`pytest`](https://docs.pytest.org/en/latest/) together with the Pylbo-Legolas interface to set up and run everything.
These tests mainly compare eigenvalues and eigenfunctions between new and previously stored results using an image-based comparison and the root-mean-squared difference.

### How the regression tests work
During the regression tests we check both the eigenvalue spectrum and (only in some cases) the eigenfunctions.
#### Spectrum tests
These are done for _all_ pre-implemented equilibria, and every spectrum test typically consists of four steps:
1. Generate the datfile for the corresponding equilibrium
2. For both the "new" datfile and the stored datfile:
  - load the file
  - for every region in the complex plane generate a spectrum (usually 3-4 images)
  - save the images to a temporary folder
3. Do a root-mean-squared (RMS) image comparison to check spectral similarity, fail the test if
  the RMS is above the tolerance (which is 2 by default).
4. Save the difference between the two images as a new image on failure, otherwise delete
  the images.

An example of a failing test is shown here, where the where the background flow parameter for the Suydam cluster modes
was slightly changed:
![image-center]({% link /assets/images/example_test_fail.png %}){: .align-center}
The left image is the expected baseline result, the center image the test result and the right image the difference
between both. In all these images there is no background alpha present, so all deviations in the RMS is purely due to shifting eigenvalues.

We also test multispectra, which is a composed spectrum for the same equilibrium but run with varying parameters. The testing principle remains
the same: once the composed spectrum is created the baseline and test result are again compared on an RMS-based criterion.

#### Eigenfunction tests
These tests use a similar background as the spectrum tests and are also RMS-image-based. We do not check the eigenfunctions for all pre-implemented equilibria,
only those that have clear, isolated eigenvalues and associated eigenfunctions that are properly resolved at low resolution.
When the test datfile is generated some pre-set eigenvalues are selected and the eigenfunctions are queried. Then, for each eigenvalue we plot the sum of the real
and imaginary parts of each eigenfunction as a function of the grid and save the images on a 3x3 grid (since we have 8 panels per image, per eigenvalue).
The same principle as above is then applied, where the RMS value of the curves is checked against each other.

### Running the tests
To run the regression tests make sure you have pytest installed, along with Pylbo and all its prerequisites. Then, navigate to the `tests` folder and move down
one level to the regression test directory. There you execute pytest as follows
```bash
pytest -v regression.py test*
```
It is important that the `regression.py` file comes first, since this one is responsible for generating the images and properly setting up everything.
The verbose flag `-v` is not necessary but is recommended to better track the progress.

<i class="fas fa-lightbulb" aria-hidden="true"></i>
**Note**: tests that are passing will clean up after themselves by deleting the generated datfile and images.
You can ask the test suite to keep the images instead (for visual inspections, for example) by simply adding the `--keep-files` flag when running the tests.
{: .notice--info}


### Adding new tests
To add a new regression test create a new file in the regression test folder, and make sure to prepend the name with `test_` so it can be picked up
by pytest. There you create a dictionary similar to when you generate a parfile with Pylbo which is used to configure and set up the test, take a look
at the already present files for more details.

For example, create a file called `test_example.py` with the following:
```python
my_test_setup = {
  "name": "my_test",
  "config": {
      "gridpoints": 51,
      "equilibrium_type": "adiabatic_homo",
      "parameters": {
          "k2": 0,
          "k3": np.pi,
          "cte_rho0": 1.0,
          "cte_T0": 1.0,
          "cte_B02": 0.0,
          "cte_B03": 1.0,
      },
      "logging_level": 0,
      "show_results": False,
      "write_eigenfunctions": True,
  },
  "all_eigenvalues_real": True,
  "image_limits": [
      {"xlims": (-600, 600), "ylims": (-0.05, 0.05), "RMS_TOLERANCE": 2.5},
      {"xlims": (-50, 50), "ylims": (-0.05, 0.05)},
  ],
  "eigenfunctions": [
      {"eigenvalue": 2.67131},
      {"eigenvalue": 2.54724},
  ],
}
```
Next you open the file `regression.py`, import `my_test_setup` and add it to the appropriate list of tests to run.
All the rest will be taken care of automatically.

The tests are configured in such a way that everything is controlled through this dictionary. The `"name"` key is just the name
of the tests, the `"config"` key represents the usual dictionary to generate a parfile.
The key `"image_limits"` should be a list, with every item a dictionary itself with the following keys:
- `"xlims"`: this is a tuple representing the range in Re$(\omega)$ of the spectrum
- `"ylims"`: this is a tuple representing the range in Im$(\omega)$ of the spectrum
- `"RMS_TOLERANCE"`: [_optional_], this is a float representing a custom RMS tolerance to be used for this specific test. If not given
  the default value is used.

The key `"eigenfunctions"` is similar: this is also a list with the `"eigenvalue"` key, representing the (approximate) eigenvalue which eigenfunctions
you want to test. We again have two tests here, one per eigenvalue. These list items also accept the optional `"RMS_TOLERANCE"` key.

Finally we have the _optional_ `"all_eigenvalues_real"` flag, which is a special test for when you are sure that the spectrum should only yield real eigenvalues.
If this key is present and its value is `True`, an extra test will run which will fail if any one of the eigenvalues has a non-zero imaginary part.

The above example hence corresponds to five tests: one to check that all eigenvalues are real, two spectrum tests and two eigenfunction tests.
Multispectrum tests are similar, take a look at the files with `test_multi` prepended to their name.

## The Pylbo unit tests
To run the Pylbo unit tests you also need pytest, along with the [`timeout`](https://pypi.org/project/pytest-timeout/) plugin which you can install as
```bash
pip install pytest-timeout
```
Then, ro tun the tests navigate to the `tests` folder, down one level to the Pylbo unit tests and run
```bash
pytest -v
```