---
title: Code coverage
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-10-27
---

Since we would like to know how much of the code base is actually checked when running the tests, we've included
code coverage options to the build process. If this flag is enabled during building, running the tests will generate
additional files in the `build` directory which are processed using coverage tools to generate detailed reports.
You'll have to compile using gfortran, since we use [gcovr](https://gcovr.com/en/stable/) for code coverage.

To locally generate the report you have to remove the `build` folder and executable in the main Legolas directory.
Next you recreate the `build` directory and compile with coverage flags enabled:
```bash
rm -rf build
mkdir build
cd build
cmake -DCoverage=ON ..
make
```
You then navigate to the `tests/core_tests` folder and do the exact same thing, this ensures that the coverage
tools are linked to the `test_legolas` executable during the compilation process.

Next you simply run all tests and the reports will be automatically generated.
To process the results, navigate to the main legolas directory and do
```bash
mkdir coverage
cd build
gcovr . -r ../src --print-summary --html --html-details --output ../coverage/coverage.html
```
This will generate the report in the `coverage` folder, after which you can open `coverage.html` in your preferred browser.
In some cases it's possible that you have to specify the `gcov` executable in order for this to work, you can do that
through the additional `--gcov-executable` flag. More info [here](https://gcovr.com/en/stable/guide.html#the-gcovr-command).
