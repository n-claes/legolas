---
title: Testing the Pylbo framework
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-10-27
---

The Pylbo framework is also tested, we use [pytest](https://docs.pytest.org/en/stable/) for this.
To run these locally you'll need pytest installed, which is explained in their
[installation instructions](https://docs.pytest.org/en/stable/getting-started.html#install-pytest).
You will also need additional pytest plugins, more specifically [pytest-timeout](https://pypi.org/project/pytest-timeout/)
which is used to fail tests if they take too long, and [pytest-mpl](https://pypi.org/project/pytest-mpl/) which
is used for image comparisons to a given baseline.
You can install all of these through
```bash
pip install pytest pytest-mpl pytest-timeout
```

As mentioned in the [code formatting](../../pylbo/about_pylbo/#code-formatting) section we use
Black as a formatting guide and flake8 as style guide.
These are both **strictly** enforced, meaning that the automated tests for
Pylbo will fail if Black or flake8 exit with a non-zero integer.

## Running the style checks
Installing Black and flake8 can be done through
```bash
pip install black flake8
```
To check if your modifications satisfy the `flake8` style you can do
```bash
flake8 --count --exclude=__init__.py --extend-ignore=E203,W503 --max-line-length=88 --show-source --statistics your_file_or_folder
```
which are the flags that are used in the automated tests.
Note that we use a line length of 88 (in contrast to the default 79 that PEP8 uses) for more flexibility.
This is also the default setting for Black.
After `flake8` you can check your `black` formatting through
```bash
black --diff --check your_file_or_folder
```
Currently all files in these folders are checked:
- `post_processing`
- `tests/regression_tests`
- `tests/pylbo_tests`

**Tip**: You can automatically format your files by removing the options for black:
`black your_file_or_folder`
{: .notice--success}

## Running the pylbo tests
Once all style checks are satisfied, you navigate to `tests/pylbo_tests` and execute
```bash
pytest --mpl --mpl-results-path=results --verbose
```
This will automatically search for all unit tests and execute them, and you'll see a detailed overview during runtime.

These tests are run automatically for each commit and pull request to the `master` and `develop` branches.
