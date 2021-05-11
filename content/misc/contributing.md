---
title: Contributing
layout: single
sidebar:
  nav: "leftcontents"
toc: true
last_modified_at: 2020-10-27
---

Legolas is developed as an open-source code, which means that we welcome contributions that improve the code in any
way. The ideal way to do this is following a proper
[GitHub workflow](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
This means that you first fork the main Legolas repository, make the edits you want and then submit a
[pull request](https://github.com/n-claes/legolas/pulls).
Since everything is contained in the main repository you can edit the Legolas source code, the docs or the Pylbo
source code this way.

Once the pull request (PR) is submitted it will be up for review by one of the main developers. Please make use of the dedicated
[labels](https://github.com/n-claes/legolas/labels) to label your pull request (or issue), this allows us to easily keep track
of the various types of additions that are being made.

**Note:** all PR's should be directed against the [`develop`](https://github.com/n-claes/legolas/tree/develop) branch,
and NOT against the [`master`](https://github.com/n-claes/legolas/tree/master) branch.
We keep the `master` branch dedicated to stable releases, and keep "new" features on the `develop` branch until we
release a new major/minor version. Doc updates are an exception to this, since the website is built from the `master` branch.
{: .notice--info}

# Adding to the Legolas code
For contributions to the Legolas source code we have the following guidelines:
- Variable names should be self-explanatory but consistent (with the exception of parameters)
- When printing information to the console, make use of the dedicated [logging module](../../src-docs/ford/module/mod_logging.html),
in particular the [`log_message()`](../../src-docs/ford/proc/log_message.html) subroutine.
- New features should go in dedicated (sub)modules if possible.
- Consistent indentation across the entire code and free-format syntax.
- In the case of new subroutines: add a test to the [core tests](../../testing/test_core)
- In the case of new equilibria: add a test to the [regression tests](../../testing/test_regression)

<i class="fas fa-exclamation-triangle"></i> **Important:**  
Every new bit of code _has_ to be supplemented by corresponding documentation, which should be clear
and concise. For Legolas docstrings we use [FORD syntax](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation),
take a look at the already documented modules for more information and examples.
{: .notice--danger}

# Adding to the Pylbo code
For contributions to the Pylbo source code there are no "real" guidelines, besides the regular Python guidelines.
This means that all new or modified code should follow the [PEP8 guidelines](https://www.python.org/dev/peps/pep-0008/) (roughly).
An exception to this is the maximum line length, since Black uses 88 (instead of the 79 suggested by PEP).
Method names are also quite flexible, if they are concise and clear it's fine.

Additionally, whenever you add code you should also add a dedicated unit test to the [pylbo tests](../../testing/test_pylbo).

<i class="fas fa-exclamation-triangle"></i> **Important:**  
Every new bit of code _has_ to be supplemented by corresponding documentation, which should be clear
and concise. For Pylbo docstrings we use [Sphinx syntax](https://www.sphinx-doc.org/en/master/usage/quickstart.html#documenting-objects),
supplemented by the [NumPy docstring convention](https://numpydoc.readthedocs.io/en/latest/format.html).
Take a look at the already documented methods and classes for more information and examples
{: .notice--danger}

# Adding to the documentation
## Edits to the website
The documentation (this website) is built from the [`docs`](https://github.com/n-claes/legolas/tree/master/docs)
folder on the `master` branch, using [GitHub Pages](https://pages.github.com).
We make use of [Jekyll](https://jekyllrb.com) and the [minimal-mistakes](https://github.com/mmistakes/minimal-mistakes) theme.
For instructions and detailed guides we refer back to those pages.

When making edits to the documentation it is a good idea to compile it locally before pushing to the repository,
that way you avoid small mistakes like broken links and incorrect formatting. You need a Ruby development environment
for this, instructions are given on the [Jekyll website](https://jekyllrb.com/docs/installation/). We'll sum them up below:
- Ruby version 2.5.0 or higher, you can check through `ruby -v`
- RubyGems, you can check through `gem -v`
- GCC and Make, you can check through `gcc -v`, `g++ -v` and `make -v`.

Once you have this set up, you install the `Jekyll` and `Bundler` Gems
```bash
gem install jekyll bundler
```
and that's it! Now you simply navigate to the `docs` folder and type
```bash
bundle exec jekyll serve
```
after which the website builds in a folder `_site`, and it will tell you the local URL that you can paste in your
browser of choice (denoted by "_Server address_). Due to the `serve` statement you can make additional edits
to the documentation and refresh the webpage to immediately look at the result.

## Edits to the Pylbo docs
For the pylbo docs you'll need [Sphinx](https://www.sphinx-doc.org/en/master/) configured,
together with its dependencies. You can install all of them through
the `install_dependencies.sh` script in the `docs` folder:
```bash
sh install_dependencies.sh
```
or, if you want to do it manually:
```
pip install -U sphinx
pip install numpydoc
pip install sphinx-rtd-theme
pip install lazy-object-proxy==1.4
pip install sphinx-autoapi
```
Once this is done, you can execute
```bash
sh build_sphinx.sh
```
to build the Sphinx documentation, which is configured in such a way that Jekyll finds the generated html files.

## Edits to the Legolas docs
For the Legolas docs you'll need [Ford](https://github.com/Fortran-FOSS-Programmers/ford).
The `install_dependencies.sh` script installs it together with the Sphinx
dependencies, but just in case you want to do it manually:
```bash
pip install ford
```
Once this is done, you can execute
```bash
ford build_ford.md
```
to build the Ford documentation, which is configured in such a way that Jekyll finds the generated html files.
