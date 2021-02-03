---
title: The Python package Pylbo
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-02-02
---
When Legolas finishes solving the eigenvalue problem it writes all information to a specially
tailored datfile. First a header is written, containing information such as the geometry, chosen equilibrium,
parameters, unit normalisations and names of the equilibrium variables. The header is followed by the actual data
(eigenvalues, equilibrium arrays etc.), supplemented with eigenfunctions and matrix data if those are requested as well.
As all of this data is written out in a binary format for efficient data-storage, it is not so straightforward to
actually read its information since you'll have to keep track of which variables are in there and in which order
they are stored.

Since this is far from user-friendly we developed the post-processing Python package, `Pylbo`, short for
"**Py**thon for **L**egolas **B**inary **O**utput". Pylbo enables you to easily load in Legolas datfiles
and access all information stored within. We even developed special classes to do post-processing, for more information
and examples we refer to [this page](../using_pylbo). At some point in the future Pylbo will move to its own dedicated
repository and be included as a submodule in the legolas repository, but for now it is included in the
[`post-processing`](https://github.com/n-claes/legolas/tree/master/post_processing) folder.

Note that if you do Legolas runs on high resolution and save the eigenfunctions, the files can easily be a few
gigabytes in size. It may be useful to know that Pylbo _never_ loads the file into memory. Instead we keep track
of the various offsets of the datablocks, which are in turn used to read in data on a query-basis.
Doing it like this means a huge boost in performance and decreases memory usage significantly, and implies that you can
easily do analysis on larger-than-RAM datasets.

For example, say you did a huge resolution run where you saved the matrices and eigenfunctions, which resulted in a
10 Gb datfile. Loading this into Pylbo will be instantaneous, since the only thing that Pylbo "calculates" are the
data offsets. Whenever you query data, Pylbo will seek the corresponding offset in the binary stream
(which is a fast operation) and read only that specific datachunk into memory.
To put some numbers on it, you can load a series of datfiles each a few Gb in size simultaneously, and Pylbo will use
a few Mb of memory, tops.

## Code formatting
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

The entire Pylbo source code is formatted using [Black](https://github.com/psf/black), which automatically takes
care of trailing whitespaces, improper indents etc. We explicitly check the formatting during our
[automated tests](../../testing/test_pylbo/). Hence, if you want to make changes to the Pylbo source code,
make sure that everything is properly formatted or the automated tests will fail.
You can either run Black yourself after you have finished editing,
or you can do this automatically using [pre-commit](https://pre-commit.com) hooks. Up to you.

We also run [flake8](https://flake8.pycqa.org/en/latest/) as style enforcement.

More information on how these style and formatting tools are configured can be found [here](../../testing/test_pylbo/#running-the-style-checks).
