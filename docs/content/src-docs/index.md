---
layout: default
title: Source code docs
nav_order: 7
last_modified_date: 2020-07-07
---

Since everyone loves a properly documented code, we added extensive docstrings
to the internal source codes of both Legolas and the post-processing framework Pylbo.
These are in turn used to automatically generate documentation.

However, since Legolas is written in Fortran and Pylbo is written in Python,
it is difficult to decently document both codes using the same tool. That's why
we opted to use language-specific tools: [FORD](https://github.com/Fortran-FOSS-Programmers/ford) for
the generated Legolas documentation, and [Sphinx](https://www.sphinx-doc.org/en/master/) for the 
generated Pylbo documentation.

Both of these have their own generated webpages, accessible through this page using the buttons below.


<span class="fs-6">
[Legolas docs](ford/index.html){: .btn .btn-green } [Pylbo docs](sphinx/index.html){: .btn .btn-green }
</span>



