---
title: Source code documentation
layout: single
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-08-28
---

Since everyone loves a properly documented code (we sure do!), we added extensive docstrings
to the internal source codes of both Legolas and the post-processing framework Pylbo.
These are in turn used to automatically generate documentation.

However, since Legolas is written in Fortran and Pylbo is written in Python,
it is difficult to decently document both codes using the same tool. That's why
we opted to use language-specific tools: [FORD](https://github.com/Fortran-FOSS-Programmers/ford) for
the generated Legolas documentation, and [Sphinx](https://www.sphinx-doc.org/en/master/) for the 
generated Pylbo documentation.

Both of these have their own generated webpages, accessible through this page using the buttons below.
These are generated completely automatic and are hence kept independent of the main website.
We really want to keep these as up-to-date as possible, so if you notice anything missing please let us know,
preferably through an [issue on GitHub](https://github.com/n-claes/legolas/issues).
Also tag it with the yellow-green **#docs** label while you're at it, so we can easily keep track of them! 

[Legolas docs](../ford/index.html){: .btn .btn--primary .btn--large } 
[Pylbo docs](../sphinx/index.html){: .btn .btn--primary .btn--large }



