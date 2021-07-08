---
title: Running Legolas
layout: single
sidebar:
  nav: "leftcontents"
toc: true
last_modified_at: 2021-02-02
---

Once you have compiled Legolas it will have placed an executable in your current working directory (or in the
source directory if you chose an [in-repository build](../installation/#1-in-repository-build)).
The executable can be run in two ways, either through direct calling or through the Pylbo interface.

# 1. Running Legolas directly
Passing parameters to Legolas is done through the use of a parameter file.
If you've set up your current running directory using `setuplegolas.py`, the script will have asked to copy over the
parfile. Edit this file to your needs, take a look at the different [possibilities](../../general/parameter_file).
Afterwards you simply call the executable from the command line, passing the parfile as an option:
```bash
./legolas -i parfile.par
```
The output is one single file, placed by default in a subdirectory `output` which is created during compilation.

# 2. Interfacing with Pylbo
A second way to run Legolas is through the Pylbo-interface, which is extremely handy when you are doing parametric studies 
or want to run multiple setups at once. This obviously requires you to have Pylbo installed, explained [here](../../pylbo/installing_pylbo).
More information on how to run Legolas like this can be found on the [Interfacing with Legolas](../../pylbo/legolas_interface) page.

