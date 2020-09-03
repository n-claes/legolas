---
title: Using Pylbo
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_label: "Pylbo examples"
toc_icon: "laptop-code"
last_modified_at: 2020-09-03
---

Using Pylbo is quite straightforward, for a detailed guide on the API we refer to the 
[Pylbo documentation](../../src-docs/sphinx/autoapi/pylbo/index.html). This page contains a 
few examples on how to use the package. In what follows we assume that Pylbo has been installed 
(see [installing Pylbo](../installing_pylbo)) and that
```python
import pylbo
```
is present and does not return an error.

## Loading datfiles
You can choose to load a single datfile or a series of datfiles. Simply do
```python 
ds = pylbo.load("path_to_datfile")
```
to load a single datfile, or
```python 
files = ["path1", "path2", "path3"]
datasets = pylbo.load(files)
```
to load multiple files at once. The method `load` returns a 
[`LegolasDataContainer`](../../src-docs/sphinx/autoapi/pylbo/index.html#pylbo.LegolasDataContainer)
instance (or in the case of `datasets`, every element in the array is an instance),
which has a lot of helpful attributes and methods to calculate anything you may need.
Additional things are added on an as-needed basis.

## Analysing a single file
### Plotting the spectrum
Say you have loaded a single datfile, whose information is contained in the `ds` instance as above, and want to
take a look at the spectrum itself. 
Create a [`SingleSpectrum`](../../src-docs/sphinx/autoapi/pylbo/visualisations/spectrum/index.html#pylbo.SingleSpectrum) 
instance using `ds`, and query for a spectrum plot:
```python 
spec = pylbo.SingleSpectrum(ds)
spec.plot_spectrum(annotate_continua=True)
pylbo.show()
```
where we call `pylbo.show()` afterwards, this is a simple call to `matplotlib.pyplot.show()` so you don't have
to explicitly import the matplotlib package. This will immediately create a plotting window, and since
we set the `annotate_continua` kwarg to `True` Pylbo enables an interactive legend. Clicking on the legend items will 
enable/disable an overlay of the different continuum regions for your equilibrium, which can be either bands
(if genuine continua are present) or single points (if the continua are collapsed). Look at the 
[documentation](../../src-docs/sphinx/autoapi/pylbo/index.html#pylbo.SingleSpectrum.plot_spectrum) to see
what this method returns and what the possible kwargs are.

### Plotting the equilibrium/continua
You can also take a look at the equilibrium itself, useful for when you're implementing one yourself and doing
dry runs to iterate over it.
```python 
spec = pylbo.SingleSpectrum(ds)
spec.plot_equilibria()
spec.plot_continua()
pylbo.show()
```
which will create two figures, one contains the equilibrium and the other contains the slow, Alfvén and thermal continua,
all plotted versus the grid. You can pass a matplotlib figure and axis instance to this method as well, in that case
those will be used for plotting instead of creating a new figure.
See the [documentation](../../src-docs/sphinx/autoapi/pylbo/index.html#pylbo.SingleSpectrum.plot_continua).

### Interactively plotting eigenfunctions
Pylbo is equipped with a completely interactive eigenfunction plotter. This means that you can literally select a point
on the spectrum and plot the corresponding eigenfunctions. To do so:
```python 
spec = pylbo.SingleSpectrum(ds)
spec.plot_eigenfunction(merge_figs=True)
pylbo.show()
```
If the `merge_figs` kwarg is set to `True`, one figure with two subplots will be created. One figure (left) will contain
the spectrum, the right one contains the eigenfunction window. If it's `False`, two separate figures will be created.
When you select a point a red cross will be drawn on top of it, deselecting the point will remove the cross.
Below is an overview of the interactive controls:

**Note:** for ease of use the interactive controls are disabled when the matplotlib toolbar is enabled. This means that
you have to disable e.g. the zoom tool before you can select points. This is done to prevent points from being selected
during zooming, for example, which is not desired behaviour.
{: .notice--info}

| Key / mouse    | Description     |
| :---:         | :----          |
| _left click_ | Select a spectrum point. Is picker-based, so if the points are too closely clustered you'll have to zoom in a bit. You can select multiple points at once. |
| _right click_ | Remove point from selection. |
| _enter_ | Plots the eigenfunctions of the currently selected points. |
| _i_ | Swaps between the real and imaginary parts of the eigenfunctions. |
| _up arrow_ | Cycles to the next variable in the list. |
| _down arrow_ | Cycles to the previous variable in the list. |
| _d_ | Clears current selection. |
| _x_ | Closes figure(s). |
| _t_ | Retransforms the eigenfunctions using the scale factor, hence has no effect if `geometry = "Cartesian"`. You can use this to cycle between the $rv_r$ and $v_r$ eigenfunctions in cylindrical geometry, for example. |
| _w_ | Prints out the eigenvalues and indices of the selected points. Useful in combination with [`ds.get_nearest_eigenvalues`](../../src-docs/sphinx/autoapi/pylbo/index.html#pylbo.LegolasDataContainer.get_nearest_eigenvalues). |

### Manually getting eigenfunctions
If you're creating a figure yourself, you may want to retrieve one or more eigenfunctions corresponding to specific
eigenvalues. Assume you looked at the spectrum, and you found interesting eigenvalues at approximately `1 + 3j` 
and `-2 + 8j`. You can ask Pylbo to use those as "guesses", and find the "actual" eigenvalues in the dataset.
Another possibility is selecting them during interactive eigenvalue plotting and pressing the `w`-key to write the
actual values to the console. In general, your workflow looks like this:
```python
import pylbo 

ds = pylbo.load("datfile.dat")
# locate the indices and nearest eigenvalues to your guess
idxs, evs = ds.get_nearest_eigenvalues([1 + 3j, -2 + 8j])
# create a handler
efh = pylbo.EigenfunctionHandler(ds)
# retrieve the eigenfunctions based on the indices
eigenfuncs = efh.get_eigenfunctions(idxs)
```
Now, `eigenfuncs` will be a Numpy-array of size 2 (since you provided 2 eigenvalue guesses), and every index corresponds
to the index of your guess. Meaning, `eigenfuncs[0]` corresponds to the eigenfunctions of `1 + 3j`, and so on.
Every element of the `eigenfuncs` array is a dictionary, containing all eigenfunctions as well as the eigenvalues.
To retrieve what you need, simply do
```python 
# for the first eigenvalue:
rho_ef1 = eigenfuncs[0].get("rho")
# for the second eigenvalue:
rho_ef2 = eigenfuncs[1].get("rho")
# to see which keys are in there:
print(eigenfuncs[0].keys())
>> "rho", "v1", "v2", "v3", "T", "a1", "a2", "a3", "eigenvalue"
```
The names of the keys are self-explanatory.


## Analysing a series of files
### Plotting the spectrum
Say you loaded a series of datfiles, which are contained in an array called `datasets` (as above).
If you did a multirun of one of the pre-implemented equilibrium submodules you can create a combined 
spectrum by doing
```python 
multispec = pylbo.MultiSpectrum(datasets)
multispec.plot_precoded_run(annotate_continua=False)
pylbo.show()
```
For user-defined equilibria you'll have to do this manually. You can annotate the continuum ranges for all
datasets on the figure by setting `annotate_continua=True`. Note that this is not (yet) interactive.

**Under development** <i class="fas fa-hard-hat"></i> <i class="fas fa-hammer"></i>  
Due to Pylbo's beta status the [`MultiSpectrum`](../../src-docs/sphinx/autoapi/pylbo/visualisations/spectrum/index.html#pylbo.visualisations.spectrum.MultiSpectrum)
class is still under development, and will be extended with additional functionality in the future.
{: .notice--danger}

### Plotting eigenfunctions
**Under development** <i class="fas fa-hard-hat"></i> <i class="fas fa-hammer"></i>  
For the moment not possible on a multispectrum, but will come in the future.
{: .notice--danger}