import numpy as np
import matplotlib.pyplot as plt

class IVPSolution:
    """
    A container for IVP data from Legolas, with built-in plotting methods.

    Attributes
    ----------
    times : np.ndarray (shape = (n_snap,))
        Physical times at each snapshot.
    data : np.ndarray (shape = (n_snap, n_comp, n_points))
        The solution array for each snapshot, each component, and each spatial point.
    component_names : dict or list
        Map from integer component index -> name. E.g. {1: "rho", 2: "v1", 3: "T1"}.
        If you prefer zero-based indexing, store it accordingly.
    x_domain : np.ndarray or None
        Optional array for the spatial axis, shape = (n_points,).
        Or store x_start / x_end if you prefer a continuous domain.
    """

    def __init__(self, times, data, component_names=None, x_domain=None):
        """
        Parameters
        ----------
        times : np.ndarray of shape (n_snap,)
            Physical times for each snapshot.
        data : np.ndarray of shape (n_snap, n_comp, n_points)
            The 3D snapshot data.
        component_names : dict or list, optional
            Mapping from integer->string. 
        x_domain : np.ndarray, optional
            If provided, a 1D array with the spatial coordinate for each point.
        """
        self.times = times                  # shape: (n_snap,)
        self.data = data                    # shape: (n_snap, n_comp, n_points)
        self.component_names = component_names or {}
        self.x_domain = x_domain            # shape: (n_points,) or None

    def _get_component_index(self, component):
        """
        Internal helper to interpret 'component' as either an int or a string name.
        If it's a string (like "rho"), return the matching integer index.
        If it's already an integer, just return it (with optional checks).
        """
        if isinstance(component, int):
            return component
        elif isinstance(component, str):
            # search in self.component_names
            for idx, nm in self.component_names.items():
                if nm == component:
                    return idx
            raise ValueError(f"No component named '{component}' found in {self.component_names}")
        else:
            raise TypeError("component must be either an int or a str")

    def get_component(self, component):
        """
        Returns (n_snap, n_points) array for the requested 'component'.
        'component' can be int or str.

        Example:
        --------
        rho_data = ivp_sol.get_component("rho")
        # shape(rho_data) = (n_snap, n_points)
        """
        comp_idx = self._get_component_index(component)
        return self.data[:, comp_idx, :]

    # -------------------------------------------------------------------------
    # PLOTTING METHODS
    # -------------------------------------------------------------------------

    def plot_space_time_heatmap(self, component, ax=None, cmap='plasma', time_range=None, **imshow_kwargs):
        """
        Plot a 2D heatmap (time vs. space) for the given component.

        The vertical axis is time, the horizontal axis is the domain coordinate.
        If x_domain is None, the domain axis will just be 0..n_points-1.

        Parameters
        ----------
        component : int or str
            Which component to plot.
        ax : matplotlib.axes.Axes, optional
            If provided, uses that axes; else creates a new figure.
        cmap : str
            Colormap to use for imshow.
        time_range : tuple (t_min, t_max), optional
            Time window to plot. If None, plots the full time span.
        imshow_kwargs : dict
            Additional arguments to pass to imshow (e.g. vmin, vmax, etc.).
        """
        comp_array = self.get_component(component)  # shape (n_snap, n_points)
        times = np.array(self.times)

        # Apply time slicing if a time range is specified
        if time_range is not None:
            t_min, t_max = time_range
            time_mask = (times >= t_min) & (times <= t_max)
            comp_array = comp_array[time_mask, :]
            times = times[time_mask]
        else:
            t_min, t_max = times[0], times[-1]

        if ax is None:
            fig, ax = plt.subplots()

        n_points = comp_array.shape[1]

        if self.x_domain is not None and len(self.x_domain) == n_points:
            x_min, x_max = self.x_domain[0], self.x_domain[-1]
        else:
            x_min, x_max = 0, n_points - 1

        im = ax.imshow(
            comp_array,
            origin='lower',
            aspect='auto',
            cmap=cmap,
            extent=[x_min, x_max, t_min, t_max],
            **imshow_kwargs
        )
        plt.colorbar(im, ax=ax, label=str(component))

        ax.set_xlabel("x coordinate" if self.x_domain is not None else "x index")
        ax.set_ylabel("Time")
        ax.set_title(f"Space-Time Heatmap of {component}")
        return ax
    
    def plot_spatial_slices(self, component, snap_indices, ax=None, **plot_kwargs):
        """
        For each time index in snap_indices, plot the entire x range of the
        chosen component as a line plot. The horizontal axis is the physical x
        (if self.x_domain is provided), otherwise it's 0..n_points-1.

        Parameters
        ----------
        component : int or str
            Which component to plot (e.g. "rho").
        snap_indices : list of int
            Which snapshot indices to display. Indices in [0..n_snap-1].
        ax : matplotlib.axes.Axes, optional
            If provided, use this axes to plot. Otherwise, create a new figure.
        plot_kwargs : dict
            Extra keyword args to pass to matplotlib's plot() function.

        Example
        -------
        ivp_sol.plot_spatial_slices("rho", snap_indices=[0, 5, 10])
        => This draws 3 lines (time=times[0], time=times[5], time=times[10]),
           each line is x vs. rho.
        """
        comp_array = self.get_component(component)   # shape (n_snap, n_points)
        n_snap, n_points = comp_array.shape

        # pick the x axis
        if self.x_domain is not None and len(self.x_domain) == n_points:
            x_vals = self.x_domain
            xlabel = "x coordinate"
        else:
            x_vals = np.arange(n_points)
            xlabel = "x index"

        if ax is None:
            fig, ax = plt.subplots()

        for snap_idx in snap_indices:
            if snap_idx < 0 or snap_idx >= n_snap:
                raise IndexError(f"Snapshot index {snap_idx} out of range (0..{n_snap-1}).")

            # data at this time step -> shape (n_points,)
            y_vals = comp_array[snap_idx, :]
            t_str = f"{self.times[snap_idx]:.3f}"
            ax.plot(x_vals, y_vals, label=f"time={t_str}", **plot_kwargs)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(str(component))
        ax.set_title(f"{component} vs. x for selected time steps")
        ax.legend()
        return ax
