import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

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

    def __init__(self, times, data, component_names=None, x_domain=None, units=None):
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
        self.units = units

    def _scale_component_array(self, component, comp_array):
        """
        Returns comp_array scaled into physical (cgs) units,
        based on the component name and the self.units dict.
        """
        if self.units is None:
            return comp_array  # no scaling if no units
        
        # Identify which scale factor to use
        comp_name = str(component).lower()
        if "rho" in comp_name:
            factor = self.units["unit_density"]             # g cm^-3
            label = r"[g cm$^{-3}$]"
        elif "v" in comp_name:
            factor = self.units["unit_velocity"]            # cm s^-1
            label = r"[cm s$^{-1}$]"
        elif "temp" in comp_name or "t" in comp_name:
            factor = self.units["unit_temperature"]  # K
            label = r"[K]"
        elif "p" in comp_name:
            factor = self.units["unit_pressure"]            # dyn cm^-2
        else:
            factor = 1.0  # fallback/no scaling
        
        return comp_array * factor, label

    def _scale_time_array(self, times):
        if self.units is None:
            return times
        # dimensionless * unit_time => seconds
        return times * self.units["unit_time"]

    def _scale_x_domain(self, x_vals):
        if self.units is None:
            return x_vals
        # dimensionless * unit_length => cm => * 10^8 => Mm
        return x_vals * self.units["unit_length"] * 1e-8

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

        # 1) Scale times -> physical units
        times_phys = self._scale_time_array(times)

        # 2) Clip times if time_range is given
        if time_range is not None:
            t_min, t_max = time_range
            mask = (times_phys >= t_min) & (times_phys <= t_max)
            comp_array = comp_array[mask, :]
            times_phys = times_phys[mask]
            t_min, t_max = times_phys[0], times_phys[-1]
        else:
            t_min, t_max = times_phys[0], times_phys[-1]

        # 3) Scale the solution array for the chosen component
        comp_array_phys, label = self._scale_component_array(component, comp_array)

        # 4) Prepare axes
        if ax is None:
            fig, ax = plt.subplots()

        # 5) Scale spatial domain
        n_points = comp_array_phys.shape[1]
        if self.x_domain is not None and len(self.x_domain) == n_points:
            x_vals_phys = self._scale_x_domain(self.x_domain)
            x_min, x_max = x_vals_phys[0], x_vals_phys[-1]
        else:
            x_min, x_max = 0, n_points - 1
            x_vals_phys = None

        im = ax.imshow(
            comp_array_phys,
            origin='lower',
            aspect='auto',
            cmap=cmap,
            extent=[x_min, x_max, t_min, t_max],
            **imshow_kwargs
        )
        cbar = plt.colorbar(im, ax=ax, format=ticker.FuncFormatter(fmt))
        cbar.set_label(f"{component} {label}")

        # 6) Axis labels with units
        if x_vals_phys is not None:
            ax.set_xlabel("x [Mm]")
        else:
            ax.set_xlabel("x index")
        ax.set_ylabel("Time [s]")
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
        comp_array_phys, label = self._scale_component_array(component, comp_array)

        times_phys = self._scale_time_array(np.array(self.times))

        n_snap, n_points = comp_array_phys.shape

        # Scale x-axis
        if self.x_domain is not None and len(self.x_domain) == n_points:
            x_vals_phys = self._scale_x_domain(self.x_domain)
            xlabel = "x [Mm]"
        else:
            x_vals_phys = np.arange(n_points)
            xlabel = "x index"

        if ax is None:
            fig, ax = plt.subplots()

        for snap_idx in snap_indices:
            if snap_idx < 0 or snap_idx >= n_snap:
                raise IndexError(f"Snapshot index {snap_idx} out of range (0..{n_snap-1}).")
            
            y_vals = comp_array_phys[snap_idx, :]
            t_val = times_phys[snap_idx]

            ax.plot(x_vals_phys, y_vals, label=f"t={t_val:.3f} s", **plot_kwargs)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"{component} {label}")
        ax.set_title(f"{component} vs. x for selected time steps")
        ax.legend()
        return ax

