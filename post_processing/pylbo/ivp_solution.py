import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import SymLogNorm, TwoSlopeNorm

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
            return comp_array, ""  # no scaling if no units
        
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
            label  = r"[dyn cm$^{-2}$]"
        else:
            factor = 1.0
            label  = ""

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

    # -----------------------------------------------------------------
    # NEW METHOD:  growth curve with optional analytic e-fold overlay
    # -----------------------------------------------------------------
    def plot_growth_vs_time(
        self,
        component,
        mode="centre",          # "centre" | "max" | "integral"
        centre_idx=None,        # if None, use n_points//2
        region=None,            # (i_min, i_max) for "integral" or "max"
        logy=True,
        tau=None,               # analytic e-folding time in *seconds* (float)
        ax=None,
        data_kw=None,           # kwargs for data curve
        fit_kw=None,            # kwargs for analytic curve
    ):
        """
        Plot amplitude of a perturbation vs. time. If `tau` is provided,
        overlay an analytic ±exp(t/tau) curve for comparison.

        Parameters
        ----------
        component : str | int
            Primitive variable ("rho", "T", "v1", ...).
        mode : {"centre","max","integral"}
            How to reduce spatial data to a scalar amplitude.
        centre_idx : int
            Spatial index for mode="centre" (default midpoint).
        region : tuple (i_min, i_max)
            Slice for mode="max" or "integral". Default: whole domain.
        logy : bool
            Use semilog-y axis.
        tau : float, optional
            Analytic e-folding time [s]. Positive => growth; negative => decay.
        ax : matplotlib Axes, optional
            Supply your own axes.
        data_kw / fit_kw : dict
            Extra kwargs forwarded to `plot` / `semilog` for data and fit.
        """
        data_kw = {} if data_kw is None else data_kw
        fit_kw  = {} if fit_kw  is None else fit_kw

        comp_arr = self.get_component(component)             # (n_snap,n_pts)
        times_phys = self._scale_time_array(np.asarray(self.times))

        # --- spatial reduction -------------------------------------------------
        if region is not None:
            i_min, i_max = region
            sl = slice(i_min, i_max)
        else:
            sl = slice(None)

        if mode == "centre":
            if centre_idx is None:
                centre_idx = comp_arr.shape[1] // 2
            amp = comp_arr[:, centre_idx]
        elif mode == "max":
            amp = np.max(np.abs(comp_arr[:, sl]), axis=1)
        elif mode == "integral":
            amp = np.sqrt(np.sum(comp_arr[:, sl] ** 2, axis=1))
        else:
            raise ValueError("mode must be 'centre', 'max' or 'integral'.")

        amp_phys, label = self._scale_component_array(component, amp)

        # --- plotting ----------------------------------------------------------
        if ax is None:
            fig, ax = plt.subplots()

        plot_fn = ax.semilogy if logy else ax.plot
        plot_fn(
            times_phys,
            amp_phys,
            label=f"{component}",
            **data_kw
        )

        # --- analytic overlay ---------------------------------------------------
        if tau is not None and tau != 0.0:
            a0 = amp_phys[0]
            t0 = times_phys[0]
            fit = a0 * np.exp((times_phys - t0) / tau)
            plot_fn(
                times_phys,
                fit,
                linestyle="--",
                label=f"analytic  exp($t/{tau:.0f}$ s)",
                **fit_kw,
            )

        # --- cosmetics ----------------------------------------------------------
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(f"Amplitude {label}")
        ax.set_title(f"{component} amplitude vs. time")

        if mode == "centre":
            ax.annotate(
                f"centre idx = {centre_idx}",
                xy=(0.02, 0.95), xycoords="axes fraction",
                fontsize=8, va="top"
            )
        elif region is not None:
            ax.annotate(
                f"region = [{sl.start}:{sl.stop}]",
                xy=(0.02, 0.95), xycoords="axes fraction",
                fontsize=8, va="top"
            )

        ax.legend()
        return ax

    # ------------------------------------------------------------------
    #   ENERGY diagnostic :  E_kin , E_int , E_tot  vs time
    # ------------------------------------------------------------------
    def plot_energy_timeseries(
            self,
            rho0, T0,                       # 1-D background arrays
            p0=None, gamma=5/3,
            v_comp="v1", rho_comp="rho", T_comp="T",
            region=None,
            logy=True,
            energies=("kin", "int", "tot"),
            ax=None,
            plot_kw=None,
            return_data=False,
            to_physical=False,              # ← NEW
            energy_unit=None):              # ← NEW
        """
        Plot (or return) kinetic, internal, and/or total perturbed energies.

        Parameters
        ----------
        to_physical : bool, optional
            If True, convert energies to cgs ergs using `energy_unit`
            or the units stored in `self.units`.
        energy_unit : float, optional
            Conversion factor (erg per code-unit energy).  Overrides auto-derivation.
        """
        # ------------------------------------------------------------------
        # 0. sanity
        valid = {"kin", "int", "tot"}
        energies = tuple(e for e in energies if e in valid)
        if not energies:
            raise ValueError(f"`energies` must contain at least one of {valid}")
        plot_kw = {} if plot_kw is None else dict(plot_kw)
        sl = slice(*region) if region else slice(None)

        # ------------------------------------------------------------------
        # 1. primitive perturbations
        v1   = self.get_component(v_comp)[:, sl]
        rho1 = self.get_component(rho_comp)[:, sl]
        T1   = self.get_component(T_comp)[:, sl]

        rho0_sl, T0_sl = rho0[sl][None, :], T0[sl][None, :]

        # ------------------------------------------------------------------
        # 2. specific heat
        c_v = 1.0 / (gamma - 1.0)

        # ------------------------------------------------------------------
        # 3. metric ds  (physical cm)
        if self.x_domain is not None:                 # non-uniform grid
            x_cm = self._scale_x_domain(self.x_domain)[sl] * 1e8
            ds   = np.gradient(x_cm)
        else:                                         # uniform grid
            ds   = np.ones_like(rho0_sl[0])
        ds = ds[None, :]                              # broadcast

        # ------------------------------------------------------------------
        # 4. energies  (still code-unit values for now)
        Ekin = 0.5 * np.sum(rho0_sl * v1**2                 * ds, axis=1)
        Eint =       np.sum(c_v * (rho0_sl*T1 + T0_sl*rho1) * ds, axis=1)
        Etot = Ekin + Eint

        # ------------------------------------------------------------------
        # 5. unit conversion (optional)
        if to_physical:
            if energy_unit is None:
                # --- try to derive from self.units ------------------------
                try:
                    rho_ref = self.units["unit_density"]          # g cm⁻³
                    T_ref   = self.units["unit_temperature"]      # K
                    L_ref   = self.units["unit_length"]           # cm
                    mu      = self.units.get("mean_molecular_weight", 1.0)
                    # physical gas constant per gram
                    k_B   = 1.380649e-16    # erg K⁻¹
                    m_H   = 1.6735575e-24   # g
                    R_phys = k_B / (mu * m_H)
                    v_ref = (R_phys * T_ref) ** 0.5               # cm s⁻¹
                    energy_unit = rho_ref * v_ref**2 * L_ref      # erg
                except Exception as exc:
                    raise RuntimeError(
                        "Cannot auto-derive `energy_unit`; "
                        "provide it explicitly."
                    ) from exc
            # apply conversion
            Ekin *= energy_unit
            Eint *= energy_unit
            Etot *= energy_unit
            y_label = "Energy [erg]"
        else:
            y_label = "Energy [code units]"

        # ------------------------------------------------------------------
        # 6. time array (physical seconds)
        t_phys = self._scale_time_array(np.asarray(self.times))

        # ------------------------------------------------------------------
        # 7. plotting
        if ax is None:
            _, ax = plt.subplots()

        if logy:
            pos_vals = np.concatenate([arr[arr > 0] for arr in (Ekin, Eint, Etot)])
            if pos_vals.size > 0:
                linth = 0.01 * pos_vals.min()
                ax.set_yscale("symlog", linthresh=linth, linscale=1.0)
                from matplotlib.ticker import NullFormatter
                ax.yaxis.set_minor_formatter(NullFormatter())

        _style = dict(
            kin=dict(color="tab:blue",  ls="--", label=r"$E_{\mathrm{kin}}$"),
            int=dict(color="tab:red",   ls="-.", label=r"$E_{\mathrm{int}}$"),
            tot=dict(color="tab:green", ls="-",  label=r"$E_{\mathrm{tot}}$"),
        )

        for key in energies:
            ax.plot(t_phys, {"kin": Ekin, "int": Eint, "tot": Etot}[key],
                    **{**_style[key], **plot_kw})

        ax.set_xlabel("Time [s]")
        ax.set_ylabel(y_label)
        if ax.get_title() == "":
            ax.set_title("Perturbed energy vs. time")
        if any(line.get_label() for line in ax.lines):
            ax.legend()

        if return_data:
            return ax, {"t": t_phys, "Ekin": Ekin, "Eint": Eint, "Etot": Etot}
        return ax


    def plot_derived_heatmap(
        self,
        kind,                          # "pressure" | "entropy"
        rho0, T0, p0=None,             # 1-D background (code units)
        gamma=5.0/3.0,
        cmap="plasma",
        ax=None,
        time_range=None,
        **imshow_kw,
    ):
        """
        Draw a space-time heat-map for a derived perturbation field.
        """
        # primitive perturbations
        rho1 = self.get_component("rho")
        T1   = self.get_component("T")

        rho0_b, T0_b = rho0.reshape(1, -1), T0.reshape(1, -1)

        if kind.lower() == "pressure":
            field   = rho0_b * T1 + T0_b * rho1          # p1
            label   = r"[code-pressure]" if self.units is None else r"[dyn cm$^{-2}$]"

        elif kind.lower() == "entropy":
            if p0 is None:
                p0 = rho0 * T0
            p0_b  = p0.reshape(1, -1)
            p1    = rho0_b * T1 + T0_b * rho1
            field = p1 / p0_b - gamma * rho1 / rho0_b    # S1 (dimensionless)
            label = ""
        else:
            raise ValueError("kind must be 'pressure' or 'entropy'")

        # unit scaling only if a unit dict is present
        if self.units is not None and kind.lower() == "pressure":
            field, _ = self._scale_component_array("p", field)

        # ------------------------------------------------------------------
        # put the derived field temporarily at the end of self.data
        idx_tmp = self.data.shape[1]
        self.data = np.concatenate([self.data, field[:, None, :]], axis=1)
        self.component_names[idx_tmp] = f"derived_{kind}"

        vmin = imshow_kw.pop("vmin", field.min())
        vmax = imshow_kw.pop("vmax", field.max())

        if ax is None:
            fig, ax = plt.subplots()

        try:
            self.plot_space_time_heatmap(
                idx_tmp, ax=ax, cmap=cmap, time_range=time_range,
                vmin=vmin, vmax=vmax,
                **imshow_kw
            )
            ax.set_title(f"{kind.capitalize()} perturbation")
            # replace auto colour-bar label
            ax.images[0].colorbar.set_label(f"{kind} {label}")
        finally:
            # clean-up: remove temporary column
            self.data = self.data[:, :-1, :]
            del self.component_names[idx_tmp]

        return ax

    def plot_heat_loss_function(
            self,
            rho0, T0,
            Lambda_func,
            H_func=None,
            cmap="coolwarm",
            ax=None,
            time_range=None,
            logscale=False,
            **imshow_kw
        ):
        r"""
        Plot a space–time heatmap of the volumetric net heat-loss term
            ρ L(ρ, T) = ρ² Λ(T) − H(ρ, T).

        Parameters
        ----------
        rho0, T0 : 1-D arrays
            Background profiles (code units).
        Lambda_func : callable
            Cooling curve Λ(T) [erg cm³ s⁻¹] (or in code units).
            For tabulated Colgan_DM data, pass an interpolator created
            with scipy.interpolate.interp1d.
        H_func : callable, optional
            Heating function H(ρ, T) with the same units as ρ² Λ(T).
            If None, assumes static equilibrium
            (H = ρ0² Λ(T0)), so that ρ L = 0 initially.
        cmap : str, default "coolwarm"
            Colormap to use for deviations (e.g. blue = net heating,
            red = net cooling, depending on sign convention).
        ax : matplotlib.axes.Axes, optional
            Plot onto existing axes.
        time_range : (t_min, t_max), optional
            Restrict the plotted time interval.
        logscale : bool, default False
            If True, use a symmetric logarithmic colour normalization
            around zero (SymLogNorm). Otherwise, linear scaling with
            midpoint at zero (TwoSlopeNorm).
        **imshow_kw :
            Additional arguments passed to imshow (e.g. vmin/vmax, alpha).
        """

        # --- Reconstruct total fields --------------------------------------------
        rho1 = self.get_component("rho")
        T1   = self.get_component("T")

        rho0_b = rho0.reshape(1, -1)
        T0_b   = T0.reshape(1, -1)
        rho = rho0_b + rho1
        T   = T0_b + T1

        # --- Heating & Cooling ---------------------------------------------------
        if H_func is None:
            # equilibrium heating → balances initial cooling
            H = (rho0_b**2) * Lambda_func(T0_b)
        else:
            H = H_func(rho, T)

        # Volumetric net heat-loss term: ρ L = ρ² Λ(T) − H
        rhoL = (rho**2) * Lambda_func(T) - H     # shape (n_snap, n_points)

        # --- Optional unit scaling ----------------------------------------------
        if self.units is not None:
            t_unit = self.units.get("unit_time", 1.0)
            e_unit = self.units.get("unit_energy_density", 1.0)
            rhoL *= e_unit / t_unit

        # --- Axes setup ----------------------------------------------------------
        times_phys = self._scale_time_array(np.asarray(self.times))
        if self.x_domain is not None:
            x_phys = self._scale_x_domain(self.x_domain)
        else:
            x_phys = np.arange(rhoL.shape[1])

        if time_range is not None:
            tmin, tmax = time_range
            mask = (times_phys >= tmin) & (times_phys <= tmax)
            rhoL = rhoL[mask, :]
            times_phys = times_phys[mask]

        # --- Determine symmetric colour limits ----------------------------------
        vabs = np.nanmax(np.abs(rhoL))
        if logscale:
            linth = imshow_kw.pop("linthresh", 1e-6)
            norm = SymLogNorm(linthresh=linth, vmin=-vabs, vmax=vabs)
        else:
            norm = TwoSlopeNorm(vmin=-vabs, vcenter=0.0, vmax=vabs)

        # --- Plot ----------------------------------------------------------------
        if ax is None:
            fig, ax = plt.subplots()

        im = ax.imshow(
            rhoL,
            origin="lower",
            aspect="auto",
            cmap=cmap,
            norm=norm,
            extent=[x_phys[0], x_phys[-1], times_phys[0], times_phys[-1]],
            **imshow_kw,
        )

        cbar = plt.colorbar(im, ax=ax, format=ticker.FuncFormatter(fmt))
        cbar.set_label(r"$\rho\mathcal{L}$ [erg cm$^{-3}$ s$^{-1}$]")

        ax.set_xlabel(r"$s$ [Mm]")
        ax.set_ylabel("Time [s]")
        # ax.set_title(r"$\rho \mathcal{L}(\rho,T)$")

        return ax
