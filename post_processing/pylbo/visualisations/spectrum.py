import matplotlib.pyplot as plt
import matplotlib.lines as mpl_lines
import matplotlib.patches as mpl_patches
import numpy as np
from ..utilities.exceptions import InconsistentMultirunFile
from ..visualisations.eigenfunctions import EigenfunctionHandler

class SingleSpectrum:
    def __init__(self, ds):
        if isinstance(ds, list):
            raise TypeError('SingleSpectrum needs a single LegolasDatacontainer instance, not a list.')
        self.ds = ds
        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.annotate_continua = True

    def plot_eigenfunctions(self, merge_figs=False, annotate_continua=True):
        if not annotate_continua == self.annotate_continua:
            self.annotate_continua = annotate_continua
        efh = EigenfunctionHandler(self, merge_figs)
        efh.connect_figure_events()

    def plot_spectrum(self, fig=None, ax=None, annotate_continua=True):
        if fig is not None and ax is not None:
            self.fig = fig
            self.ax = ax
        self.annotate_continua = annotate_continua
        self.ax.plot(self.ds.eigenvalues.real, self.ds.eigenvalues.imag, '.b', alpha=0.8)
        self.ax.axhline(y=0, linestyle='dotted', color='grey', alpha=0.3)
        self.ax.axvline(x=0, linestyle='dotted', color='grey', alpha=0.3)
        self.ax.set_xlabel(r"Re($\omega$)")
        self.ax.set_ylabel(r"Im($\omega$)")
        self.ax.set_title(self.ds.eq_type)
        if annotate_continua:
            self._annotate_continua(self.fig, self.ax)
        self.fig.canvas.draw()
        self.fig.tight_layout()
        return self.fig, self.ax

    def _annotate_continua(self, fig, ax):
        legend_items = []
        alpha_box = 0.3
        alpha_point = 0.5

        def on_legend_pick(event):
            artist = event.artist
            # do nothing if clicking outside of legend
            if not artist in regions:
                return
            item = regions[artist]
            visible = not item.get_visible()
            item.set_visible(visible)
            if visible:
                if isinstance(artist, mpl_lines.Line2D):
                    artist.set_alpha(alpha_point)
                else:
                    artist.set_alpha(alpha_box)
            else:
                artist.set_alpha(0.1)
            fig.canvas.draw()

        def draw_region(array, facecolor, legend_lbl=None, draw_horizontal=False):
            lb = np.min(array)
            rb = np.max(array)
            if np.abs(lb - rb) > 1e-10:
                if draw_horizontal:
                    legend_item = ax.axhspan(np.min(array), np.max(array), facecolor=facecolor,
                                             alpha=alpha_box, label=legend_lbl)
                else:
                    legend_item = ax.axvspan(np.min(array), np.max(array), facecolor=facecolor,
                                             alpha=alpha_box, label=legend_lbl)
            else:  # handles the case where continua are collapsed to a single point
                legend_item, = ax.plot([lb], 0, marker='p', markersize=8, color=facecolor, linestyle='none',
                                       alpha=alpha_point, label=legend_lbl)
            legend_items.append(legend_item)

        # draw regions on plot
        draw_region(self.ds.continua['wS+'], facecolor='red',
                    legend_lbl=r'$\Omega_S^+$ slow continuum')
        draw_region(self.ds.continua['wS-'], facecolor='red',
                    legend_lbl=r'$\Omega_S^-$ slow continuum')
        draw_region(self.ds.continua['wA+'], facecolor='cyan',
                    legend_lbl=r'$\Omega_A^+$ Alfvén continuum')
        draw_region(self.ds.continua['wA-'], facecolor='cyan',
                    legend_lbl=r'$\Omega_A^-$ Alfvén continuum')
        # thermal continuum is imaginary, so draw horizontally
        draw_region(self.ds.continua['wth'], facecolor='green',
                    legend_lbl=r'$\Omega_{th}$ thermal continuum', draw_horizontal=True)
        legend = ax.legend(loc='best')

        # retrieve legend patches and items for interactive enabling of continuum regions
        regions = {}
        patches = legend.get_patches()
        lines = legend.get_lines()
        for region in legend_items:
            if isinstance(region, mpl_patches.Polygon):
                l_item = patches.pop(0)
            elif isinstance(region, mpl_lines.Line2D):
                l_item = lines.pop(0)
            else:
                raise ValueError
            l_item.set_picker(10)
            l_item.set_alpha(0.1)
            region.set_visible(False)  # make continuum regions invisible by default
            regions[l_item] = region
        fig.canvas.mpl_connect('pick_event', on_legend_pick)

    def plot_continua(self, fig=None, ax=None):
        if fig is None or ax is None:
            fig, ax = plt.subplots(3, 1, figsize=(12, 8))
        ax[0].plot(self.ds.grid_gauss, self.ds.continua['wS-'], color='red',
                   linestyle='solid', label=r'$\Omega_S^-$ slow continuum')
        ax[0].plot(self.ds.grid_gauss, self.ds.continua['wS+'], color='red',
                   linestyle='dashed', label=r'$\Omega_S^+$ slow continuum')
        ax[1].plot(self.ds.grid_gauss, self.ds.continua['wA-'], color='blue',
                   linestyle='solid', label=r'$\Omega_A^-$ Alfvén continuum')
        ax[1].plot(self.ds.grid_gauss, self.ds.continua['wA+'], color='blue',
                   linestyle='dashed', label=r'$\Omega_A^+$ Alfvén continuum')
        ax[2].plot(self.ds.grid_gauss, self.ds.continua['wth'], color='green',
                   linestyle='solid', label=r'Im$(\Omega_{th})$ thermal continuum')
        ax[2].set_xlabel('Grid')
        for ax_panel in ax:
            ax_panel.legend(loc='best')
            ax_panel.set_xlim([self.ds.grid_gauss[0], self.ds.grid_gauss[-1]])
        fig.subplots_adjust(hspace=0.05)
        fig.tight_layout()
        return fig, ax

    def plot_equilibria(self, fig=None, ax=None):
        if fig is None or ax is None:
            fig, ax = plt.subplots(1, figsize=(12, 8))
        for idx, var in enumerate(self.ds.header.get('equil_names')):
            values = self.ds.equilibria.get(var)
            if (values == 0).all():
                continue
            ax.plot(self.ds.grid_gauss, values, '-o', markersize=2, label=var, alpha=0.8)
        ax.set_xlabel('grid')
        ax.set_ylabel('values')
        ax.set_title('Equilibrium settings')
        ax.legend(loc='best')
        fig.canvas.draw()
        return fig, ax


class MultiSpectrum:
    def __init__(self, datasets):
        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.datasets = datasets
        self.slow_cont_line = {'xvals': [], 'min': [], 'max': []}
        self.alfv_cont_line = {'xvals': [], 'min': [], 'max': []}

    def _check_datasets(self, equilibrium):
        for ds in self.datasets:
            if not equilibrium == ds.header.get('eq_type'):
                raise InconsistentMultirunFile(ds.datfile, expected=equilibrium, found=ds.header.get('eq_type'))

    def plot_precoded_run(self, annotate_continua=False):
        equilibrium = self.datasets[0].header.get('eq_type')
        self._check_datasets(equilibrium)
        if equilibrium == 'gravito_acoustic':
            self._plot_gravito_acoustic()
        elif equilibrium == 'gravito_mhd':
            self._plot_gravito_mhd()
        elif equilibrium == 'interchange_modes':
            self._plot_interchange_modes()
        elif equilibrium == 'constant_current_tokamak':
            self._plot_constant_current()
        elif equilibrium in ('photospheric_flux_tube', 'coronal_flux_tube'):
            self._plot_flux_tube()
        else:
            raise ValueError('Automatically plotting the spectrum failed, should be done manually.')
        if annotate_continua:
            self._annotate_continua()
        return self.fig, self.ax

    def _annotate_continua(self):
        # annotate slow continuum
        xvals, slowmin, slowmax = self.slow_cont_line.values()
        self.ax.plot(xvals, slowmin, color='red', alpha=0.8)
        self.ax.plot(xvals, slowmax, color='red', alpha=0.8)
        self.ax.fill_between(xvals, slowmin, slowmax, color='red', alpha=0.3, label='slow_continuum')
        # annotate alfven continuum
        xvals, alfvmin, alfvmax = self.alfv_cont_line.values()
        self.ax.plot(xvals, alfvmin, color='cyan', alpha=0.8)
        self.ax.plot(xvals, alfvmax, color='cyan', alpha=0.8)
        self.ax.fill_between(xvals, alfvmin, alfvmax, color='cyan', alpha=0.3, label='Alfven continuum')
        self.ax.legend(loc='best')

    def _set_continua_lines(self, ds, x_value, prefactor=1, squared=True):
        exp = 2
        if not squared:
            exp = 1
        slow = prefactor * np.real(ds.continua['wS+']**exp)
        alfv = prefactor * np.real(ds.continua['wA+']**exp)
        for vals, cont in zip((slow, alfv), (self.slow_cont_line, self.alfv_cont_line)):
            cont['xvals'].append(x_value)
            cont['min'].append(np.min(vals))
            cont['max'].append(np.max(vals))

    def _plot_gravito_acoustic(self):
        for ds in self.datasets:
            cs = ds.get_sound_speed()[0]
            w_sq = (1 / cs**2) * np.real(ds.eigenvalues**2)
            w_sq[np.where(np.abs(w_sq) < 1e-8)] = np.nan
            k0_sq =  ds.get_k0_squared() * np.ones_like(w_sq)
            self.ax.plot(k0_sq, w_sq, '.b', markersize=3, alpha=0.8)
            self._set_continua_lines(ds, x_value=k0_sq[0], prefactor=(1 / cs**2))
        self.ax.set_ylabel(r'$\dfrac{1}{c^2}\omega^2$')
        self.ax.set_xlabel('$k_0^2$')
        self.ax.set_title('Gravito acoustic multirun')
        self.ax.set_xlim([0, 500])
        self.ax.set_ylim([0, 500])

    def _plot_gravito_mhd(self):
        for ds in self.datasets:
            # alfven speed is equal everywhere
            vA = ds.get_alfven_speed()[0]
            w_sq = (1 / vA**2) * np.real(ds.eigenvalues**2)
            w_sq[np.where(np.abs(w_sq) < 1e-8)] = np.nan
            k0_sq = ds.get_k0_squared() * np.ones_like(w_sq)
            self.ax.plot(k0_sq, w_sq, '.b', markersize=3, alpha=0.8)
            self._set_continua_lines(ds, x_value=k0_sq[0], prefactor=(1 / vA**2))
        self.ax.set_ylabel(r'$\dfrac{1}{v_A^2}\omega^2$')
        self.ax.set_xlabel('$k_0^2$')
        self.ax.set_title('Gravito mhd multirun')

    def _plot_interchange_modes(self):
        # here we plot w**2 as a function of theta/pi
        thetas = np.linspace(0, np.pi, len(self.datasets)) / np.pi
        for idx, ds in enumerate(self.datasets):
            vA = ds.get_alfven_speed()[0]
            w_sq = (1 / vA**2) * np.real(ds.eigenvalues**2)
            w_sq[np.where(np.abs(w_sq) < 1e-8)] = np.nan
            th = thetas[idx] * np.ones_like(w_sq)
            self.ax.plot(th, w_sq, '.b', markersize=2, alpha=0.8)
            self._set_continua_lines(ds, x_value=thetas[idx], prefactor=(1 / vA**2))
        self.ax.set_ylabel(r'$\dfrac{1}{v_A^2}\omega^2$')
        self.ax.set_xlabel(r'$\theta$ / $\pi$')
        self.ax.set_title('Interchange modes multirun')
        self.ax.set_ylim([-4.1, 14.4])
        self.ax.set_xlim([-0.01, 1.01])
        self.ax.set_xticks(np.arange(0, 1.2, 0.2))
        self.ax.set_yticks(np.arange(-4, 16, 2))

    def _plot_constant_current(self):
        # here we plot w**2 vs the safety factor q = 2*k3/j0
        qfactors = np.linspace(1.9, 2.1, len(self.datasets))
        for idx, ds in enumerate(self.datasets):
            w_sq = np.real(ds.eigenvalues**2)
            qfact = qfactors[idx] * np.ones_like(w_sq)
            self.ax.plot(qfact, w_sq, '.b', markersize=2, alpha=0.8)
            self.ax.set_yscale('symlog', linthreshy=1e-8)
            self._set_continua_lines(ds, x_value=qfactors[idx])
        self.ax.set_xlabel(r"Safety factor $q = \frac{2k}{j0}$")
        self.ax.set_ylabel(r"$\omega^2$")
        self.ax.set_title('Constant current m=-2 multirun')

    def _plot_flux_tube(self):
        # these values should be equal for all runs
        a = self.datasets[0].parameters['r0']
        gamma = self.datasets[0].gamma
        rho0 = self.datasets[0].parameters['cte_rho0']
        p0 = self.datasets[0].parameters['cte_p0']
        B0 = 2 * np.sqrt(gamma * p0)
        if self.datasets[0].eq_type == 'photospheric_flux_tube':
            rhoe = 8 * rho0 * (2 * gamma + 1) / (gamma + 18)
            pe = 18 * p0 * (2 * gamma + 1) / (gamma + 18)
            Be = np.sqrt(2 * gamma * p0 * (2 * gamma + 1) / (gamma + 18))
        else:
            rhoe = 4 * rho0 * (2 * gamma + 1) / (50 * gamma + 1)
            pe = p0 * (2 * gamma + 1) / (50 * gamma + 1)
            Be = 10 * np.sqrt(gamma * p0 * (2 * gamma + 1) / (50 * gamma + 1))
        cs = np.sqrt(gamma * p0 / rho0)
        cse = np.sqrt(gamma * pe / rhoe)
        cA = np.sqrt(B0 ** 2 / rho0)
        cAe = np.sqrt(Be ** 2 / rhoe)
        ct = cs * cA / np.sqrt(cs ** 2 + cA ** 2)
        cte = cse * cAe / np.sqrt(cse ** 2 + cAe ** 2)
        ck = np.sqrt((rho0 * cA ** 2 + rhoe * cAe ** 2) / (rho0 + rhoe))
        kz_values = np.linspace(0.1, 6.2, len(self.datasets))
        for idx, ds in enumerate(self.datasets):
            # phase speed c = w / kz
            c = np.abs(ds.eigenvalues.real) / kz_values[idx]
            kza = kz_values[idx] * a * np.ones_like(c)
            self.ax.plot(kza, c / cs, '.b', markersize=2)
        self.ax.axhline(cA / cs, label='$c_A$', lw=1, color='black', alpha=0.6)
        self.ax.axhline(cAe / cs, label='$c_{Ae}$', lw=1, color='black', linestyle='dotted', alpha=0.6)
        self.ax.axhline(cs / cs, label='$c_s$', lw=1, color='red', alpha=0.6)
        self.ax.axhline(cse / cs, label='$c_{se}$', lw=1, color='red', linestyle='dotted', alpha=0.6)
        self.ax.axhline(ct / cs, label='$c_t$', lw=1, color='green', alpha=0.6)
        self.ax.axhline(cte / cs, label='$c_{te}$', lw=1, color='green', linestyle='dotted', alpha=0.6)
        self.ax.axhline(ck / cs, label='$c_k$', lw=1, color='cyan', alpha=0.6)
        self.ax.legend(loc='best')
        self.ax.set_xlabel(r"$k_za$")
        self.ax.set_ylabel(r"$\omega / k_z c_s$")
        self.ax.set_title(self.datasets[0].eq_type.replace("_", " ") + ' multirun')
