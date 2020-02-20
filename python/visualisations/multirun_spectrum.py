from data_container import LEGOLASDataContainer
from utilities.parameters import PRECODED_MULTIRUNS

import matplotlib.pyplot as plt
import numpy as np


class MultiSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.gridpts_title = ' ({} gridpts / run)'.format(self.datacontainer[0].gridpts)


    def plot(self):
        mr_name = self._get_precoded_name()

        # first check custom plot requirements
        if 'gravito_mhd' in mr_name:
            for data in self.datacontainer:
                w_prefact = np.sqrt(data.params['cte_rho0'])
                omega2 = w_prefact**2 * np.real(data.omegas**2)
                k02 = (data.params['k2'] ** 2 + data.params['k3'] ** 2) * np.ones_like(omega2)
                self.ax.plot(k02, omega2, '.b', markersize=2, alpha=0.8)
            self.ax.set_ylabel(r'$\dfrac{1}{v_A^2}\omega^2$')
            self.ax.set_xlabel('$k_0^2$')
            self.ax.set_title(mr_name + self.gridpts_title)

        elif 'interchange_modes' in mr_name:
            # for interchange modes we plot w**2 as a function of theta/pi
            thetas = np.linspace(0, np.pi, len(self.datacontainer)) / np.pi
            for idx, data in enumerate(self.datacontainer):
                # omega is normalised as in gravito_mhd
                w_prefact = np.sqrt(data.params['cte_rho0'])
                omega2 = w_prefact**2 * np.real(data.omegas**2)
                th = thetas[idx] * np.ones_like(omega2)
                self.ax.plot(th, omega2, '.b', markersize=2, alpha=0.8)
            self.ax.set_ylabel(r'$\dfrac{1}{v_A^2}\omega^2$')
            self.ax.set_xlabel(r'$\theta$ / $\pi$')
            self.ax.set_title(mr_name + self.gridpts_title)
            self.ax.set_ylim([-3, 14])
            self.ax.set_xlim([0, 1])

        elif 'constant_current' in mr_name:
            # for constant current modes we plot w**2 vs the safety factor q = 2*k/j0
            qfactors = np.linspace(1.9, 2.1, len(self.datacontainer))
            for idx, data in enumerate(self.datacontainer):
                omega2 = np.real(data.omegas**2)
                qfact = qfactors[idx] * np.ones_like(omega2)

                self.ax.plot(qfact, omega2, '.b', markersize=2, alpha=0.8)
                self.ax.set_yscale('symlog', linthreshy=1e-8)
            self.ax.set_xlabel(r"Safety factor $q = \frac{2k}{j0}$")
            self.ax.set_ylabel(r"$\omega^2$")
            self.ax.set_title(mr_name + self.gridpts_title)

        else:
            for data in self.datacontainer:
                omega2 = np.real(data.omegas**2)
                k02 = (data.params['k2'] ** 2 + data.params['k3'] ** 2) * np.ones_like(omega2)
                self.ax.plot(k02, omega2, '.b', markersize=2, alpha=0.8)
            self.ax.set_ylabel(r'$\omega^2$')
            self.ax.set_xlabel('$k_0^2$')
            self.ax.set_title(self.datacontainer[0].current_eq + self.gridpts_title)


    def _get_precoded_name(self):
        for mr_name in PRECODED_MULTIRUNS:
            if mr_name in self.datacontainer[0].fname_w:
                return mr_name
        else:
            return None
