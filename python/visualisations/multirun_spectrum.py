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
        if 'gravito_acoustic' in mr_name:
            for data in self.datacontainer:
                # prefactor in this case is a/c where a is domain length (=1) and c sound speed (c2 = gamma*p0/rho0)
                w_prefact = 1.0 / np.sqrt(data.gamma * data.params['cte_p0'] / data.params['cte_rho0'])
                omega2 = w_prefact**2 * np.real(data.omegas**2)
                k02 = (data.params['k2']**2 + data.params['k3']**2) * np.ones_like(omega2)
                self.ax.plot(k02, omega2, '.b', markersize=2, alpha=0.8)
            self.ax.axhline(y=100, linestyle='dotted', color='grey', lw=1)
            self.ax.set_ylabel(r'$\dfrac{1}{c^2}\omega^2$')
            self.ax.set_xlabel('$k_0^2$')
            self.ax.set_title(mr_name + self.gridpts_title)
            self.ax.set_xlim([0, 500])
            self.ax.set_ylim([0, 500])

        elif 'gravito_mhd' in mr_name:
            for data in self.datacontainer:
                w_prefact = np.sqrt(data.params['cte_rho0'])
                omega2 = w_prefact**2 * np.real(data.omegas**2)
                k02 = (data.params['k2']**2 + data.params['k3']**2) * np.ones_like(omega2)
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
            self.ax.set_ylim([-4.1, 14.4])
            self.ax.set_xlim([-0.01, 1.01])
            self.ax.set_xticks(np.arange(0, 1.2, 0.2))
            self.ax.set_yticks(np.arange(-4, 16, 2))

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

        elif mr_name in ('photospheric_flux_tube', 'coronal_flux_tube'):
            a = self.datacontainer[0].params['r0']
            gamma = self.datacontainer[0].gamma
            rho0 = self.datacontainer[0].params['cte_rho0']
            p0 = self.datacontainer[0].params['cte_p0']
            B0 = 2 * np.sqrt(gamma * p0)
            if mr_name == 'photospheric_flux_tube':
                rhoe = 8 * rho0 * (2 * gamma + 1) / (gamma + 18)
                pe = 18 * p0 * (2 * gamma + 1) / (gamma + 18)
                Be = np.sqrt(2 * gamma * p0 * (2 * gamma + 1) / (gamma + 18))
            else:
                rhoe = 4 * rho0 * (2 * gamma + 1) / (50 * gamma + 1)
                pe = p0 * (2 * gamma + 1) / (50 * gamma + 1)
                Be = 10 * np.sqrt(gamma * p0 * (2 * gamma + 1) / (50 * gamma + 1))
            cs = np.sqrt(gamma * p0 / rho0)
            cse = np.sqrt(gamma * pe / rhoe)
            cA = np.sqrt(B0**2 / rho0)
            cAe = np.sqrt(Be**2 / rhoe)
            ct = cs*cA / np.sqrt(cs**2 + cA**2)
            cte = cse*cAe / np.sqrt(cse**2 + cAe**2)
            ck = np.sqrt((rho0 * cA**2 + rhoe * cAe**2)/(rho0 + rhoe))
            kz_values = np.linspace(0.1, 6.2, len(self.datacontainer))
            for idx, data in enumerate(self.datacontainer):
                # phase speed c = w / kz
                c = np.abs(data.omegas.real) / kz_values[idx]
                kza = kz_values[idx] * a * np.ones_like(c)
                self.ax.plot(kza, c/cs, '.b', markersize=2)
            self.ax.axhline(cA/cs, label='$c_A$', lw=1, color='black', alpha=0.6)
            self.ax.axhline(cAe/cs, label='$c_{Ae}$', lw=1, color='black', linestyle='dotted', alpha=0.6)
            self.ax.axhline(cs/cs, label='$c_s$', lw=1, color='red', alpha=0.6)
            self.ax.axhline(cse/cs, label='$c_{se}$', lw=1, color='red', linestyle='dotted', alpha=0.6)
            self.ax.axhline(ct/cs, label='$c_t$', lw=1, color='green', alpha=0.6)
            self.ax.axhline(cte/cs, label='$c_{te}$', lw=1, color='green', linestyle='dotted', alpha=0.6)
            self.ax.axhline(ck/cs, label='$c_k$', lw=1, color='cyan', alpha=0.6)
            self.ax.legend(loc='best')
            self.ax.set_xlabel(r"$k_za$")
            self.ax.set_ylabel(r"$\omega / k_z c_s$")
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
