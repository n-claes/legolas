from data_container import LEGOLASDataContainer
from utilities.parameters import DEFAULT_PARAMS

import matplotlib.pyplot as plt
import numpy as np

class MultiSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.varied_param = self._get_varied_param()


    def plot(self):
        for data in self:
            param_array = data.params[self.varied_param] * np.ones(len(data.omegas))
            omega2 = np.real(data.omegas)**2 + np.imag(data.omegas)**2

            self.ax.plot(param_array, omega2, '.b', alpha=0.8)
        self.ax.set_ylabel('$\omega^2$')
        self.ax.set_xlabel(self.varied_param)
        self.ax.set_title(self[0].current_eq)



    def _get_varied_param(self):
        varied_param = None

        for param in DEFAULT_PARAMS[self[0].current_eq]:
            param_string = '_{}_'.format(param)
            for data in self:
                if param_string in data.fname_w:
                    if varied_param is None:
                        varied_param = param
                    elif param != varied_param:
                        raise ValueError('Error with varied parameter, wrong config file in list? \n'
                                         'Previously detected varied parameter: {} \n'
                                         'Encountered different varied parameter: {}'.format(varied_param, param))
        return varied_param

