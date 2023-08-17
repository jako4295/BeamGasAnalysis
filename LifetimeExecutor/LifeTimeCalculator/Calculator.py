import numpy as np
import pandas as pd
import scipy
from scipy.optimize import curve_fit

from DataHandler import DataObject
from .ElectronMethods import ElectronMethods, ElectronEnum
from .Tools import Tools
from .ResidualGasConstantType import ResidualGasConstantType as res_gas


class Calculator(Tools, ElectronMethods):
    # TODO: Add method for calculating sigmas based on the lifetime.
    #
    # TODO: Also add method for getting the lifetime based on data from the ring (which can the use the aforementioned
    #       method)
    #       - Perhaps this method will be good for fitting the data to a lifetime (first fit, then calculate when the
    #       curve is zero):
    #          https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    beta: float = None
    Z_p: float = None
    q: float = None
    e_kin: float = None
    I_p: float = None
    n_0: float = None

    def __init__(self, data: DataObject):
        self.data = data
        self.pressure_status = data.pressure_data * 1e2  # mbar to Pa
        self.gas_fractions = data.gas_fractions
        # setting attributes from projectile_data
        for name in data.projectile_data.columns:
            if "beta" in name:
                nam = "beta"
            elif "Kin" in name:
                nam = "e_kin"
            elif "q" in name:
                nam = "q"
            elif "Z" in name:
                nam = "Z_p"
            else:
                nam = name
            setattr(self, nam, data.projectile_data[name])

    def calculate_sigma_electron_loss_parser(self) -> callable:
        sigma_el = self.get_method(ElectronEnum.DuBois_Shevelko)
        return sigma_el

    def calculate_sigma_electron_capture_parser(self):
        sigma_ec = self.get_method(ElectronEnum.Schlachter)
        return sigma_ec

    def get_all_molecular_sigmas(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        sigmas_el = {}
        sigmas_ec = {}
        for attr, Z_val in res_gas.__dict__.items():
            if attr.startswith("__"):
                continue
            sigmas_el[attr[2:]] = self.calculate_sigma_electron_loss_parser()(
                Z_val,
                self.Z_p,
                self.q,
                self.e_kin,
                self.I_p,
                self.n_0,
                self.beta,
            )
            # I haven't looked at the units, but Elias converts ekin from MeV to keV, so I'm doing the same here
            sigmas_ec[attr[2:]] = self.calculate_sigma_electron_capture_parser()(
                self.q, Z_val, self.e_kin * 1e3
            )

        sigma_molecular_el = self.get_molecular_cross_sections(sigmas_el)
        sigma_molecular_ec = self.get_molecular_cross_sections(sigmas_ec)

        return sigma_molecular_el, sigma_molecular_ec

    def calculate_full_lifetime(self) -> pd.Series:
        sigma_molecular_el, sigma_molecular_ec = self.get_all_molecular_sigmas()

        molecular_density_n = self.get_molecular_densities(
            self.gas_fractions, self.pressure_status
        )

        tau_el = self.get_molecular_lifetimes(
            molecular_density_n, sigma_molecular_el, self.beta
        )
        tau_ec = self.get_molecular_lifetimes(
            molecular_density_n, sigma_molecular_ec, self.beta
        )

        tau = 1 / ((1 / tau_el).sum() + (1 / tau_ec).sum())

        return tau

    def get_lifetime_from_data(self, injection_idx, extraction_idx, fitting_statistics=True) -> pd.Series:
        tau_series = pd.Series(dtype=float)
        for i, col in enumerate(self.data.elements.columns):
            test_element = self.data.elements[col][injection_idx:extraction_idx]
            xdata = np.array(np.arange(len(test_element)), dtype=float)
            ydata = np.array(test_element.values, dtype=float).reshape(-1)

            fun = lambda x_variable, a, b: a * np.exp(-b * x_variable)

            optimal_params, pcov = curve_fit(fun, xdata, ydata, maxfev=2000)

            tau = 1 / optimal_params[1]
            tau_series[col] = tau

            if fitting_statistics:
                self.get_lifetime_fit_statistics(
                    fun,
                    optimal_params,
                    pcov,
                    xdata,
                    ydata,
                    self.data,
                    i
                )

        return tau_series
