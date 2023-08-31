import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import constants
from typing import Union

from BeamGas.DataHandler import DataObject
from .ElectronMethods import ElectronMethods, ElectronEnum
from .Tools import Tools
from .ResidualGasConstantType import ResidualGasConstantType as res_gas
from .RegressionExp import ExpRegression, StatisticalSummary


class Calculator(Tools, ElectronMethods):
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

    def calculate_sigma_electron_capture_parser(self) -> callable:
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

    def get_lifetime_from_data(
        self,
        injection_idx: int = None,
        extraction_idx: int = None,
        fitting_statistics=True,
        return_lifetime_band=True,
    ) -> pd.Series:
        if injection_idx is None:
            injection_idx = 0
        if extraction_idx is None:
            extraction_idx = len(self.data.elements) + 1
        tau_series = pd.Series(dtype=float)
        a_arr = np.zeros(len(self.data.elements.columns))
        tau_arr = np.zeros(len(self.data.elements.columns))
        a_pm_arr = np.zeros(len(self.data.elements.columns))
        tau_pm_arr = np.zeros(len(self.data.elements.columns))
        for i, col in enumerate(self.data.elements.columns):
            test_element = self.data.elements[col][injection_idx:extraction_idx]
            xdata = np.array(np.arange(len(test_element)), dtype=float)
            ydata = np.array(test_element.values, dtype=float).reshape(-1)
            ydata_smooth = pd.Series(ydata).rolling(10, win_type="triang").mean().values

            fun = lambda x_variable, a, b: a * np.exp(-x_variable / b)

            optimal_params, pcov = curve_fit(fun, xdata, ydata, maxfev=2000)
            # optimal_params = ExpRegression(xdata, ydata).optimal_parameters

            # tau = 1 / optimal_params[1]
            tau = optimal_params[1]
            tau_series[col] = tau

            a_arr[i], tau_arr[i] = optimal_params
            if fitting_statistics:
                a_pm, tau_pm = StatisticalSummary.get_lifetime_fit_statistics(
                    fun, optimal_params, pcov, xdata, ydata, self.data, i
                )
                a_pm_arr[i], tau_pm_arr[i] = a_pm, tau_pm

        if return_lifetime_band:
            tau_center, tau_deviation = StatisticalSummary.plot_confidence_ellipse(
                a_arr, tau_arr
            )
            StatisticalSummary.plot_exponential_regression_summary(
                a_arr,
                tau_arr,
                a_pm_arr,
                tau_pm_arr,
                self.data.elements.iloc[injection_idx:extraction_idx].reset_index(
                    drop=True
                ),
            )
            tau_series["tau_lower"] = tau_center + tau_deviation
            tau_series["tau_upper"] = tau_center - tau_deviation

        return tau_series

    def get_sigma_from_lifetime(
        self, tau: pd.Series, projectile: str = "Pb54"
    ) -> Union[pd.Series, pd.DataFrame]:
        molecular_density_n = self.get_molecular_densities(
            self.data.gas_fractions, self.data.pressure_data
        )
        beta = self.beta
        if ("tau_upper" and "tau_lower") in tau.index:
            sigma_from_tau = pd.DataFrame(
                {
                    name: 1
                    / (molecular_density_n * tau_val * beta[projectile] * constants.c)
                    for name, tau_val in tau.iloc[-2:].items()
                }
            )
        else:
            sigma_from_tau = 1 / (
                molecular_density_n * tau.mean() * beta[projectile] * constants.c
            )
        return sigma_from_tau
