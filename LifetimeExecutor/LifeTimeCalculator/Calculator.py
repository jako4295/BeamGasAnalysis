import pandas as pd

from .ElectronMethods import ElectronMethods, ElectronEnum
from .Tools import Tools
from .ResidualGasConstantType import ResidualGasConstantType as res_gas


class Calculator(Tools, ElectronMethods):
    # TODO: Add method for calculating sigmas based on the lifetime.
    #       Also add method for getting the lifetime based on data from the ring (which can the use the aforementioned
    #       method)
    beta: float = None
    Z_p: float = None
    q: float = None
    e_kin: float = None
    I_p: float = None
    n_0: float = None

    def __init__(
        self,
        pressure_status: float,
        gas_fractions: pd.Series,
        projectile_data: pd.DataFrame,
    ):
        self.pressure_status = pressure_status * 1e2  # mbar to Pa
        self.gas_fractions = gas_fractions
        # setting attributes from projectile_data
        for name in projectile_data.columns:
            if 'beta' in name:
                nam = 'beta'
            elif 'Kin' in name:
                nam = 'e_kin'
            elif 'q' in name:
                nam = 'q'
            elif 'Z' in name:
                nam = 'Z_p'
            else:
                nam = name
            setattr(self, nam, projectile_data[name])

    def calculate_sigma_electron_loss_parser(self) -> callable:
        sigma_el = self.get_method(ElectronEnum.Shevelko)
        return sigma_el

    def calculate_sigma_electron_capture_parser(self):
        sigma_ec = self.get_method(ElectronEnum.Schlachter)
        return sigma_ec

    @property
    def calculate_full_lifetime(self):
        sigmas_el = {}
        sigmas_ec = {}
        for attr, Z_val in res_gas.__dict__.items():
            if attr.startswith('__'):
                continue
            sigmas_el[attr[2:]] = self.calculate_sigma_electron_loss_parser()(
                self.beta,
                Z_val,
                self.I_p,
                self.n_0,
                self.q
            )
            # I haven't looked at the units, but Elias converts ekin from MeV to keV, so I'm doing the same here
            sigmas_ec[attr[2:]] = self.calculate_sigma_electron_capture_parser()(self.q, Z_val, self.e_kin * 1e3)

        sigma_molecular_el = self.get_molecular_cross_sections(sigmas_el)
        sigma_molecular_ec = self.get_molecular_cross_sections(sigmas_ec)

        molecular_density_n = self.get_molecular_densities(self.gas_fractions, self.pressure_status)

        tau_el = self.get_molecular_lifetimes(molecular_density_n, sigma_molecular_el, self.beta)
        tau_ec = self.get_molecular_lifetimes(molecular_density_n, sigma_molecular_ec, self.beta)

        tau = 1/((1/tau_el).sum() + (1/tau_ec).sum())

        return tau


