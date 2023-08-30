import numpy as np
import pandas as pd
from scipy import constants


class Tools:
    @staticmethod
    def get_molecular_cross_sections(sigmas: dict) -> pd.DataFrame:
        sigmas_mol = {
            # Molecular cross sections using additivity rule from:
            # Eq (5) in https://journals.aps.org/pra/abstract/10.1103/PhysRevA.67.022706
            "H2": 2 * sigmas["H"],
            "He": sigmas["He"],
            "H2O": 2 * sigmas["H"] + sigmas["O"],
            "CO": sigmas["C"] + sigmas["O"],
            "CH4": sigmas["C"] + 4 * sigmas["H"],
            "CO2": sigmas["C"] + 2 * sigmas["O"],
            "O2": 2 * sigmas["O"],
            "Ar": sigmas["Ar"],
        }
        return pd.DataFrame(sigmas_mol)

    @staticmethod
    def get_molecular_densities(
        gas_fractions: pd.Series, pressure_status: float
    ) -> pd.Series:
        # define constants:
        K = constants.Boltzmann
        T = 298

        molecular_density_n = (gas_fractions * pressure_status) / (K * T)
        return molecular_density_n

    @staticmethod
    def get_molecular_lifetimes(
        molecular_density_n: pd.Series, sigmas_molecular: pd.DataFrame, beta: pd.Series
    ) -> pd.DataFrame:
        c_light = constants.c

        tau_df = pd.DataFrame()
        for projectile, projectile_row in pd.DataFrame(sigmas_molecular).iterrows():
            # Note that pandas automatically treats 1/0 as inf, so we don't need to worry about that
            tau_df[projectile] = 1 / (
                molecular_density_n * projectile_row * beta[projectile] * c_light
            )

        return tau_df
