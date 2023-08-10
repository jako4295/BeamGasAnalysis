from enum import Enum

import numpy as np
from scipy.constants import c as c_light

import pandas as pd


class ElectronEnum(Enum):
    Schlachter = 0
    Shevelko = 1


class ElectronMethods:
    @classmethod
    def get_method(cls, method: ElectronEnum) -> callable:
        match method:
            case ElectronEnum.Schlachter:
                return cls.schlachter
            case ElectronEnum.Shevelko:
                return cls.shevelko
            case _:
                raise NotImplementedError(f"Method {method} not implemented")

    @classmethod
    def schlachter(cls, q: pd.Series, Z: pd.Series, E_kin: pd.Series) -> pd.Series:
        """
        EC
        Calculates the electron capture cross section based on the Schlachter formula:
        (3), (5) in: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.3372
        """

        E_tilde = E_kin / ((Z**1.25)*(q**0.7))
        sigma_tilde = ((1.1 * (10 ** (-8))) / (E_tilde ** 4.8)) * (1 - np.exp(-0.037 * E_tilde ** 2.2)) * (
                    1 - np.exp(-2.44 * (10 ** (-5)) * E_tilde ** 2.6))

        return (sigma_tilde * (q ** 0.5)) / (Z ** 1.8)

    @classmethod
    def shevelko(
            cls,
            beta: pd.Series,
            Z: pd.Series,
            I_p: pd.Series,
            n_0: pd.Series,
            q: pd.Series,
    ) -> pd.Series:
        """
        EL
        Calculates the electron loss cross section based on the Shevelko formula:
        (12)-(13) in: https://www.sciencedirect.com/science/article/pii/S0168583X11003272
        """
        Ry = 13.606 / 1e3  # Rydberg energy in keV
        u = (beta ** 2 * c_light ** 2 * Ry) / I_p
        return 0.88 * (10 ** (-16)) * ((Z + 1) ** 2) * (u / (u ** 2 + 3.5)) * (Ry / I_p) ** (1 + 0.01 * q) * (
                4 + (1.31 / n_0) * np.log(4 * u + 1))
