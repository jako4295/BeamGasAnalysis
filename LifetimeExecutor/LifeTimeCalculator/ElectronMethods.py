from enum import Enum

import numpy as np
from scipy.constants import c as c_light

import pandas as pd


class ElectronEnum(Enum):
    Schlachter = 0
    Shevelko = 1
    DuBois_Shevelko = 2


class ElectronMethods:
    @classmethod
    def get_method(cls, method: ElectronEnum) -> callable:
        match method:
            case ElectronEnum.Schlachter:
                return cls.schlachter
            case ElectronEnum.Shevelko:
                return cls.shevelko
            case ElectronEnum.DuBois_Shevelko:
                return cls.dubois_shevelko
            case _:
                raise NotImplementedError(f"Method {method} not implemented")

    @classmethod
    def dubois_shevelko(
        cls,
        Z: int,
        Z_p: pd.Series,
        q: pd.Series,
        E_kin: pd.Series,
        I_p: pd.Series,
        n_0: pd.Series,
        beta: pd.Series,
    ) -> pd.Series:
        par = np.array(
            [
                10.88,
                0.95,
                2.5,
                1.1137,
                -0.1805,
                2.64886,
                1.35832,
                0.80696,
                1.00514,
                6.13667,
            ]
        )
        dubois_shevelko_comb = cls._dubois_comb(
            Z, E_kin, I_p, beta, par
        ) * cls._shevelko_comb(Z_p, q, E_kin, I_p, n_0, beta, par)
        dubois_shevelko_comb *= 1e-4  # conversion to SI units

        dubois_shevelko_comb[
            q == Z_p
        ] = 0  # Fully stripped ion (i.e. q = Z_p) cannot lose electrons
        return dubois_shevelko_comb

    @classmethod
    def _dubois_comb(
        cls, Z: int, E_kin: pd.Series, I_p: pd.Series, beta: pd.Series, par: np.ndarray
    ) -> pd.Series:
        alpha = 7.2973525376 * (10 ** (-3))
        AU = 931.5016  # MeV
        m_e = 510.998  # keV
        g = 1 / (np.sqrt(1 - beta**2))  # Lorentz factor

        N_eff = min(10 ** (par[3] * np.log10(Z) + par[4] * np.log10(Z) ** 2), Z)
        F1 = (
            (N_eff + par[0] * (Z - N_eff) * (g - 1))
            .abs()
            .apply(lambda x: x if x < Z else Z)
        )  # min(Z, F1) for Series
        F2 = Z * (1 - np.sqrt(1 - F1 / Z))
        F3 = (-(F2 ** par[1])) / (
            (np.sqrt(2 * E_kin / AU) + np.sqrt(2 * I_p / m_e)) / alpha
        )
        return F1 + (F2 * np.exp(F3)) ** 2

    @classmethod
    def _shevelko_comb(
        cls,
        Z_p: pd.Series,
        q: pd.Series,
        E_kin: pd.Series,
        I_p: pd.Series,
        n_0: pd.Series,
        beta: pd.Series,
        par: np.ndarray,
    ) -> pd.Series:
        alpha = 7.2973525376 * (10 ** (-3))
        Ry = 13.606 / 1e3  # Rydberg energy in keV
        g = 1 / (np.sqrt(1 - beta**2))  # Lorentz factor

        # AU = 931.5016  # MeV
        # g = 1.0 + E_kin / AU  # gamma factor
        # beta = np.sqrt(1.0 - 1.0 / g**2)  # beta factor
        u = ((beta / alpha) ** 2) / (I_p / Ry)
        return (
            ((par[5] * (10 ** (-16)) * u) / (u**2 + par[6]))
            * (Ry / I_p)
            ** (par[7] + ((q + par[2]) / Z_p) * (1 - np.exp(-par[9] * (u**2))))
            * (1 + (par[8] / n_0) * np.log((u + 1) * g))
        )

    @classmethod
    def schlachter(cls, q: pd.Series, Z: int, E_kin: pd.Series) -> pd.Series:
        """
        EC
        Calculates the electron capture cross section based on the Schlachter formula:
        (3), (5) in: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.3372
        """

        E_tilde = E_kin / ((Z**1.25) * (q**0.7))
        sigma_tilde = (
            ((1.1 * (10 ** (-8))) / (E_tilde**4.8))
            * (1 - np.exp(-0.037 * E_tilde**2.2))
            * (1 - np.exp(-2.44 * (10 ** (-5)) * E_tilde**2.6))
        )
        sigma = (sigma_tilde * (q**0.5)) / (Z**1.8)

        return sigma * 1e-4  # conversion to SI units (m^2)

    @classmethod
    def shevelko(
        cls,
        beta: pd.Series,
        Z: int,
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
        u = (beta**2 * c_light**2 * Ry) / I_p
        return (
            0.88
            * (10 ** (-16))
            * ((Z + 1) ** 2)
            * (u / (u**2 + 3.5))
            * (Ry / I_p) ** (1 + 0.01 * q)
            * (4 + (1.31 / n_0) * np.log(4 * u + 1))
        )
