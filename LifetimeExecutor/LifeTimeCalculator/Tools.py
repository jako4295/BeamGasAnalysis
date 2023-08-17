import numpy as np
import pandas as pd
from scipy import constants, stats

from DataHandler import DataObject


class Tools:
    @staticmethod
    def get_beta(atomic_mass_in_u, Z_p, q, e_kin) -> float:
        # Mass of ion in eV
        mass_u_stripped = atomic_mass_in_u - (Z_p - q) * constants.physical_constants['electron mass in u'][0]
        mass_in_eV = mass_u_stripped * constants.physical_constants['atomic mass unit-electron volt relationship'][0]
        # Kinetic energy of ion in eV
        e_tot = mass_in_eV * 1e6 * e_kin * mass_u_stripped
        # Relativistic factors
        gamma = e_tot / mass_in_eV
    
        return np.sqrt(1 - 1 / gamma ** 2)

    @staticmethod
    def get_molecular_cross_sections(sigmas: dict) -> pd.DataFrame:
        sigmas_mol = {
            # Molecular cross sections using additivity rule from:
            # Eq (5) in https://journals.aps.org/pra/abstract/10.1103/PhysRevA.67.022706
            'H2': 2 * sigmas['H'],
            'He': sigmas['He'],
            'H2O': 2 * sigmas['H'] + sigmas['O'],
            'CO': sigmas['C'] + sigmas['O'],
            'CH4': sigmas['C'] + 4 * sigmas['H'],
            'CO2': sigmas['C'] + 2 * sigmas['O'],
            'O2': 2 * sigmas['O'],
            'Ar': sigmas['Ar'],
        }
        return pd.DataFrame(sigmas_mol)

    @staticmethod
    def get_molecular_densities(gas_fractions, pressure_status) -> pd.Series:
        # define constants:
        K = constants.Boltzmann
        T = 298

        molecular_density_n = (gas_fractions * pressure_status) / (K * T)
        return molecular_density_n

    @staticmethod
    def get_molecular_lifetimes(
            molecular_density_n: pd.Series,
            sigmas_molecular: pd.DataFrame,
            beta: pd.Series
    ) -> pd.DataFrame:
        c_light = constants.c

        tau_df = pd.DataFrame()
        for projectile, projectile_row in pd.DataFrame(sigmas_molecular).iterrows():
            # Note that pandas automatically treats 1/0 as inf, so we don't need to worry about that
            tau_df[projectile] = 1 / (molecular_density_n * projectile_row * beta[projectile] * c_light)

        return tau_df

    @staticmethod
    def get_lifetime_fit_statistics(
            fun: callable,
            optimal_params: np.ndarray,
            pcov: np.ndarray,
            xdata: np.ndarray,
            ydata: np.ndarray,
            data: DataObject,
            idx: int
    ):
        fun_optimal = lambda x_variable: fun(x_variable, *optimal_params)

        y_estimate = fun_optimal(xdata)
        res = ydata - y_estimate
        ss_res = np.sum(res**2)
        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        # for 95% confidence interval:
        # https://stats.stackexchange.com/questions/72047/when-fitting-a-curve-how-do-i-calculate-the-95-confidence-interval-for-my-fitt
        SE_a = np.sqrt(((ss_res / (len(xdata) - 2)) * pcov[0, 0]))
        t_quantile_a = stats.t(len(xdata) - 2).ppf(0.95)
        a_pm = t_quantile_a * SE_a

        SE_b = np.sqrt(((ss_res / (len(xdata) - 2)) * pcov[1, 1]))
        t_quantile_b = stats.t(len(xdata) - 2).ppf(0.95)
        b_pm = t_quantile_b * SE_b

        print("#"*47 + "\nColumn nr: " + str(idx) + "\n" + "#"*47)
        print(f"a is estimated to: {optimal_params[0]:.2e} +- {a_pm:.2e}\nb is" +
              f" estimated to: {optimal_params[1]:.2e} +- {b_pm:.2e} \nR^2={r_squared:.5f}\nLifetime is " +
              " "*10 + f"{(1/optimal_params[1])*data.sample_frequency[idx]} \nLifetime upper bound: " +
              f"{(1/(optimal_params[1] - b_pm))*data.sample_frequency[idx]} \nLifetime lower bound: " +
              f"{(1/(optimal_params[1] + b_pm))*data.sample_frequency[idx]}\n")

        # plt.plot(xdata, ydata, ".")
        # plt.plot(xdata, y_estimate)
        # lifetime_sec = ((1/optimal_params[1])*data.sample_frequency[idx]).total_seconds()
        # plt.title(f"{idx} - $R^2$: {r_squared:.3f} - lifetime: {lifetime_sec:.3f} s")
        #
        # plt.show()
