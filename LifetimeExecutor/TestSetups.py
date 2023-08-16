import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from DataHandler import DataObject
from LifeTimeCalculator import Calculator


class Tests:
    @staticmethod
    def get_lifetimes_from_sigma_estimates():
        ring_type = 'PS'
        data = DataObject(ring_type=ring_type)
        ps_object = Calculator(
            data.pressure_data,
            data.gas_fractions,
            data.projectile_data,
        )
        tau = ps_object.calculate_full_lifetime()

        fig, ax = plt.subplots()
        bar = tau.plot.bar(ax=ax, logy=True)
        ax.set_ylabel(r"Lifetime $\tau$ [s]")
        ax.set_xlabel("Projectiles")
        for i, container in enumerate(bar.containers):
            bar.bar_label(
                container,
                labels=[f"{e:,.1e}" for e in tau]
            )
        fig.tight_layout()
        fig.show()

    @staticmethod
    def get_sigma_estimates_with_varying_i_p():
        """
        This method will only vary one projectile parameter and plot the resulting sigma estimates for different I_p.

        Note: I think I_p is correct except for Pb. Reference for this is:
            https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
            First column on this page is the energy it takes to remove the first electron from the atom. This is
            reffered to as the forst ionization energy. I_p is first ionization potential of projectile.

            It may be in wrong units and I_p for Pb54 in our data is not equal to the one on the website.
        """
        ring_type = "PS"
        projectile = "O4"

        sig_el_list = {}
        sig_ec_list = {}
        # I_p_list = [0.01361806, 13.61806, 13.67]  # corresponds to [website KeV, website eV, our data]
        I_p_list = np.linspace(0.01, 20, 50)
        for I_p in I_p_list:
            data = DataObject(ring_type=ring_type)
            ps_object = Calculator(
                data.pressure_data,
                data.gas_fractions,
                data.projectile_data.loc[[projectile]],
            )
            ps_object.I_p = I_p
            sigma_el, sigma_ec = ps_object.get_all_molecular_sigmas()

            sig_el_list[I_p] = sigma_el.loc[projectile]
            sig_ec_list[I_p] = sigma_ec.loc[projectile]

        sig_el_list = pd.DataFrame(sig_el_list).T.rename_axis("I_p")
        sig_ec_list = pd.DataFrame(sig_ec_list).T.rename_axis("I_p")

        fig, ax = plt.subplots(1, 3)
        chosen_gasses = ['H2', 'Ar', 'CO2']
        for i, gas in enumerate(chosen_gasses):
            ax[i].plot(sig_el_list[gas], label="EL")
            ax[i].plot(sig_ec_list[gas], label="EC")
            ax[i].set_title(gas)
            ax[i].set_xlabel(r"$I_p$ [eV]")
            if i == 0:
                ax[i].set_ylabel(r"$\sigma$ [m$^2$]")
            #ax[i].set_xscale('log')
            ax[i].set_yscale('log')
            ax[i].legend()
        fig.tight_layout()
        fig.show()

    @staticmethod
    def get_sigma_estimates_with_varying_n_0():
        ring_type = "PS"
        projectile = "O4"

        sig_el_list = {}
        sig_ec_list = {}
        n_0_list = np.arange(1, 11)  # corresponds to [website KeV, website eV, our data]
        for n_0 in n_0_list:
            data = DataObject(ring_type=ring_type)
            ps_object = Calculator(
                data.pressure_data,
                data.gas_fractions,
                data.projectile_data.loc[[projectile]],
            )
            ps_object.n_0 = n_0
            sigma_el, sigma_ec = ps_object.get_all_molecular_sigmas()

            sig_el_list[n_0] = sigma_el.loc[projectile]
            sig_ec_list[n_0] = sigma_ec.loc[projectile]

        sig_el_list = pd.DataFrame(sig_el_list).T.rename_axis("I_p")
        sig_ec_list = pd.DataFrame(sig_ec_list).T.rename_axis("I_p")

        fig, ax = plt.subplots(1, 3)
        chosen_gasses = ['H2', 'Ar', 'CO2']
        for i, gas in enumerate(chosen_gasses):
            ax[i].plot(sig_el_list[gas].values, '.', label="EL")
            ax[i].plot(sig_ec_list[gas].values, '.', label="EC")
            ax[i].set_title(gas)
            ax[i].set_xlabel(r"$I_p$ [eV]")
            if i == 0:
                ax[i].set_ylabel(r"$\sigma$ [m$^2$]")
            #ax[i].set_xscale('log')
            ax[i].set_yscale('log')
            ax[i].legend()
        fig.tight_layout()
        fig.show()
