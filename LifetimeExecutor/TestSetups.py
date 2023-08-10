from matplotlib import pyplot as plt

from DataHandler import DataObject
from LifeTimeCalculator import Calculator


class Tests:
    @staticmethod
    def get_lifetimes_from_sigma_estimates():
        ring_type = 'SPS'
        data = DataObject(ring_type=ring_type)
        ps_object = Calculator(
            data.pressure_data,
            data.gas_fractions,
            data.projectile_data,
        )
        tau = ps_object.calculate_full_lifetime

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
