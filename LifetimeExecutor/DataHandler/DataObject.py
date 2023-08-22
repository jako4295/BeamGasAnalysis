import pandas as pd

from .BeamDataHandler import BeamData


class DataObject(BeamData):
    def __init__(self, ring_type: str):
        super().__init__()
        self.gas_fractions = pd.read_csv("DataHandler/Gas_fractions.csv", index_col=0)
        self.pressure_data = pd.read_csv("DataHandler/Pressure_data.csv", index_col=0).T
        self.projectile_data = pd.read_csv("DataHandler/Projectile_data2.csv", index_col=0)
        self.get_ring_sorted_data(ring_type)

    def get_ring_sorted_data(self, ring_type: str):
        col_names_ = ['Z', 'I_p', 'n_0']
        col_names = [col_nam for col_nam in self.projectile_data.columns if ring_type in col_nam[:len(ring_type)]]
        col_names = col_names_ + col_names
        self.projectile_data = self.projectile_data[col_names]
        self.projectile_data[["I_p", "n_0"]] = self.projectile_data[["I_p", "n_0"]][
            self.projectile_data[["I_p", "n_0"]] != 0
        ]
        self.projectile_data = self.projectile_data.dropna()
        self.gas_fractions = self.gas_fractions[ring_type]
        self.pressure_data = self.pressure_data[ring_type][0]
