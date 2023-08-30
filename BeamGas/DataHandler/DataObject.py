from typing import Union
import pandas as pd

from .BeamDataHandler import BeamData


class DataObject(BeamData):
    def __init__(
        self,
        ring_type: str,
        gas_fraction: Union[str, pd.DataFrame],
        pressure: Union[str, pd.DataFrame],
        projectile: Union[str, pd.DataFrame],
        beam_data: str = None,
    ):
        """
        :param ring_type: str, either 'LEIR', 'PS', or 'SPS'
        :param gas_fraction: str or pd.DataFrame, path to file or dataframe containing gas fractions
        :param pressure: str or pd.DataFrame, path to file or dataframe containing pressure data
        :param projectile: str or pd.DataFrame, path to file or dataframe containing projectile data
        :param beam_data: str, path to file containing beam data. If None, no beam data will be used.
        """
        super().__init__()
        dataframes = []
        for item in [gas_fraction, pressure, projectile]:
            if type(item) == str:
                try:
                    dataframes.append(pd.read_csv(item, index_col=0))
                except FileNotFoundError:
                    raise FileNotFoundError(f"File {item} not found.")
            elif type(item) == pd.DataFrame:
                dataframes.append(item)
            else:
                raise TypeError(
                    f"Expected type str or pd.DataFrame, got {type(item)} instead."
                )
        self.gas_fractions, self.pressure_data, self.projectile_data = dataframes
        self.pressure_data = self.pressure_data.T

        self.get_ring_sorted_data(ring_type)

        if beam_data:
            self.get_data(beam_data)

    def get_ring_sorted_data(self, ring_type: str) -> None:
        col_names_ = ["Z", "I_p", "n_0"]
        col_names = [
            col_nam
            for col_nam in self.projectile_data.columns
            if ring_type in col_nam[: len(ring_type)]
        ]
        col_names = col_names_ + col_names
        self.projectile_data = self.projectile_data[col_names]
        self.projectile_data.loc[:, ("I_p", "n_0")] = self.projectile_data.loc[
            :, ("I_p", "n_0")
        ][self.projectile_data.loc[:, ("I_p", "n_0")] != 0]
        self.projectile_data = self.projectile_data.dropna()
        self.gas_fractions = self.gas_fractions[ring_type]
        self.pressure_data = self.pressure_data[ring_type][0]
