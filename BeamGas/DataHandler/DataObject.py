from typing import Union
import pandas as pd

from .BeamDataHandler import BeamData


class DataObject(BeamData):
    def __init__(
        self,
        ring_type: str = None,
        gas_fraction: Union[str, pd.DataFrame] = None,
        pressure: Union[str, pd.DataFrame] = None,
        projectile: Union[str, pd.DataFrame] = None,
    ):
        """
        If ring_type is given and gas_fraction, pressure, and projectile are not given, then the
        purpose would be to edit the current data. If ring_type is not given, then gas_fraction,
        pressure, and projectile must be given.

        :param ring_type: str, either 'LEIR', 'PS', or 'SPS'
        :param gas_fraction: str or pd.DataFrame, path to file or dataframe containing gas fractions
        :param pressure: str or pd.DataFrame, path to file or dataframe containing pressure data
        :param projectile: str or pd.DataFrame, path to file or dataframe containing projectile data
        :param beam_data: str, path to file containing beam data. If None, no beam data will be used.
        """
        super().__init__()
        # Check validity of input
        if all(
            [
                ring_type is None,
                gas_fraction is None,
                pressure is None,
                projectile is None,
            ]
        ):
            raise ValueError(
                "No data given. Please provide either ring_type, or gas_fraction, pressure, and projectile."
            )
        if ring_type is None and not all(
            [gas_fraction is not None, pressure is not None, projectile is not None]
        ):
            raise ValueError(
                "ring_type is None, but not all gas_fraction, pressure, projectile are given."
            )

        # Get data using the two options
        if all(
            [gas_fraction is not None, pressure is not None, projectile is not None]
        ):
            dataframes = []
            for item in [gas_fraction, pressure, projectile]:
                if type(item) == str:
                    try:
                        dataframes.append(pd.read_csv(item, index_col=0))
                    except FileNotFoundError:
                        raise FileNotFoundError(f"File {item} not found.")
                elif (
                    type(item) is pd.DataFrame
                    or type(item) is pd.Series
                    or type(item) is float
                ):
                    dataframes.append(item)
                else:
                    raise TypeError(
                        f"Expected type str or pd.DataFrame, got {type(item)} instead."
                    )

            self.gas_fractions, self.pressure_data, self.projectile_data = dataframes

            if ring_type is not None:
                if ring_type.upper() not in ("LEIR", "PS", "SPS"):
                    raise ValueError("ring_type must be either 'LEIR', 'PS', or 'SPS'.")

                self.pressure_data = self.pressure_data.T
                ring_type = ring_type.upper()
                self.get_ring_sorted_data(ring_type)

        else:
            if ring_type.upper() not in ("LEIR", "PS", "SPS"):
                raise ValueError("ring_type must be either 'LEIR', 'PS', or 'SPS'.")

            ring_type = ring_type.upper()
            self.get_ring_info(ring_type)
            self.get_ring_sorted_data(ring_type)

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
        if type(self.gas_fractions) == pd.DataFrame:
            self.gas_fractions = self.gas_fractions[ring_type]
        if type(self.pressure_data) == pd.DataFrame:
            self.pressure_data = self.pressure_data[ring_type].iloc[0]

    def get_ring_info(self, ring_type: str):
        if ring_type == "LEIR":
            gas_dict = {
                "H2": 0.83,
                "H2O": 0.02,
                "CO": 0.04,
                "CH4": 0.05,
                "CO2": 0.06,
                "He": 0.0,
                "O2": 0.0,
                "Ar": 0.0,
            }
            self.gas_fractions = pd.Series(gas_dict)
            self.pressure_data = 1e-11

        elif ring_type == "PS":
            gas_dict = {
                "H2": 0.9,
                "H2O": 0.1,
                "CO": 0.0,
                "CH4": 0.0,
                "CO2": 0.0,
                "He": 0.0,
                "O2": 0.0,
                "Ar": 0.0,
            }
            self.gas_fractions = pd.Series(gas_dict)
            self.pressure_data = 1.2e-09

        elif ring_type == "SPS":
            gas_dict = {
                "H2": 0.905,
                "H2O": 0.035,
                "CO": 0.025,
                "CH4": 0.025,
                "CO2": 0.01,
                "He": 0.0,
                "O2": 0.0,
                "Ar": 0.0,
            }
            self.gas_fractions = pd.Series(gas_dict)
            self.pressure_data = 1e-08

        data_rows = [
            [
                "He1",
                4.2,
                67.02,
                5722.74,
                0.094612,
                0.360286,
                0.99,
                1,
                1,
                2,
                2,
                0.0544177655282,
                1,
            ],
            [
                "He2",
                4.2,
                245.39,
                12279.29,
                0.094612,
                0.61131,
                0.99,
                2,
                2,
                2,
                2,
                0,
                0,
            ],
            [
                "O4",
                4.2,
                67.08,
                5723.5,
                0.094657,
                0.36059,
                0.99,
                4,
                4,
                8,
                8,
                0.113899,
                2,
            ],
            [
                "O8",
                4.2,
                245.58,
                12280.12,
                0.094657,
                0.61168,
                0.99,
                8,
                8,
                8,
                8,
                0,
                0,
            ],
            [
                "Mg6",
                4.2,
                67.09,
                5723.75,
                0.094672,
                0.360686,
                0.99,
                6,
                6,
                12,
                12,
                0.22502,
                2,
            ],
            [
                "Mg7",
                4.2,
                90.24,
                6812.68,
                0.094672,
                0.411251,
                0.99,
                7,
                7,
                12,
                12,
                0.265924,
                2,
            ],
            [
                "Pb54",
                4.2,
                72.13,
                5974.37,
                0.094647,
                0.372499,
                0.99,
                54,
                54,
                82,
                82,
                5.414,
                5,
            ],
        ]
        data_columns = [
            "Projectile",
            "LEIR_Kinj",
            "PS_Kinj",
            "SPS_Kinj",
            "LEIR_beta",
            "PS_beta",
            "SPS_beta",
            "LEIR_q",
            "PS_q",
            "SPS_q",
            "Z",
            "I_p",
            "n_0",
        ]
        self.projectile_data = pd.DataFrame(data_rows, columns=data_columns).set_index(
            "Projectile"
        )
