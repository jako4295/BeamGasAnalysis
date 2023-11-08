import pandas as pd
import os


class BeamData:
    nxcals_timestamp = None
    name = None
    elements = None

    def __init__(self):
        pass

    def get_data(self, path: str, date: str) -> None:
        """
        :param path: Path to the folder containing the data.
        :param date: Format: 'YYYY-MM-DD HH:MM-HH:MM'. Any split character
                             is allowed (also no split characters). The HH:MM
                             -HH:MM is from start to end time of the data.
        """
        data_number = "".join([n for n in date if n.isdigit()])
        element_path = None
        info_path = None
        for file in os.listdir(path):
            numbers_in_file = "".join([num for num in file if num.isdigit()])
            if numbers_in_file == data_number:
                element_path = file if "element" in file else element_path
                info_path = file if "info" in file else info_path

        if element_path is None or info_path is None:
            raise FileNotFoundError("No data found for the given date.")

        self.elements = pd.read_csv(path + element_path, index_col=0)
        self.df_info = pd.read_csv(path + info_path, index_col=0)

        self.nxcals_timestamp = self.df_info["nxcals_timestamp"]
        # self.nxcals_timestamp2 = self.df_info["nxcals_timestamp2"]
        # self.cycle_times = pd.to_datetime(self.nxcals_timestamp2) - pd.to_datetime(
        #     self.nxcals_timestamp
        # )
        self.sample_frequency = len(self.elements) / (
            self.elements.index[-1] + self.elements.index[1]
        )
        self.name = self.df_info["nxcals_variable_name"]
