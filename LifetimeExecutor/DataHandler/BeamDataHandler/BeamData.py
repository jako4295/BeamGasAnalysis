import pandas as pd


class BeamData:
    nxcals_timestamp = None
    name = None
    elements = None

    def __init__(self):
        path = "DataHandler/BeamDataHandler/"
        self.elements = pd.read_csv(path+'df_elements1-long.csv', index_col=0)
        self.df_info = pd.read_csv(path+'df_info1-long.csv', index_col=0)

        self.nxcals_timestamp = self.df_info["nxcals_timestamp"]
        self.nxcals_timestamp2 = self.df_info["nxcals_timestamp2"]
        self.cycle_times = pd.to_datetime(self.nxcals_timestamp2) - pd.to_datetime(
            self.nxcals_timestamp
        )
        self.sample_frequency = self.cycle_times / len(self.elements)
        self.name = self.df_info["nxcals_variable_name"]
