from .Plots import Plots


def main():
    Plots.get_lifetimes_from_sigma_estimates()
    Plots.get_sigma_estimates_with_varying_i_p()
    Plots.get_sigma_estimates_with_varying_n_0()
    Plots.get_lifetime_from_data()


if __name__ == "__main__":
    main()
