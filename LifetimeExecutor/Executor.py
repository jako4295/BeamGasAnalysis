from TestSetups import Tests


def main():
    Tests.get_lifetimes_from_sigma_estimates()
    Tests.get_sigma_estimates_with_varying_i_p()
    Tests.get_sigma_estimates_with_varying_n_0()


if __name__ == "__main__":
    main()
