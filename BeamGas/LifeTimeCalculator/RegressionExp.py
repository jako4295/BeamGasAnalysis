from typing import Any
import matplotlib

import numpy as np
from matplotlib import transforms, pyplot as plt
from matplotlib.patches import Ellipse
from numpy import ndarray
import pandas as pd
from scipy.optimize import newton
from scipy import stats

from BeamGas.DataHandler import DataObject


class ExpRegression:
    # Currently this method is not used, but it is kept for future reference.
    # code and derivation from:
    # https://scipython.com/blog/least-squares-fitting-to-an-exponential-function/#comments
    optimal_parameters: tuple[float, float] = None

    def __init__(self, x: np.ndarray, y: np.ndarray):
        a, b = self.log_transformed_fit(
            x, y
        )  # used for initial guess for least squares problem
        a1, b1 = self.nonlinear_one_dimension_fit(x, y, (a, b))
        self.optimal_parameters = (a1, -b1)

    def log_transformed_fit(self, x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
        """Ordinary linear least-squares fit to ln(y) = ln(a) + bx."""
        n = len(x)
        lny = np.log(y)
        Sx = np.sum(x)
        Sy = np.sum(lny)
        Sxx = np.sum(x**2)
        Sxy = np.sum(x * lny)
        den = n * Sxx - Sx**2
        a = (Sy * Sxx - Sxy * Sx) / den
        b = (n * Sxy - Sy * Sx) / den
        return np.exp(a), b

    def nonlinear_one_dimension_fit(
        self, x: np.ndarray, y: np.ndarray, prms0: tuple[float, float]
    ) -> tuple[float, float]:
        """Indirect nonlinear fit to y = a.exp(bx), treating a = a(b)."""

        b0 = prms0[1]
        # Use Newton-Raphson to find the root of dr2/db.
        b = newton(self._dr2db, b0, args=(x, y))
        a, _ = self._get_a_and_dadb(b, x, y)
        return a, b

    def _dr2db(self, b: float, x: np.ndarray, y: np.ndarray) -> float:
        a, dadb = self._get_a_and_dadb(b, x, y)
        fac1 = y - a * np.exp(b * x)
        fac2 = dadb + a * x
        return -2 * np.sum(fac1 * np.exp(b * x) * fac2)

    def _get_a_and_dadb(
        self, b: float, x: np.ndarray, y: np.ndarray
    ) -> tuple[float, float]:
        S1 = np.sum(y * np.exp(b * x))
        S2 = np.sum(np.exp(2 * b * x))
        dS1db = np.sum(x * y * np.exp(b * x))
        dS2db = np.sum(2 * x * np.exp(2 * b * x))
        a = S1 / S2
        dadb = (S2 * dS1db - S1 * dS2db) / S2**2
        return a, dadb


class StatisticalSummary:
    @classmethod
    def plot_confidence_ellipse(
        cls, a_arr: np.ndarray, tau_arr: np.ndarray
    ) -> tuple[ndarray, float | Any]:
        """
        Plots a confidence ellipse of the two parameters from the exponential
        fit given as: I(t) = a*exp(-t/tau), where I(t) is the beam intensity
        at time t, a is the initial intensity and tau is the lifetime.

        The confidence ellipse is based on the standard deviation (n_std) from
        the cls._confidence_ellipse function. Remember 2 std is approximately
        a 95% confidence interval.

        :param x: np.ndarray. The x-axis values for the data.
        :param a_arr: np.ndarray. The a values from the exponential fit.
        :param tau_arr: np.ndarray. The tau values from the exponential fit.
        """

        fig, ax = plt.subplots()
        for a_, b_ in zip(a_arr, tau_arr):
            ax.scatter(a_, b_, c="b", marker="x")
        ax_patch, b_center, b_deviation = cls._confidence_ellipse(
            a_arr, tau_arr, ax, 2, edgecolor="b"
        )
        ax.hlines(
            b_center - b_deviation,
            xmin=np.min(a_arr),
            xmax=np.max(a_arr),
            color="gray",
            alpha=0.5,
        )
        ax.hlines(
            b_center + b_deviation,
            xmin=np.min(a_arr),
            xmax=np.max(a_arr),
            color="gray",
            alpha=0.5,
        )
        ax.set_xlabel("a")
        ax.set_ylabel(r"$\tau$")
        fig.show()

        print(
            f"With 95.5% confidence b is in the interval: {b_center - b_deviation} to {b_center + b_deviation}."
            + f"\nThis is based on {len(tau_arr)} cycles of the experiment."
        )

        return b_center, b_deviation

    @classmethod
    def _confidence_ellipse(
        cls,
        x: np.ndarray,
        y: np.ndarray,
        ax: plt.Axes,
        n_std: int = 2,
        facecolor: str = "none",
        **kwargs,
    ) -> tuple[matplotlib.patches.Ellipse, float, float]:
        """
        Credit: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

        Create a plot of the covariance confidence ellipse of *x* and *y*.
        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.
        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.
        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.
        **kwargs
            Forwarded to `~matplotlib.patches.Ellipse`
        Returns
        -------
        matplotlib.patches.Ellipse, center in Y axis, deviation in Y axis
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")
        cov = np.cov(x, y)
        pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensional dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse(
            (0, 0),
            width=ell_radius_x * 2,
            height=ell_radius_y * 2,
            facecolor=facecolor,
            **kwargs,
        )
        # Calculating the standard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)
        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)
        transf = (
            transforms.Affine2D()
            .rotate_deg(45)
            .scale(scale_x, scale_y)
            .translate(mean_x, mean_y)
        )
        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse), mean_y, scale_y

    @staticmethod
    def get_lifetime_fit_statistics(
        fun: callable,
        optimal_params: np.ndarray,
        pcov: np.ndarray,
        xdata: np.ndarray,
        ydata: np.ndarray,
        data: DataObject,
        idx: int,
    ) -> tuple[float, float]:
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

        SE_tau = np.sqrt(((ss_res / (len(xdata) - 2)) * pcov[1, 1]))
        t_quantile_tau = stats.t(len(xdata) - 2).ppf(0.95)
        tau_pm = t_quantile_tau * SE_tau

        print("#" * 47 + "\nColumn nr: " + str(idx) + "\n" + "#" * 47)
        print(
            f"a is estimated to: {optimal_params[0]:.2e} +- {a_pm:.2e}\nb is"
            + f" estimated to: {optimal_params[1]:.2e} +- {tau_pm:.2e} \nR^2={r_squared:.5f}\nLifetime is "
            + " " * 10
            + f"{optimal_params[1]*(1 / data.sample_frequency)} \nLifetime upper bound: "
            + f"{(optimal_params[1] - tau_pm)*(1 / data.sample_frequency)} \nLifetime lower bound: "
            + f"{(optimal_params[1] + tau_pm)*(1 / data.sample_frequency)}\n"
        )

        return a_pm, tau_pm

    @staticmethod
    def plot_exponential_regression_summary(
        a: np.ndarray, tau: np.ndarray, a_pm, tau_pm, elements: pd.Series
    ):
        x_axis = np.arange(len(elements))
        fits_upper = np.array(
            list(
                map(
                    lambda a_, tau_: a_ * np.exp(-x_axis / tau_),
                    list(a + a_pm),
                    list(tau + tau_pm),
                )
            )
        )
        fits_lower = np.array(
            list(
                map(
                    lambda a_, tau_: a_ * np.exp(-x_axis / tau_),
                    list(a - a_pm),
                    list(tau - tau_pm),
                )
            )
        )
        fits = np.array(
            list(map(lambda a_, tau_: a_ * np.exp(-x_axis / tau_), list(a), list(tau)))
        )
        fig, ax = plt.subplots()
        elements.plot(
            ax=ax, style=".", markersize=0.5, alpha=0.15, color="red", legend=False
        )
        for fit_l, fit_u, fit in zip(fits_lower, fits_upper, fits):
            ax.plot(x_axis, fit_l, color="gray", alpha=0.05)
            ax.plot(x_axis, fit_u, color="gray", alpha=0.05)
            ax.plot(x_axis, fit, color="k", alpha=0.05)
        ax.set_xlabel("Time [ms]")
        ax.set_ylabel("Intensity (x 1e10 charges)")
        fig.show()
