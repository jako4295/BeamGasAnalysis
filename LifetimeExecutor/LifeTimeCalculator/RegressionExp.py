from typing import Any

import numpy as np
from matplotlib import transforms, pyplot as plt
from matplotlib.patches import Ellipse
from numpy import ndarray
from scipy.optimize import newton
from scipy import stats

from DataHandler import DataObject


class ExpRegression:
    # code and derivation from:
    # https://scipython.com/blog/least-squares-fitting-to-an-exponential-function/#comments
    optimal_parameters: tuple[float, float] = None

    def __init__(self, x: np.ndarray, y: np.ndarray):
        a, b = self.log_transformed_fit(
            x, y
        )  # used for initial guess for least squares problem
        a1, b1 = self.nonlinear_one_dimension_fit(x, y, (a, b))
        self.optimal_parameters = (a1, -b1)

    def log_transformed_fit(self, x: np.ndarray, y: np.ndarray):
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

    def nonlinear_one_dimension_fit(self, x, y, prms0):
        """Indirect nonlinear fit to y = a.exp(bx), treating a = a(b)."""

        b0 = prms0[1]
        # Use Newton-Raphson to find the root of dr2/db.
        b = newton(self._dr2db, b0, args=(x, y))
        a, _ = self._get_a_and_dadb(b, x, y)
        return a, b

    def _dr2db(self, b, x, y):
        a, dadb = self._get_a_and_dadb(b, x, y)
        fac1 = y - a * np.exp(b * x)
        fac2 = dadb + a * x
        return -2 * np.sum(fac1 * np.exp(b * x) * fac2)

    def _get_a_and_dadb(self, b, x, y):
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
        cls, x: np.ndarray, a_arr: np.ndarray, b_arr: np.ndarray
    ) -> tuple[ndarray, float | Any]:
        # The confidence region size is based on the standard deviation (n_std) from the cls.confidence_ellipse function
        # Remember 2 std = 95% confidence interval
        fig, ax = plt.subplots()
        da = 0.1
        A, B = np.meshgrid(
            np.linspace(np.min(a_arr) * (1 - da), np.max(a_arr) * (1 + da), 50),
            np.linspace(np.min(b_arr) * (1 - da), np.max(b_arr) * (1 + da), 50),
        )
        y = np.mean(a_arr) * np.exp(-x / np.mean(b_arr))
        r2 = np.sum((y - A[:, :, None] * np.exp(-x / B[:, :, None])) ** 2, axis=2)
        r2 = np.log(r2)
        ax.contourf(A, B, r2, alpha=0.5)
        for a_, b_ in zip(a_arr, b_arr):
            ax.scatter(a_, b_, c="b", marker="x")
        ax_patch, b_center, b_deviation = cls.confidence_ellipse(
            a_arr, b_arr, ax, 2, edgecolor="b"
        )
        ax.hlines(
            b_center - b_deviation,
            xmin=np.min(A[:, :, None]),
            xmax=np.max(A[:, :, None]),
            color="gray",
            alpha=0.5,
        )
        ax.hlines(
            b_center + b_deviation,
            xmin=np.min(A[:, :, None]),
            xmax=np.max(A[:, :, None]),
            color="gray",
            alpha=0.5,
        )
        ax.set_xlabel("a")
        ax.set_ylabel(r"$\tau$")
        fig.show()

        print(
            f"With 95.5% confidence b is in the interval: {b_center - b_deviation} to {b_center + b_deviation}."
            + f"\nThis is based on {len(b_arr)} cycles of the experiment."
        )

        return b_center, b_deviation

    @classmethod
    def confidence_ellipse(cls, x, y, ax, n_std=3.0, facecolor="none", **kwargs):
        """
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
    ):
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

        SE_b = np.sqrt(((ss_res / (len(xdata) - 2)) * pcov[1, 1]))
        t_quantile_b = stats.t(len(xdata) - 2).ppf(0.95)
        b_pm = t_quantile_b * SE_b

        print("#" * 47 + "\nColumn nr: " + str(idx) + "\n" + "#" * 47)
        print(
            f"a is estimated to: {optimal_params[0]:.2e} +- {a_pm:.2e}\nb is"
            + f" estimated to: {optimal_params[1]:.2e} +- {b_pm:.2e} \nR^2={r_squared:.5f}\nLifetime is "
            + " " * 10
            + f"{optimal_params[1]*data.sample_frequency[idx]} \nLifetime upper bound: "
            + f"{(optimal_params[1] - b_pm)*data.sample_frequency[idx]} \nLifetime lower bound: "
            + f"{(optimal_params[1] + b_pm)*data.sample_frequency[idx]}\n"
        )

        # plt.plot(xdata, ydata, ".")
        # plt.plot(xdata, y_estimate)
        # lifetime_sec = ((1/optimal_params[1])*data.sample_frequency[idx]).total_seconds()
        # plt.title(f"{idx} - $R^2$: {r_squared:.3f} - lifetime: {lifetime_sec:.3f} s")
        #
        # plt.show()
