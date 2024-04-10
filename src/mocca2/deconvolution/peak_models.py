"""
Definitions of functions that model peak shape. All models should inherit from PeakModel.

See nice summary of peak models is at:
https://doi.org/10.1016/S0922-3487(98)80024-5
"""

from typing import List, Tuple
from numpy.typing import NDArray

import numpy as np


class PeakModel:
    """Abstract class for functions that model peak shape"""

    def __call__(self, t: float | NDArray, *params: float) -> float | NDArray:
        """
        Returns height of peak specified by params at given time `t`.
        """
        return self.val(t, *params)

    def val(self, t: float | NDArray, *params: float) -> float | NDArray:
        """
        Returns height of peak specified by params at given time `t`.
        """
        raise NotImplemented

    def grad(self, t: float | NDArray, *params: float) -> NDArray:
        """
        Returns gradient w.r.t params. If t is array, the return shape is [parameter, t]
        """
        raise NotImplemented

    def init_guess(
        self, height: float, maximum: float, width_left: float, width_right: float
    ) -> NDArray:
        """
        Creates initial guess of the model parameters from peak descriptors (height, location of maximum, left and right widths of peak)
        """
        raise NotImplemented

    def get_bounds(self, max_t: float) -> List[Tuple[float, float]]:
        """
        Returns bounds for individual parameters in the shape [(p1 min, p1 max), (p2 min, p2 max), ...].

        `max_t` is the maximum `t` that will be passed into the model, and also the number of time points in the peak region.
        """
        raise NotImplemented

    def n_params(self) -> int:
        """Returns number of parameters of the model"""
        raise NotImplemented


def test_gradients(model: PeakModel, t: NDArray, params: List[float]) -> bool:
    """
    Compares gradient from the model with numerical gradients. Returns True is all gradients are close.

    This is useful when making sure that gradients in new peak model are all correct.
    """
    eps = 1e-5
    grad = model.grad(t, *params)

    for pidx, pgrad in enumerate(grad):
        pl = params.copy()
        pr = params.copy()
        pl[pidx] -= eps
        pr[pidx] += eps
        num_grad = (model.val(t, *pr) - model.val(t, *pl)) / (2 * eps)

        if not np.allclose(pgrad, num_grad, rtol=1e-5):
            print("Gradient does not match for parameter", pidx + 1)
            print(pgrad)
            print(num_grad)
            return False
    return True


class BiGaussian(PeakModel):
    """BiGaussian peak model"""

    def n_params(self) -> int:
        return 4

    def val(self, t: float | NDArray, *params: float) -> float | NDArray:
        h, mu, s1, s2 = params
        v = (t < mu) * np.exp(-((t - mu) ** 2) / (2 * s1**2))
        v += (t > mu) * np.exp(-((t - mu) ** 2) / (2 * s2**2))
        return v * h

    def grad(self, t: float | NDArray, *params: float) -> NDArray:
        h, mu, s1, s2 = params

        left_gauss = (t < mu) * np.exp(-((t - mu) ** 2) / (2 * s1**2))
        right_gauss = (t > mu) * np.exp(-((t - mu) ** 2) / (2 * s2**2))

        grad_h = left_gauss + right_gauss
        grad_mu = h * (t - mu) * (left_gauss / s1**2 + right_gauss / s2**2)
        grad_s1 = h * (t - mu) ** 2 / s1**3 * left_gauss
        grad_s2 = h * (t - mu) ** 2 / s2**3 * right_gauss

        return np.array([grad_h, grad_mu, grad_s1, grad_s2])

    def init_guess(
        self, height: float, maximum: float, width_left: float, width_right: float
    ) -> NDArray:
        x0 = [height, maximum, width_left, width_right]

        return np.array(x0)

    def get_bounds(self, max_t: float) -> List[Tuple[float, float]]:
        const = [(0.0, np.inf), (0.0, max_t), (1.0, max_t / 4.0), (1.0, max_t / 4.0)]

        return const


class BiGaussianTailing(PeakModel):
    """BiGaussian with exponential tailing peak model"""

    def n_params(self) -> int:
        return 6

    def val(self, t: float | NDArray, *params: float) -> float | NDArray:
        h, mu, s1, s2, tail_h, tail_w = params
        v = (t < mu) * np.exp(-((t - mu) ** 2) / (2 * s1**2))
        v += (t > mu) * np.exp(-((t - mu) ** 2) / (2 * s2**2)) * (1.0 - tail_h)
        v += (t > mu) * np.exp(-np.clip(t - mu, 0, None) / tail_w) * tail_h
        return v * h

    def grad(self, t: float | NDArray, *params: float) -> NDArray:
        h, mu, s1, s2, tail_h, tail_w = params

        left_gauss = (t < mu) * np.exp(-((t - mu) ** 2) / (2 * s1**2))
        right_gauss = (t > mu) * np.exp(-((t - mu) ** 2) / (2 * s2**2))
        tail = (t > mu) * np.exp(-np.clip(t - mu, 0, None) / tail_w)

        grad_h = left_gauss + right_gauss * (1.0 - tail_h) + tail * tail_h
        grad_mu = (
            h * (t - mu) * (left_gauss / s1**2 + right_gauss * (1.0 - tail_h) / s2**2)
        )
        grad_mu += h * tail * tail_h / tail_w
        grad_s1 = h * (t - mu) ** 2 / s1**3 * left_gauss
        grad_s2 = h * (t - mu) ** 2 / s2**3 * right_gauss * (1.0 - tail_h)
        grad_tail_h = -h * right_gauss + h * tail
        grad_tail_w = h * tail * tail_h * (t - mu) / tail_w**2

        return np.array([grad_h, grad_mu, grad_s1, grad_s2, grad_tail_h, grad_tail_w])

    def init_guess(
        self, height: float, maximum: float, width_left: float, width_right: float
    ) -> NDArray:
        x0 = [height, maximum, width_left, width_right, 0.0, 10.0]

        return np.array(x0)

    def get_bounds(self, max_t: float) -> List[Tuple[float, float]]:
        const = [
            (0.0, np.inf),
            (0.0, max_t),
            (1.0, max_t / 4.0),
            (1.0, max_t / 4.0),
            (0.0, 0.2),
            (1.0, max_t),
        ]

        return const


class FraserSuzuki(PeakModel):
    """Fraser-Suzuki peak model"""

    def n_params(self) -> int:
        return 4

    def val(self, t: float | NDArray, *params: float) -> float | NDArray:
        h, mu, sigma, a = params

        z = 1 + a * (t - mu) / sigma
        defined = z > 0
        z = np.clip(z, 1e-100, None)

        ln_z = np.log(z)

        f = defined * np.exp(-1 / 2 / a**2 * ln_z**2)

        return f * h

    def grad(self, t: float | NDArray, *params: float) -> NDArray:
        h, mu, sigma, a = params

        z = 1 + a * (t - mu) / sigma
        defined = z > 0
        z = np.clip(z, 1e-100, None)

        ln_z = np.log(z)

        f = defined * np.exp(-1 / 2 / a**2 * ln_z**2)

        grad_h = f

        grad_mu = h * f * ln_z / (z * sigma * a)

        grad_sigma = grad_mu * (t - mu) / sigma

        grad_a = grad_mu / a**2 * (-a * (t - mu) + z * sigma * ln_z)

        return np.array([grad_h, grad_mu, grad_sigma, grad_a])

    def init_guess(
        self, height: float, maximum: float, width_left: float, width_right: float
    ) -> NDArray:
        x0 = [
            height,
            maximum,
            max((width_left + width_right) / 2, 1),
            0.1,
        ]

        return np.array(x0)

    def get_bounds(self, max_t: float) -> List[Tuple[float, float]]:
        const = [
            (0.0, np.inf),
            (0.0, max_t),
            (1.0, max_t / 4.0),
            (1e-10, 10.0),
        ]

        return const


from scipy.special import log_ndtr


def log_erfc(x):
    return log_ndtr(-x * np.sqrt(2)) + np.log(2)


class Bemg(PeakModel):
    """Bi-Exponentially Modified Gaussian peak model"""

    def n_params(self) -> int:
        return 5

    def val(self, t: float | NDArray, *params: float) -> float | NDArray:
        h, m, s, a, b = params
        tm = t - m
        as2 = a * s**2
        bs2 = b * s**2
        sq2s2 = np.sqrt(2 * s**2)

        erfc_a_arg = (as2 + tm) / sq2s2
        exp_a_arg = a * as2 / 2 + a * tm

        erfc_b_arg = (bs2 - tm) / sq2s2
        exp_b_arg = b * bs2 / 2 - b * tm

        part_a = np.exp(log_erfc(erfc_a_arg) + exp_a_arg)
        part_b = np.exp(log_erfc(erfc_b_arg) + exp_b_arg)

        scale = 0.5 * (s + 1)

        val = scale * (part_a + part_b)

        return val * h

    def grad(self, t: float | NDArray, *params: float) -> NDArray:
        h, m, s, a, b = params
        tm = t - m
        ai = a
        bi = b
        as2 = ai * s**2
        bs2 = bi * s**2
        sq2s2 = np.sqrt(2 * s**2)
        r2 = np.sqrt(2)
        rpi = np.sqrt(np.pi)

        erfc_a_arg = (as2 + tm) / sq2s2
        exp_a_arg = ai * as2 / 2 + ai * tm

        erfc_b_arg = (bs2 - tm) / sq2s2
        exp_b_arg = bi * bs2 / 2 - bi * tm

        part_a = np.exp(log_erfc(erfc_a_arg) + exp_a_arg)
        part_b = np.exp(log_erfc(erfc_b_arg) + exp_b_arg)

        scale = 0.5 * (s + 1)

        val = scale * (part_a + part_b)

        d_erfc_a = -2 / rpi * np.exp(-(erfc_a_arg**2) + exp_a_arg)
        d_erfc_b = -2 / rpi * np.exp(-(erfc_b_arg**2) + exp_b_arg)

        d_m = -d_erfc_a / sq2s2 - a * part_a
        d_m += d_erfc_b / sq2s2 + b * part_b
        d_m *= scale

        d_s = d_erfc_a * (ai - tm / s**2) / r2 + s * ai**2 * part_a
        d_s += d_erfc_b * (bi + tm / s**2) / r2 + s * bi**2 * part_b
        d_s *= scale
        d_s += 0.5 * (part_a + part_b)

        d_a = s / r2 * d_erfc_a + (as2 + tm) * part_a
        d_a *= scale

        d_b = s / r2 * d_erfc_b + (bs2 - tm) * part_b
        d_b *= scale

        return np.array([val, d_m * h, d_s * h, d_a * h, d_b * h])

    def init_guess(
        self, height: float, maximum: float, width_left: float, width_right: float
    ) -> NDArray:
        width = (width_left + width_right) / 2.0
        x0 = [height, maximum, width, 1.0, 1.0]

        return np.array(x0)

    def get_bounds(self, max_t: float) -> List[Tuple[float, float]]:
        const = [
            (0.0, np.inf),
            (0.0, max_t),
            (1.0, max_t),
            (0.01, max_t),
            (0.01, max_t),
        ]

        return const
