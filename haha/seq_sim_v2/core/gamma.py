"""Gamma distribution: random variate generation and discrete rate categories.

Uses scipy for probability functions (gammaln, gammainc, gamma.ppf, chi2.ppf, norm.ppf)
and the project's seeded numpy RNG for uniform/gamma random variates.
"""

import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma as _gamma_dist
from ..utils import gamma as _rng_gamma


def rndgamma(shape):
    """Generate a Gamma(shape, 1) random variate using the seeded RNG."""
    return _rng_gamma(shape)


def _gamma_ppf(prob, alpha, beta):
    """Gamma quantile: P(X <= x) = prob for Gamma(alpha, rate=beta).

    Uses scipy scaling convention: Gamma(a, scale=1/beta).
    """
    return _gamma_dist.ppf(prob, alpha, scale=1.0 / beta)


def discrete_gamma(alpha, k, median=False):
    """Return k equiprobable rate categories from Gamma(alpha, alpha) with mean 1.

    Returns (freqK, rK) where freqK are equiprobable bin probabilities (all 1/k)
    and rK are the representative rate values for each category.

    Three methods, applied in order of preference:
      1. mean  — integrate conditional means per bin via incomplete gamma
      2. median — use conditional medians per bin
      3. equal — fallback with zero-handling for extreme parameters

    Parameters
    ----------
    alpha : float
        Gamma shape parameter (also acts as rate, so mean = alpha/alpha = 1).
    k : int
        Number of discrete categories (2–32).
    median : bool
        If True, use median method directly instead of mean.
    """
    beta = alpha
    factor = alpha / beta * k  # = k
    rK = [0.0] * k
    freqK = [0.0] * k

    def _is_valid():
        return all(v >= 0 and np.isfinite(v) for v in rK) and any(v > 1e-12 for v in rK)

    def _compute_mean():
        """Integrate per-bin means:
            rK[i] = k * [I_avg(i+1) - I_avg(i)]
        where I_avg(i) = P(X < q(i), shape=alpha+1) is computed via
        the regularized incomplete gamma."""
        for i in range(k - 1):
            q = _gamma_ppf((i + 1.0) / k, alpha, beta)
            v = gammainc(alpha + 1.0, q * beta)
            if not np.isfinite(v) or v < 0:
                raise ValueError
            freqK[i] = v
        rK[0] = freqK[0] * factor
        rK[k - 1] = (1.0 - freqK[k - 2]) * factor
        for i in range(1, k - 1):
            rK[i] = (freqK[i] - freqK[i - 1]) * factor
        t = sum(rK)
        for i in range(k):
            rK[i] *= factor / t

    def _compute_median():
        """Use conditional medians of equiprobable intervals."""
        step = 1.0 / (2.0 * k)
        for i in range(k):
            rK[i] = _gamma_ppf((2 * i + 1) * step, alpha, beta)
        t = sum(rK)
        for i in range(k):
            rK[i] *= factor / t

    def _fallback_equal():
        """Equal-spacing with zero-value replacement for extreme parameters."""
        step = 1.0 / (2.0 * k)
        for i in range(k):
            v = _gamma_ppf((2 * i + 1) * step, alpha, beta)
            rK[i] = v if (v >= 0 and np.isfinite(v)) else 0.0
        positives = [v for v in rK if v > 0]
        if positives:
            min_pos = min(positives)
            for i in range(k):
                if rK[i] == 0.0:
                    rK[i] = min_pos * 0.5
        t = sum(rK)
        for i in range(k):
            rK[i] *= factor / t

    try:
        if median:
            _compute_median()
        else:
            _compute_mean()
        if not _is_valid():
            raise ValueError("rates invalid after primary method")
    except (ValueError, ZeroDivisionError, OverflowError):
        try:
            _compute_median()
            if not _is_valid():
                raise ValueError("rates invalid after median fallback")
        except (ValueError, ZeroDivisionError, OverflowError):
            _fallback_equal()

    for i in range(k):
        freqK[i] = 1.0 / k
    return freqK, rK
