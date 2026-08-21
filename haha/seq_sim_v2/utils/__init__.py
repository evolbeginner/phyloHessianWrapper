"""RNG utilities: seed management and uniform/gamma random number generation.

Not thread-safe: uses a module-level numpy Generator singleton.
"""

import time
import numpy as np

_rng = np.random.default_rng()

__all__ = ["set_seed", "create_seed", "random", "gamma"]


def set_seed(seed):
    """Initialize the random number generator with a given seed."""
    global _rng
    _rng = np.random.default_rng(seed)


def create_seed():
    """Generate a seed from system time."""
    t = int(time.time() * 1e6)
    return t & 0xFFFFFFFF


def random():
    """Return a uniform random number in [0, 1)."""
    return _rng.random()


def gamma(shape):
    """Return a Gamma(shape, 1.0) random variate using the seeded RNG."""
    return _rng.gamma(shape, 1.0)
