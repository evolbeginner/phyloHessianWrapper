"""Shared validation utilities for substitution models."""

import numpy as np


def validate_nuc_freqs(freqs):
    """Validate nucleotide frequencies: A,C,G,T order, positive, sum to 1."""
    if len(freqs) != 4:
        raise ValueError(
            f"freqs must contain exactly four values in A,C,G,T order, "
            f"got {len(freqs)}"
        )
    for i, pi in enumerate(freqs):
        if not np.isfinite(pi) or pi <= 0.0:
            raise ValueError(
                f"all base frequencies must be finite and positive, "
                f"got freq[{i}]={pi}"
            )
    total = sum(freqs)
    if not np.isclose(total, 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError(f"base frequencies must sum to 1.0, got {total}")


def validate_aa_freqs(freqs, aa_order=None):
    """Validate amino acid frequencies: positive, sum to 1."""
    n_aa = len(freqs)
    for i, pi in enumerate(freqs):
        if not np.isfinite(pi) or pi <= 0.0:
            label = f" ({aa_order[i]})" if aa_order else ""
            raise ValueError(
                f"all AA frequencies must be finite and positive, "
                f"got freq[{i}]={pi}{label}"
            )
    total = sum(freqs)
    if not np.isclose(total, 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError(f"AA frequencies must sum to 1.0, got {total}")


def validate_state(state, num_states=4):
    """Validate that state is an integer in [0, num_states)."""
    if isinstance(state, float) and not float(state).is_integer():
        raise ValueError(f"state must be an integer, got {state}")
    state_int = int(state)
    if state_int < 0 or state_int >= num_states:
        raise ValueError(
            f"state must be in [0, {num_states}), got {state_int}"
        )
    return state_int


def validate_length(length):
    """Validate that branch length is a finite non-negative real number."""
    try:
        length = float(length)
    except (TypeError, ValueError):
        raise ValueError(f"length must be a real number, got {length!r}")
    if not np.isfinite(length):
        raise ValueError(f"length must be finite, got {length}")
    if length < 0.0:
        raise ValueError(f"length must be non-negative, got {length}")
    return length


def validate_output_array(arr, expected_size):
    """Validate that output array has correct size, floating dtype, and is writable."""
    if arr.size != expected_size:
        raise ValueError(
            f"output array must contain exactly {expected_size} elements, "
            f"got {arr.size}"
        )
    if not np.issubdtype(arr.dtype, np.floating):
        raise TypeError(
            f"output array must have floating dtype, got {arr.dtype}"
        )
    if not arr.flags.writeable:
        raise ValueError("output array must be writable")


def normalize_freqs(arr):
    """Normalize a list of floats to sum to 1.0 in place."""
    s = sum(arr)
    if s == 0.0:
        raise ValueError("cannot normalize frequencies: sum is zero")
    for i in range(len(arr)):
        arr[i] /= s
