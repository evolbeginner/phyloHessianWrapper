"""Model registry: dispatch substitution model by name and return configured instance."""

import numpy as np
from ._validators import normalize_freqs
from .nuc_models import HKYModel, F84Model, GTRModel, JC69Model
from .aa_models import create_aa_model

NUC_MODELS = {"HKY", "F84", "GTR", "REV", "JC69"}
AA_MODELS = {"JTT", "WAG", "PAM", "BLOSUM", "LG", "GENERAL"}


def setup_model(name, tstv=None, equal_exchangeability=True, freqs=None, equal_freqs=True,
                relative_rates=None):
    """Create and return a configured substitution model instance."""
    name = name.upper()

    if name == "REV":
        name = "GTR"

    if name in NUC_MODELS:
        if freqs is None:
            freqs = [0.25, 0.25, 0.25, 0.25]
        else:
            freqs = list(freqs)
            normalize_freqs(freqs)

        if name == "JC69":
            return JC69Model()
        elif name == "HKY":
            return HKYModel(freqs, tstv, equal_exchangeability)
        elif name == "F84":
            return F84Model(freqs, tstv, equal_exchangeability)
        elif name == "GTR":
            if relative_rates is None:
                relative_rates = [1.0] * 6
            return GTRModel(freqs, relative_rates)

    elif name in AA_MODELS:
        return create_aa_model(name, user_freqs=freqs, user_rates=relative_rates,
                               equal_freqs=equal_freqs)
    else:
        raise ValueError(f"Unknown model: {name}")


def build_model_registry(name, freqs, tstv, equal_exchangeability, equal_freqs,
                          relative_rates, alt_freqs, alt_tstv, alt_rates):
    """Build a dict mapping (comp_idx, rmat_idx) to model instances for node-heterogeneous simulation."""
    name = name.upper()
    if name == "REV":
        name = "GTR"

    all_freqs = [freqs]
    if alt_freqs:
        all_freqs.extend(alt_freqs)

    all_rates = [relative_rates]
    if alt_rates:
        all_rates.extend(alt_rates)

    all_tstv_vals = [tstv]
    if alt_tstv:
        all_tstv_vals.extend(alt_tstv)

    registry = {}
    n_comps = len(all_freqs)
    n_rmats = max(len(all_rates), len(all_tstv_vals), 1)

    if name == "JC69":
        for c in range(n_comps):
            for r in range(n_rmats):
                registry[(c, r)] = JC69Model()
    elif name in ("HKY", "F84"):
        for c in range(n_comps):
            f = all_freqs[c] if c < len(all_freqs) else None
            eq_f = equal_freqs if c == 0 else (f is None)
            for r in range(n_rmats):
                tv = all_tstv_vals[r] if r < len(all_tstv_vals) else None
                eq_e = equal_exchangeability if (c == 0 and r == 0) else (tv is None)
                registry[(c, r)] = setup_model(name, tstv=tv,
                    equal_exchangeability=eq_e, freqs=f, equal_freqs=eq_f)
    elif name == "GTR":
        for c in range(n_comps):
            f = all_freqs[c] if c < len(all_freqs) else None
            eq_f = equal_freqs if c == 0 else (f is None)
            for r in range(n_rmats):
                rr = all_rates[r] if r < len(all_rates) else None
                registry[(c, r)] = setup_model(name, freqs=f, equal_freqs=eq_f,
                    relative_rates=rr)
    elif name in AA_MODELS:
        for c in range(n_comps):
            f = all_freqs[c] if c < len(all_freqs) else None
            eq_f = equal_freqs if c == 0 else (f is None)
            for r in range(n_rmats):
                rr = all_rates[r] if r < len(all_rates) else None
                registry[(c, r)] = setup_model(name, freqs=f, equal_freqs=eq_f,
                    relative_rates=rr)
    else:
        raise ValueError(f"Unknown model: {name}")

    return registry
