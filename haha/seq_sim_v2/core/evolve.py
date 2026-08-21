"""Monte Carlo evolution engine: root-sequence generation and recursive tree traversal."""

import sys
import numpy as np
from ..utils import random


def evolve_sequences(tree, from_site, num_sites, scale, model_registry, rate_info,
                     ancestor=None):
    """Evolve sequences from root to all tips for one partition.

    If ancestor is provided, its characters are used as the root sequence
    instead of generating a random one from equilibrium frequencies.
    """
    if rate_info.get('invariable'):
        scale /= (1.0 - rate_info['invariable_prop'])

    root_model = model_registry[(tree.root.comp_idx, tree.root.rmat_idx)]
    root_seq = tree.root.sequence
    if ancestor is not None:
        _set_sequence(root_seq, ancestor, from_site, num_sites, root_model)
    else:
        _random_sequence(root_seq, num_sites, root_model.add_freq)

    _set_site_rates(rate_info, from_site, num_sites, scale)

    _evolve_node(tree.root, tree.root.branch1, from_site, num_sites,
                 scale, model_registry, rate_info)
    _evolve_node(tree.root, tree.root.branch2, from_site, num_sites,
                 scale, model_registry, rate_info)
    if not tree.rooted:
        _evolve_node(tree.root, tree.root.branch0, from_site, num_sites,
                     scale, model_registry, rate_info)


def _raise_recursion_limit(num_nodes):
    """Ensure recursion limit exceeds tree depth."""
    min_limit = num_nodes + 100
    if sys.getrecursionlimit() < min_limit:
        sys.setrecursionlimit(max(min_limit, 2000))


def _random_sequence(seq, num_sites, add_freq):
    """Generate a random root sequence drawn from equilibrium frequencies."""
    for i in range(num_sites):
        r = random()
        for s, af in enumerate(add_freq):
            if r < af:
                seq[i] = s
                break


def _set_sequence(seq, source, from_site, num_sites, model):
    """Copy an encoded ancestral sequence into node buffer."""
    chars = model.state_chars
    char_map = {c: i for i, c in enumerate(chars)}
    for i in range(num_sites):
        ch = source[from_site + i]
        idx = char_map.get(ch.upper())
        if idx is None:
            raise ValueError(f"Invalid state '{ch}' in ancestral sequence")
        seq[i] = idx


def _evolve_node(anc, des, from_site, num_sites, scale, model_registry, rate_info):
    """Recursive pre-order: copy parent, mutate, then descend."""
    branch_model = model_registry[(des.comp_idx, des.rmat_idx)]
    length = des.length0 * scale * des.param
    des.sequence[:] = anc.sequence[:]
    _mutate_sequence(des.sequence, from_site, num_sites, length,
                     branch_model, rate_info)
    if des.tip_no == -1:
        _evolve_node(des, des.branch1, from_site, num_sites,
                     scale, model_registry, rate_info)
        _evolve_node(des, des.branch2, from_site, num_sites,
                     scale, model_registry, rate_info)


def _mutate_sequence(seq, from_site, num_sites, length, model, rate_info):
    """Mutate each site independently using the model's P(t).

    Computes per-site effective branch lengths from rate heterogeneity settings
    and applies stochastic state transitions in a single site loop.
    """
    htype = rate_info['type']
    n_states = model.num_states
    num_cats = rate_info.get('num_cats', 1)
    inv = None
    if rate_info.get('invariable_sites') is not None:
        inv = rate_info['invariable_sites'][from_site:from_site + num_sites]

    # --- discrete gamma: precompute category matrices ---
    if htype == 'discrete_gamma':
        cat_rates = rate_info['cat_rates']
        matrices = [np.zeros(n_states * n_states) for _ in range(num_cats)]
        for c in range(num_cats):
            model.set_matrix(matrices[c], cat_rates[c] * length)
        cats = rate_info['categories'][from_site:from_site + num_sites]
        for i in range(num_sites):
            if inv is not None and inv[i]:
                continue
            row = int(seq[i]) * n_states
            seq[i] = _pick_state_from_matrix(matrices[cats[i]], row, n_states)
        return

    # --- uniform / codon / continuous gamma: use set_vector per site ---
    if htype == 'none':
        site_lengths = np.full(num_sites, length)
    elif htype == 'codon':
        site_lengths = np.array(rate_info['cat_rates'], dtype=np.float64)
        site_lengths = site_lengths[np.arange(num_sites) % 3] * length
    elif htype == 'continuous_gamma':
        gr = rate_info['gamma_rates'][from_site:from_site + num_sites]
        site_lengths = gr * length
    else:
        raise ValueError(f"unknown rate heterogeneity type: {htype}")

    vec = np.zeros(n_states)
    for i in range(num_sites):
        if inv is not None and inv[i]:
            continue
        model.set_vector(vec, int(seq[i]), float(site_lengths[i]))
        seq[i] = _pick_state(vec)


def _pick_state(cumulative_vec):
    """Select a state from a cumulative probability vector by random draw."""
    r = random()
    for s, v in enumerate(cumulative_vec):
        if r < v:
            return s
    return len(cumulative_vec) - 1


def _pick_state_from_matrix(matrix, row_start, n_states=4):
    """Select a state from a cumulative row within a probability matrix."""
    r = random()
    for j in range(n_states):
        if r < matrix[row_start + j]:
            return j
    return n_states - 1


def _set_site_rates(rate_info, from_site, num_sites, scale):
    """Precompute per-site relative rates for -wr output."""
    site_rates = rate_info.get('site_rates')
    if site_rates is None:
        return
    htype = rate_info['type']
    sr = site_rates[from_site:from_site + num_sites]
    inv = None
    if rate_info.get('invariable_sites') is not None:
        inv = rate_info['invariable_sites'][from_site:from_site + num_sites]

    if htype == 'continuous_gamma':
        gr = rate_info['gamma_rates'][from_site:from_site + num_sites]
        if inv is not None:
            for i in range(num_sites):
                sr[i] = 0.0 if inv[i] else gr[i] * scale
        else:
            sr[:] = gr * scale
    elif htype == 'discrete_gamma':
        cat_rates = rate_info['cat_rates']
        cats = rate_info['categories'][from_site:from_site + num_sites]
        if inv is not None:
            for i in range(num_sites):
                sr[i] = 0.0 if inv[i] else cat_rates[cats[i]] * scale
        else:
            for i in range(num_sites):
                sr[i] = cat_rates[cats[i]] * scale
    elif htype == 'codon':
        cat_rates = rate_info['cat_rates']
        for i in range(num_sites):
            sr[i] = cat_rates[i % 3] * scale
    else:
        if inv is not None:
            for i in range(num_sites):
                sr[i] = 0.0 if inv[i] else scale
        else:
            sr[:] = scale
