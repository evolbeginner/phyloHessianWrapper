"""Nucleotide substitution models: JC69, HKY, F84, GTR."""

import numpy as np
from ._validators import (validate_nuc_freqs, validate_state, validate_length,
                           validate_output_array)

NUM_NUC = 4

A, C, G, T = 0, 1, 2, 3

_GTR_RATE_LABELS = ["AC", "AG", "AT", "CG", "CT", "GT"]


class JC69Model:
    """JC69 (Jukes-Cantor 1969) model: equal base frequencies and equal substitution rates.

    P_ii(t) = 1/4 + 3/4 * exp(-4t/3)
    P_ij(t) = 1/4 - 1/4 * exp(-4t/3)    i != j
    """

    def __init__(self):
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

    def set_matrix(self, matrix, length):
        """Fill a flat array with cumulative transition probabilities.
        Each row's last element is exactly 1.0."""
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(self.num_states, self.num_states)
        a = np.exp(-4.0 * length / 3.0)
        p_same = 0.25 + 0.75 * a
        p_diff = 0.25 - 0.25 * a

        m[:] = p_diff
        for i in range(self.num_states):
            m[i, i] = p_same
        np.cumsum(m, axis=1, out=m)
        m[:, -1] = 1.0

    def set_vector(self, vector, state, length):
        """Fill a flat array with a single cumulative row of P(t).
        The last element is exactly 1.0."""
        validate_output_array(vector, self.num_states)
        length = validate_length(length)
        state = validate_state(state, self.num_states)
        a = np.exp(-4.0 * length / 3.0)
        p_same = 0.25 + 0.75 * a
        p_diff = 0.25 - 0.25 * a

        v = vector[:self.num_states]
        v[:] = p_diff
        v[state] = p_same
        np.cumsum(v, out=v)
        v[-1] = 1.0


class HKYModel:
    """HKY model: closed-form cumulative transition probabilities."""

    def __init__(self, freqs, tstv, equal_exchangeability):
        validate_nuc_freqs(freqs)
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        pi_A, pi_C, pi_G, pi_T = pi
        pi_R = pi_A + pi_G
        pi_Y = pi_C + pi_T
        pi_AG = pi_A * pi_G
        pi_CT = pi_C * pi_T

        if not equal_exchangeability:
            if tstv is None or not np.isfinite(tstv) or tstv <= 0.0:
                raise ValueError(f"tstv must be finite and positive, got {tstv}")
            kappa = (tstv * pi_R * pi_Y) / (pi_AG + pi_CT)
        else:
            kappa = 1.0

        self._pi_A = pi_A; self._pi_C = pi_C; self._pi_G = pi_G; self._pi_T = pi_T
        self._tab1A = pi_A * pi_Y / pi_R
        self._tab2A = pi_G / pi_R
        self._tab3A = pi_A / pi_R
        self._tab1C = pi_C * pi_R / pi_Y
        self._tab2C = pi_T / pi_Y
        self._tab3C = pi_C / pi_Y
        self._tab1G = pi_G * pi_Y / pi_R
        self._tab2G = pi_A / pi_R
        self._tab3G = pi_G / pi_R
        self._tab1T = pi_T * pi_R / pi_Y
        self._tab2T = pi_C / pi_Y
        self._tab3T = pi_T / pi_Y

        self._beta = -1.0 / (2.0 * (pi_R * pi_Y + kappa * (pi_AG + pi_CT)))
        self._beta_A_R = self._beta * (1.0 + pi_R * (kappa - 1.0))
        self._beta_A_Y = self._beta * (1.0 + pi_Y * (kappa - 1.0))

    def set_matrix(self, matrix, length):
        """Fill a flat array (16 elements) with cumulative transition probabilities,
        row-major order.  Each row's last element is exactly 1.0.
        Accepts (16,) or (4,4) input."""
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(-1)
        aa = np.exp(self._beta * length)
        bbR = np.exp(self._beta_A_R * length)
        bbY = np.exp(self._beta_A_Y * length)
        _fill_matrix(m, self._pi_A, self._pi_C, self._pi_G, self._pi_T,
                     self._tab1A, self._tab2A, self._tab3A,
                     self._tab1C, self._tab2C, self._tab3C,
                     self._tab1G, self._tab2G, self._tab3G,
                     self._tab1T, self._tab2T, self._tab3T,
                     aa, bbR, bbY)

    def set_vector(self, vector, state, length):
        """Fill a flat array with a single cumulative row of P(t) for the given starting
        state.  The last element is exactly 1.0."""
        validate_output_array(vector, self.num_states)
        length = validate_length(length)
        state = validate_state(state, self.num_states)
        aa = np.exp(self._beta * length)
        bbR = np.exp(self._beta_A_R * length)
        bbY = np.exp(self._beta_A_Y * length)
        _fill_vector(vector, state, self._pi_A, self._pi_C, self._pi_G, self._pi_T,
                     self._tab1A, self._tab2A, self._tab3A,
                     self._tab1C, self._tab2C, self._tab3C,
                     self._tab1G, self._tab2G, self._tab3G,
                     self._tab1T, self._tab2T, self._tab3T,
                     aa, bbR, bbY)


class F84Model:
    """F84 model: closed-form cumulative transition probabilities."""

    def __init__(self, freqs, tstv, equal_exchangeability):
        validate_nuc_freqs(freqs)
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        pi_A, pi_C, pi_G, pi_T = pi
        pi_R = pi_A + pi_G
        pi_Y = pi_C + pi_T
        pi_AG = pi_A * pi_G
        pi_CT = pi_C * pi_T
        F84temp1 = pi_A**2 + pi_C**2 + pi_G**2 + pi_T**2
        F84temp2 = (pi_A**2 / pi_R) + (pi_C**2 / pi_Y) + (pi_G**2 / pi_R) + (pi_T**2 / pi_Y)

        if not equal_exchangeability:
            if tstv is None or not np.isfinite(tstv) or tstv <= 0.0:
                raise ValueError(f"tstv must be finite and positive, got {tstv}")
            xi = pi_R * pi_Y * (pi_R * pi_Y * tstv - pi_AG - pi_CT)
            xv = pi_CT * pi_R + pi_AG * pi_Y
            kappa = xi / xv
        else:
            kappa = 0.0

        min_pi_ry = min(pi_R, pi_Y)
        if kappa < -min_pi_ry:
            raise ValueError(
                f"tstv={tstv} produces kappa={kappa:.6f} < "
                f"-min(pi_R, pi_Y)={-min_pi_ry:.6f}, "
                f"resulting in negative instantaneous substitution rates. "
                f"Use a larger tstv ratio."
            )

        denom = (1.0 - F84temp1) + (kappa * (1.0 - F84temp2))
        if not np.isfinite(denom) or denom <= 0.0:
            raise ValueError(
                f"invalid F84 rate normalization denominator: {denom}. "
                f"Check tstv and base frequencies."
            )
        mu = -1.0 / denom
        self._mu = mu
        self._mu_kappa_1 = mu * (kappa + 1.0)

        self._pi_A = pi_A; self._pi_C = pi_C; self._pi_G = pi_G; self._pi_T = pi_T
        self._tab1A = pi_A * pi_Y / pi_R
        self._tab2A = pi_G / pi_R
        self._tab3A = pi_A / pi_R
        self._tab1C = pi_C * pi_R / pi_Y
        self._tab2C = pi_T / pi_Y
        self._tab3C = pi_C / pi_Y
        self._tab1G = pi_G * pi_Y / pi_R
        self._tab2G = pi_A / pi_R
        self._tab3G = pi_G / pi_R
        self._tab1T = pi_T * pi_R / pi_Y
        self._tab2T = pi_C / pi_Y
        self._tab3T = pi_T / pi_Y

    def set_matrix(self, matrix, length):
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(-1)
        aa = np.exp(self._mu * length)
        bb = np.exp(self._mu_kappa_1 * length)
        _fill_matrix(m, self._pi_A, self._pi_C, self._pi_G, self._pi_T,
                     self._tab1A, self._tab2A, self._tab3A,
                     self._tab1C, self._tab2C, self._tab3C,
                     self._tab1G, self._tab2G, self._tab3G,
                     self._tab1T, self._tab2T, self._tab3T,
                     aa, bb, bb)

    def set_vector(self, vector, state, length):
        validate_output_array(vector, self.num_states)
        length = validate_length(length)
        state = validate_state(state, self.num_states)
        aa = np.exp(self._mu * length)
        bb = np.exp(self._mu_kappa_1 * length)
        _fill_vector(vector, state, self._pi_A, self._pi_C, self._pi_G, self._pi_T,
                     self._tab1A, self._tab2A, self._tab3A,
                     self._tab1C, self._tab2C, self._tab3C,
                     self._tab1G, self._tab2G, self._tab3G,
                     self._tab1T, self._tab2T, self._tab3T,
                     aa, bb, bb)


class GTRModel:
    """General Time Reversible model: symmetric eigendecomposition of 4x4 rate matrix.

    Parameters
    ----------
    freqs : list of float, length 4
        Equilibrium base frequencies in A, C, G, T order.
    relative_rates : list of float, length 6
        Exchangeabilities in order: r_AC, r_AG, r_AT, r_CG, r_CT, r_GT.
    """

    def __init__(self, freqs, relative_rates):
        validate_nuc_freqs(freqs)
        if len(relative_rates) != 6:
            raise ValueError(
                f"relative_rates must contain exactly six values in order "
                f"AC, AG, AT, CG, CT, GT, got {len(relative_rates)}"
            )
        for i, r in enumerate(relative_rates):
            if not np.isfinite(r):
                raise ValueError(
                    f"relative_rates must be finite, got rate[{i}]={r}"
                )
            if r < 0.0:
                raise ValueError(
                    f"relative_rates must be non-negative, got rate[{i}]={r}"
                )
        if all(r == 0.0 for r in relative_rates):
            raise ValueError("at least one relative rate must be positive")

        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        R = np.zeros((NUM_NUC, NUM_NUC))
        k = 0
        for i in range(NUM_NUC - 1):
            for j in range(i + 1, NUM_NUC):
                R[i][j] = R[j][i] = relative_rates[k]
                k += 1

        Q = R * pi[None, :]
        for i in range(NUM_NUC):
            Q[i][i] = -np.sum(Q[i])

        mean_rate = np.sum(pi * (-np.diag(Q)))
        if not np.isfinite(mean_rate) or mean_rate <= 0.0:
            raise ValueError(f"invalid GTR mean substitution rate: {mean_rate}")
        Q /= mean_rate

        sqrt_pi = np.sqrt(pi)
        inv_sqrt_pi = 1.0 / sqrt_pi

        S = sqrt_pi[:, None] * Q * inv_sqrt_pi[None, :]
        asym = np.max(np.abs(S - S.T))
        if asym > 1e-12:
            raise FloatingPointError(
                f"GTR reversible transform is unexpectedly asymmetric: {asym}"
            )
        S = 0.5 * (S + S.T)

        roots, U = np.linalg.eigh(S)

        self._roots = roots
        self._eigvecs = U
        self._sqrt_pi = sqrt_pi
        self._inv_sqrt_pi = inv_sqrt_pi

    def set_matrix(self, matrix, length):
        """Fill a flat array with cumulative transition probabilities.
        Each row last element is exactly 1.0."""
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(self.num_states, self.num_states)
        if length == 0.0:
            m[:] = np.eye(self.num_states)
        else:
            E = self._eigvecs * np.exp(self._roots * length)[None, :]
            E = E @ self._eigvecs.T
            m[:] = E * (self._sqrt_pi[None, :] / self._sqrt_pi[:, None])
        min_prob = np.min(m)
        if min_prob < -1e-12:
            raise FloatingPointError(
                f"GTR transition matrix contains significantly negative "
                f"probability: {min_prob}"
            )
        m[m < 0.0] = 0.0
        m /= m.sum(axis=1, keepdims=True)
        np.cumsum(m, axis=1, out=m)
        m[:, -1] = 1.0

    def set_vector(self, vector, state, length):
        """Fill a flat array with a single cumulative row of P(t).
        The last element is exactly 1.0."""
        validate_output_array(vector, self.num_states)
        length = validate_length(length)
        state = validate_state(state, self.num_states)
        v = vector[:self.num_states]
        if length == 0.0:
            v[:] = 0.0
            v[state] = 1.0
        else:
            E_row = self._eigvecs[state] * np.exp(self._roots * length)
            v[:] = E_row @ self._eigvecs.T
            v[:] *= self._sqrt_pi / self._sqrt_pi[state]
        min_prob = np.min(v)
        if min_prob < -1e-12:
            raise FloatingPointError(
                f"GTR transition vector contains significantly negative "
                f"probability: {min_prob}"
            )
        v[v < 0.0] = 0.0
        s = v.sum()
        if s > 0.0:
            v /= s
        np.cumsum(v, out=v)
        v[-1] = 1.0


def _fill_matrix(mat, piA, piC, piG, piT,
                 t1A, t2A, t3A, t1C, t2C, t3C,
                 t1G, t2G, t3G, t1T, t2T, t3T,
                 aa, bbR, bbY):
    """Fill a flat 16-element array with cumulative P(t) rows.

    Builds the raw 4x4 transition-probability matrix from HKY/F84 closed-form
    expressions, then applies row-wise cumulative sum.  bbR==bbY for F84.
    """
    m = mat.reshape(4, 4)
    compl = 1.0 - aa  # transversion complement factor

    m[0, 0] = piA + t1A * aa + t2A * bbR
    m[0, 1] = piC * compl
    m[0, 2] = piG + t1G * aa - t3G * bbR
    m[0, 3] = piT * compl

    m[1, 0] = piA * compl
    m[1, 1] = piC + t1C * aa + t2C * bbY
    m[1, 2] = piG * compl
    m[1, 3] = piT + t1T * aa - t3T * bbY

    m[2, 0] = piA + t1A * aa - t3A * bbR
    m[2, 1] = piC * compl
    m[2, 2] = piG + t1G * aa + t2G * bbR
    m[2, 3] = piT * compl

    m[3, 0] = piA * compl
    m[3, 1] = piC + t1C * aa - t3C * bbY
    m[3, 2] = piG * compl
    m[3, 3] = piT + t1T * aa + t2T * bbY

    np.cumsum(m, axis=1, out=m)
    m[:, 3] = 1.0


def _fill_vector(vec, state, piA, piC, piG, piT,
                 t1A, t2A, t3A, t1C, t2C, t3C,
                 t1G, t2G, t3G, t1T, t2T, t3T,
                 aa, bbR, bbY):
    """Fill a 4-element array with a single cumulative row of P(t)."""
    v = vec[:4]
    compl = 1.0 - aa

    if state == A:
        v[0] = piA + t1A * aa + t2A * bbR
        v[1] = piC * compl
        v[2] = piG + t1G * aa - t3G * bbR
        v[3] = piT * compl
    elif state == C:
        v[0] = piA * compl
        v[1] = piC + t1C * aa + t2C * bbY
        v[2] = piG * compl
        v[3] = piT + t1T * aa - t3T * bbY
    elif state == G:
        v[0] = piA + t1A * aa - t3A * bbR
        v[1] = piC * compl
        v[2] = piG + t1G * aa + t2G * bbR
        v[3] = piT * compl
    elif state == T:
        v[0] = piA * compl
        v[1] = piC + t1C * aa - t3C * bbY
        v[2] = piG * compl
        v[3] = piT + t1T * aa + t2T * bbY
    else:
        raise ValueError(f"invalid nucleotide state: {state}")

    np.cumsum(v, out=v)
    v[3] = 1.0
