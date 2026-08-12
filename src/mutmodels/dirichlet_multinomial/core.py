"""
dirichlet multinomial models to compare mutation spectra
"""
import numpy as np
from scipy.special import gammaln
import warnings
from scipy.optimize import minimize_scalar


# support functions
def _to_props(counts):
    return counts / counts.sum(axis=1, keepdims=True)

def _get_stat(P_A, P_B, studentize=True, transform=None, return_all=False):
    Praw_A, Praw_B = P_A, P_B                      # keep untransformed copies
    if transform is not None:
        P_A, P_B = transform(P_A), transform(P_B)

    diff = P_A.mean(0) - P_B.mean(0)
    se = np.sqrt(P_A.var(0, ddof=1)/len(P_A) + P_B.var(0, ddof=1)/len(P_B))
    t = diff / (se + np.median(se)) if studentize else diff

    # calculate total variation difference
    tvd = 0.5* np.sum(np.abs(Praw_A.mean(0) - Praw_B.mean(0)))

    if not return_all:
        return t

    k = np.argmax(np.abs(t))
    return {
        't': t,'winner': k, 'stat': abs(t[k]),
        'pA': Praw_A.mean(0), 'pB': Praw_B.mean(0),
        'TVD':tvd
    }

def _onesample_get_stat(P,z0,studentize=True,transform=None):
    """
    P : per-replicate proportions, shape (n_reps, K)
    z0: reference, already on the same scale as transform(P), shape (K,)
    """
    if transform is not None:
        P = transform(P)

    diff = P.mean(axis=0) - z0

    if studentize:
        se = np.sqrt(P.var(axis=0, ddof=1) / P.shape[0])
        s0 = np.median(se)
        diff = diff / (se + s0)

    return diff


# fitting dirichlet-multinomial model

def dm_loglik_fixed_mean(log_alpha0, counts, p0):
    """
    Dirichlet-multinomial log-likelihood
    with alpha = alpha0 * p0

    the new code that gives correct results for the dm_fixed mean case

    log_alpha0 : scalar (log scale)
    counts     : (R, K) array
    p0         : (K,) fixed mean (must sum to 1)
    """

    alpha0 = np.exp(log_alpha0[0])  # ensure positivity
    alpha = alpha0 * p0
    alpha_sum = alpha0  # because p0 sums to 1

    ll = 0.0

    for x in counts:
        N = x.sum()

        ll += (
            gammaln(N + 1)
            - np.sum(gammaln(x + 1))
            + gammaln(alpha_sum)
            - gammaln(alpha_sum + N)
            + np.sum(gammaln(x + alpha) - gammaln(alpha))
        )

    return ll

def fit_dm_fixed_mean(
        counts,
        p0,
        bounds_log_alpha0=(np.log(1e-3), np.log(1e6)),
        on_failure="warn",):         # "raise" | "warn" | "ignore"
    """
    """
    counts = np.asarray(counts, dtype=float)
    p0 = np.asarray(p0, dtype=float)

    # objective: negative log-likelihood as a function of one scalar, log_alpha0
    def neg_ll(log_alpha0, counts, p0):
        return -dm_loglik_fixed_mean(np.array([log_alpha0]), counts, p0)

    result = minimize_scalar(
        neg_ll,
        bounds=bounds_log_alpha0,
        method="bounded",
        args=(counts, p0),
    )

    alpha0_hat = float(np.exp(result.x))

    # post-fit validation
    lo, hi = bounds_log_alpha0
    problems = []
    if not result.success:
        problems.append(f"optimizer reported failure: {result.message}")
    if not np.isfinite(result.fun):
        problems.append("non-finite log-likelihood at solution")
    if np.isclose(result.x, lo) or np.isclose(result.x, hi):
        problems.append(
            f"solution at a bound (alpha0_hat={alpha0_hat:.4g}); data may be "
            f"near-multinomial or extremely overdispersed"
        )

    if problems:
        msg = "DM fixed-mean fit unreliable: " + "; ".join(problems)
        if on_failure == "raise":
            raise DMFitError(msg)
        elif on_failure == "warn":
            warnings.warn(msg, RuntimeWarning)

    return alpha0_hat, result


class DMFitError(RuntimeError):
    """DM fixed-mean dispersion fit did not converge to a usable optimum."""


# two-condition test
def dm_two_condition(counts_A, counts_B, n_bootstraps=1000, sig_level=0.05,
                             dispersion_type='split', studentize=True, stat_type='max',
                             transform=None, rng=None):
    """Two-condition test using a DM null and a max-type statistic."""

    if min(counts_A.shape[0], counts_B.shape[0]) < 2:
        raise ValueError('need >= 2 replicates per group to studentize')

    pooled_counts = np.vstack([counts_A, counts_B])

    # drop categories empty in both conditions — non-informative, and they'd
    # give p0=0 -> alpha=0 -> gammaln poles in the DM fit
    col_totals = pooled_counts.sum(axis=0)
    empty = col_totals == 0
    if np.any(empty):
        empty_idx = np.where(empty)[0].tolist()
        warnings.warn(
            f"dropping {empty.sum()} categories with no counts in either "
            f"condition (indices {empty_idx}); not informative for the test",
            RuntimeWarning,
        )
        keep = ~empty
        counts_A = counts_A[:, keep]
        counts_B = counts_B[:, keep]
        pooled_counts = pooled_counts[:, keep]

    p0_null = np.mean(_to_props(pooled_counts),axis=0)   # shared mean, all branches

    if dispersion_type == 'shared':
        a0_null, _ = fit_dm_fixed_mean(pooled_counts, p0_null)
        a0_A = a0_B = a0_null
    elif dispersion_type == 'split':
        a0_A, _ = fit_dm_fixed_mean(counts_A, p0_null)
        a0_B, _ = fit_dm_fixed_mean(counts_B, p0_null)
    elif dispersion_type == 'highest':
        a0_A, _ = fit_dm_fixed_mean(counts_A, p0_null)
        a0_B, _ = fit_dm_fixed_mean(counts_B, p0_null)
        a0_A = a0_B = min(a0_A, a0_B)        # lowest alpha0 = highest dispersion
    else:
        raise ValueError(f'dispersion type: {dispersion_type} not a viable option')

    alpha_null_A = a0_A * p0_null
    alpha_null_B = a0_B * p0_null

    rep_sizes_A = counts_A.sum(axis=1)
    rep_sizes_B = counts_B.sum(axis=1)

    P_A, P_B = _to_props(counts_A), _to_props(counts_B)
    print('props A:')
    print(P_A.shape)
    print('props B:')    
    print(P_B.shape)
    obs_diff_stats = _get_stat(P_A, P_B, studentize, transform,return_all=True)
    if stat_type == 'max':
        obs_stat = np.max(np.abs(obs_diff_stats['t']))
    elif stat_type == 'sum':
        obs_stat = np.sum(np.abs(obs_diff_stats['t']))
    else:
        raise ValueError(f'stat type {stat_type} not a valid option')

    rng = np.random.default_rng(rng)
    sim_stats = np.empty(n_bootstraps)

    for b in range(n_bootstraps):
        sim_A = np.vstack([rng.multinomial(N, rng.dirichlet(alpha_null_A))
                           for N in rep_sizes_A])
        sim_B = np.vstack([rng.multinomial(N, rng.dirichlet(alpha_null_B))
                           for N in rep_sizes_B])

        sP_A, sP_B = _to_props(sim_A), _to_props(sim_B)

        cur_stat = _get_stat(sP_A, sP_B, studentize, transform)
        if stat_type == 'max':
            sim_stats[b] = np.max(np.abs(cur_stat))
        elif stat_type == 'sum':
            sim_stats[b] = np.sum(np.abs(cur_stat))

    p_value = (np.sum(sim_stats >= obs_stat) + 1) / (n_bootstraps + 1)

    # compute total variation difference

    return {
        'p_value': p_value,
        'obs_max_stat': obs_stat,
        'winner':obs_diff_stats['winner'],
        'obs_stats': obs_diff_stats['t'],
        'pA':obs_diff_stats['pA'],
        'pB':obs_diff_stats['pB'],
        'stat_cutoff': np.quantile(sim_stats, 1 - sig_level),
        'sim_stats': sim_stats,
        'a0_A': a0_A, 'a0_B': a0_B,
        'TVD':obs_diff_stats['TVD'],
    }

# one-sample test against reference set of proportions (e.g. a cosmic signature)
def dm_onesample(obs_counts, p0, n_bootstraps=1000, sig_level=0.05,
                stat_type='max',studentize=True, transform=None, rng=None):
    """Test a set of replicates against a fixed reference composition p0."""

    if obs_counts.shape[0] < 2:
        raise ValueError('need >= 2 replicates to studentize')

    p0 = np.asarray(p0, dtype=float)
    p0 = p0 / p0.sum() # guard: must sum to 1

    a0_hat, _ = fit_dm_fixed_mean(obs_counts, p0)
    alpha_hat = a0_hat * p0
    rep_sizes = obs_counts.sum(axis=1)

    # reference on both scales, transformed exactly like the data
    z0_t   = transform(p0[None, :])[0] if transform is not None else p0
    z0_raw = p0

    P_obs = _to_props(obs_counts)

    obs_stat = _onesample_get_stat(P_obs, z0_t,   studentize, transform)
    obs_raw  = _onesample_get_stat(P_obs, z0_raw, studentize=False, transform=None)

    if stat_type == 'max':
        obs_stat = np.max(np.abs(obs_stat))
        obs_raw = np.max(np.abs(obs_raw))
    elif stat_type == 'sum':
        obs_stat = np.sum(np.abs(obs_stat))
        obs_raw = np.sum(np.abs(obs_raw))
    else:
        raise ValueError(f'stat type {stat_type} not a valid option')

    rng = np.random.default_rng(rng)
    sim_stats = np.empty(n_bootstraps)
    sim_raws  = np.empty(n_bootstraps)

    for b in range(n_bootstraps):
        sim_counts = np.vstack([
            rng.multinomial(N, rng.dirichlet(alpha_hat)) for N in rep_sizes
        ])
        P_sim = _to_props(sim_counts)
        sim_stat = _onesample_get_stat(P_sim, z0_t,   studentize, transform)
        sim_raw = _onesample_get_stat(P_sim, z0_raw, studentize=False, transform=None)

        if stat_type == 'max':
            sim_stat = np.max(np.abs(sim_stat))
            sim_raw = np.max(np.abs(sim_raw))
        elif stat_type == 'sum':
            sim_stat = np.sum(np.abs(sim_stat))
            sim_raw = np.sum(np.abs(sim_raw))

        sim_stats[b] = sim_stat
        sim_raws[b] = sim_raw

    p_value = (np.sum(sim_stats >= obs_stat) + 1) / (n_bootstraps + 1)

    # compute total variation distance
    tvd = 0.5 * np.sum(np.abs(z0_raw - P_obs.mean(0)))

    return {
        'p_value': p_value,
        'obs_stat': obs_stat,
        'stat_cutoff': np.quantile(sim_stats, 1 - sig_level),
        'sim_stats': sim_stats,
        'alpha0_hat': a0_hat,
        'TVD': tvd,
    }

if __name__ == "__main__":
    pass
