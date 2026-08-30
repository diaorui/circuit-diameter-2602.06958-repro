"""Faithful implementation of Algorithm 1 from
B. Natura, "Circuit diameter of polyhedra is strongly polynomial"
(arXiv:2602.06958), Section 3.2.

reset_on_support_shrink=False  → published pseudocode
    (reset reference index r only when the trapped set T changes)
reset_on_support_shrink=True   → proposed one-line fix
    (also reset r when |supp(x_N)| shrinks)

Exact rational arithmetic via sympy.
"""
from itertools import combinations
from sympy import Matrix, Rational, S


def elementary_vectors(A, m, n):
    """All elementary vectors (circuits) of ker(A), up to sign/scale.
    Support size is in [2, m+1] since rank(A) <= m.
    """
    elem = []
    for size in range(2, m + 2):
        for cols in combinations(range(n), size):
            A_S = A[:, list(cols)]
            if A_S.rank() != size - 1:
                continue
            ns = A_S.nullspace()
            if len(ns) != 1:
                continue
            v = ns[0]
            if any(v[i] == 0 for i in range(size)):
                continue
            full = [S.Zero] * n
            for idx, c in enumerate(cols):
                full[c] = v[idx]
            elem.append(Matrix(full))
    return elem


def conformal_decompose(w, circuits, n):
    """Greedy conformal decomposition of w into elementary vectors.
    The paper: any conformal decomposition of size at most m is valid.
    """
    residual = Matrix(w)
    decomposition = []
    for _ in range(n + 5):
        if all(residual[i] == 0 for i in range(n)):
            break
        supp = {i for i in range(n) if residual[i] != 0}
        found = None
        for v in circuits:
            v_supp = [i for i in range(n) if v[i] != 0]
            if not set(v_supp).issubset(supp):
                continue
            if all((v[i] > 0) == (residual[i] > 0) for i in v_supp):
                found = v
            elif all((-v[i] > 0) == (residual[i] > 0) for i in v_supp):
                found = -v
            if found is not None:
                break
        if found is None:
            raise RuntimeError("no conformal elementary vector found")
        v_supp = [i for i in range(n) if found[i] != 0]
        beta = min(Rational(residual[i], found[i]) for i in v_supp)
        g = beta * found
        decomposition.append(g)
        residual = residual - g
    return decomposition, residual


def aug_P(x, g, n):
    """Maximum feasible step x + alpha*g, alpha = min {-x_i/g_i : g_i < 0}."""
    alphas = [Rational(-x[i], g[i]) for i in range(n) if g[i] < 0]
    if not alphas:
        raise RuntimeError("unbounded augmentation direction")
    alpha = min(alphas)
    return x + alpha * g, alpha


class StallError(Exception):
    def __init__(self, msg, state):
        super().__init__(msg)
        self.state = state


def _supp(v, idxs):
    return [j for j in idxs if v[j] != 0]


def _weighted_l1_on_N(g, x_r, N):
    return sum(abs(g[j]) / x_r[j] for j in N if x_r[j] != 0)


def run_algorithm(A, b, x0, xstar, B, N, m, n,
                  reset_on_support_shrink=False,
                  max_phase2_iters=500, verbose=False):
    """Algorithm 1.

    reset_on_support_shrink=False: published rule (reset r only when T changes).
    reset_on_support_shrink=True : also reset r when a nonbasic coordinate is zeroed.
    """
    circuits = elementary_vectors(A, m, n)
    tau = Rational(1, (2 * m) ** 3)
    lam = Rational(1, (2 * m) ** 2)

    x = Matrix(x0)
    xstar = Matrix(xstar)
    traj = [Matrix(x)]
    i = 0

    # Phase 1: reduce |supp(x_N)| to at most m
    guard = 0
    while len(_supp(x, N)) >= m + 1:
        guard += 1
        if guard > n:
            raise StallError("Phase 1 did not terminate", dict(x=x))
        supp_xN = set(_supp(x, N))
        found = None
        for v in circuits:
            v_supp = [j for j in range(n) if v[j] != 0]
            if not set(v_supp).issubset(supp_xN):
                continue
            if any(v[j] < 0 for j in v_supp):
                found = v
                break
            if any((-v)[j] < 0 for j in v_supp):
                found = -v
                break
        if found is None:
            raise StallError("Phase 1: no valid circuit found", dict(x=x))
        x, _alpha = aug_P(x, found, n)
        traj.append(Matrix(x))
        i += 1

    # Phase 2
    r = i
    T_prev = set()
    n_support_prev = len(_supp(x, N))
    log = []
    reset_count = 0

    for _ in range(max_phase2_iters):
        if all(x[j] == xstar[j] for j in range(n)):
            return dict(success=True, x=x, i=i, log=log,
                        reset_count=reset_count)

        T = {j for j in B if x[j] <= m * xstar[j]}
        did_reset = False
        if T != T_prev:
            r = i
            did_reset = True
            reset_count += 1
        elif reset_on_support_shrink and len(_supp(x, N)) < n_support_prev:
            # proposed fix: also reset r when a nonbasic coordinate was zeroed
            r = i
            did_reset = True
            reset_count += 1
        T_prev = T
        n_support_prev = len(_supp(x, N))
        x_r = traj[r]

        ratios = {}
        for j in N:
            if x_r[j] == 0:
                if x[j] != 0:
                    raise StallError(
                        "x_j^{(r)}=0 but x_j^{(i)}!=0 (invariant violated)",
                        dict(x=x, x_r=x_r, i=i, r=r))
                ratios[j] = S.Zero
            else:
                ratios[j] = x[j] / x_r[j]
        max_ratio = max(ratios.values()) if ratios else S.Zero

        if max_ratio > tau:
            w = xstar - x
            decomp, resid = conformal_decompose(w, circuits, n)
            if any(resid[k] != 0 for k in range(n)):
                raise StallError(
                    "norm-reduction: decomposition leftover nonzero",
                    dict(x=x, w=w, resid=resid))
            gstar = max(decomp, key=lambda g: _weighted_l1_on_N(g, x_r, N))
            x_new, alpha = aug_P(x, gstar, n)
            step_type = "norm-reduction"
        else:
            q = max(N, key=lambda j: ratios[j])
            varrho = ratios[q]
            if varrho == 1:
                raise StallError("rho == 1 (degenerate)", dict(x=x, r=r))
            y = x + (varrho / (1 - varrho)) * (x - x_r)
            z = y + lam * (xstar - y)
            w = z - x
            decomp, resid = conformal_decompose(w, circuits, n)
            if any(resid[k] != 0 for k in range(n)):
                raise StallError(
                    "elimination: decomposition leftover nonzero",
                    dict(x=x, w=w, resid=resid))
            gstar = max(decomp, key=lambda g: -g[q])
            x_new, alpha = aug_P(x, gstar, n)
            step_type = "elimination"

        log.append(dict(i=i, step_type=step_type, alpha=alpha,
                        reset=did_reset, T=frozenset(T)))
        if verbose:
            print(f"iter {i}: {step_type} alpha={alpha} reset={did_reset} "
                  f"T={T} x={list(x)}")

        if alpha == 0:
            raise StallError(
                "alpha=0, x unchanged, infinite loop",
                dict(x=x, x_r=x_r, i=i, r=r, T=T, step_type=step_type))

        x = x_new
        traj.append(Matrix(x))
        i += 1

    return dict(success=False, reason="max_phase2_iters exceeded",
                x=x, i=i, log=log, reset_count=reset_count)
