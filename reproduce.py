"""
Reproduction: a concrete instance on which Algorithm 1 of
    B. Natura, "Circuit diameter of polyhedra is strongly polynomial"
    (arXiv:2602.06958), Section 3.2
fails to terminate, and a one-line change to the reference-point reset
rule that fixes it.

Exact rational arithmetic (sympy): "stuck" and "zero" are not float artifacts.

THE BUG
--------
Algorithm 1 resets the reference index r only when the trapped set T
changes (`if T^(i) != T^(i-1): r <- i`). It does not reset r when a
non-basic coordinate is zeroed — the other progress event that the
O(m^2 log m) bound counts on.

Once some coordinate l in N is zeroed it stays zero. But x^(r) can still
hold a stale nonzero value for l (r was not advanced). A later elimination
step builds
    y = x + (rho/(1-rho))(x - x^(r))
    z = y + lambda(x* - y)
using that stale x_l^(r), so z_l < 0 even though x_l is already 0.
Then every conformal decomposition of z - x has a component with a
negative weight on l. Since x_l = 0 already, augmenting along any such
component has step size alpha = 0: a null step. If every circuit that
helps the intended target q also decreases some already-zero coordinate,
alpha = 0 is forced no matter which conformal decomposition is chosen,
T does not change, r is not reset, and the algorithm repeats forever.

THE FIX
--------
Reset r on either progress event: T grows, or |supp(x_N)| shrinks.
This matches the proof's own accounting (at most 2m progress events).

This script:
  1. Runs published Algorithm 1 on a concrete instance (m=3, n=7); it loops.
  2. Checks that every circuit usable for that elimination step has alpha=0.
  3. Runs the one-line fix on the same instance; it converges.
"""
from sympy import Matrix, Rational
from algorithm1 import elementary_vectors, run_algorithm, StallError

A = Matrix([
    [-1, -1,  3,  0,  2,  1, -3],
    [ 1, -3,  2, -2, -3,  3,  0],
    [ 0,  0,  3,  1,  2,  3, -3],
])
m, n = 3, 7
b = Matrix([3, -11, 9])
B = [1, 2, 3]
N = [0, 4, 5, 6]
xstar = Matrix([0, 3, 2, 3, 0, 0, 0])
x0 = Matrix([Rational(7, 4), Rational(23, 12), 1, 1, 3, Rational(2, 3), 1])

assert A * xstar == b
assert A * x0 == b
assert all(v >= 0 for v in xstar) and all(v >= 0 for v in x0)

print("=" * 70)
print("STEP 1: published Algorithm 1 (reset r only when T changes)")
print("=" * 70)
try:
    run_algorithm(A, b, x0, xstar, B, N, m, n,
                  reset_on_support_shrink=False,
                  max_phase2_iters=200, verbose=True)
    print("(unexpectedly terminated — bug not reproduced)")
    raise SystemExit(1)
except StallError as e:
    stall = e.state
    print(f"\n>>> STALL: {e}\n")

print("=" * 70)
print("STEP 2: not a bad circuit choice")
print("=" * 70)
if stall.get("step_type") != "elimination":
    print(f"stall was {stall.get('step_type')}, not elimination; Step 2 skipped")
    raise SystemExit(1)

x_i = stall["x"]
x_r = stall["x_r"]
lam = Rational(1, (2 * m) ** 2)

ratios = {j: (x_i[j] / x_r[j] if x_r[j] != 0 else 0) for j in N}
q = max(N, key=lambda j: ratios[j])
varrho = ratios[q]
y = x_i + (varrho / (1 - varrho)) * (x_i - x_r)
z = y + lam * (xstar - y)
w = z - x_i

already_zero = [j for j in N if x_i[j] == 0]
stale = {j: x_r[j] for j in already_zero if x_r[j] != 0}
neg_already_zero = [j for j in already_zero if w[j] < 0]
print(f"already-zero nonbasic coords: {already_zero}")
print(f"still nonzero in stale x^(r): {stale}")
print(f"z - x is negative on already-zero coords {neg_already_zero} "
      f"→ decreasing them forces α=0")


def can_appear_in_conformal_decomp(g, w):
    """g can appear in some conformal decomposition of w iff
    supp(g) ⊆ supp(w) and g, w have matching signs on supp(g)."""
    support = [i for i in range(n) if g[i] != 0]
    return all(w[i] != 0 and (g[i] > 0) == (w[i] > 0) for i in support)


def forces_null_step(g, x):
    return any(g[j] < 0 and x[j] == 0 for j in range(n))


ev = elementary_vectors(A, m, n)
compatible = []
for g in ev:
    if g[q] == 0:
        continue
    for cand in (g, -g):
        if can_appear_in_conformal_decomp(cand, w):
            compatible.append(cand)

n_null = sum(1 for g in compatible if forces_null_step(g, x_i))
print(f"\ncircuits usable for this elimination step: {len(compatible)}")
print(f"of those with α=0: {n_null} / {len(compatible)}")
if compatible and n_null == len(compatible):
    print(">>> every usable circuit has α=0")

print()
print("=" * 70)
print("STEP 3: fixed Algorithm 1 (reset r on either progress event)")
print("=" * 70)
res = run_algorithm(A, b, x0, xstar, B, N, m, n,
                    reset_on_support_shrink=True,
                    max_phase2_iters=200, verbose=True)
print(f"\n>>> success={res['success']}, iterations={res['i']}, "
      f"resets={res['reset_count']} (bound is 2m={2 * m})")
