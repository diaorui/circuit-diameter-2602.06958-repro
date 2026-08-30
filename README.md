# Algorithm 1 stall (arXiv:2602.06958)

Reproduction of a non-termination bug in Algorithm 1 of

> B. Natura, *Circuit diameter of polyhedra is strongly polynomial*, arXiv:2602.06958

The **proofs are fine**. The paper measures Phase 2 progress by two events: the trapped set $T$ grows, or a non-basic coordinate is zeroed. Algorithm 1 only resets the reference index $r$ on the first of these. After a coordinate in $N$ is zeroed, $x^{(r)}$ can still hold a stale nonzero value for it, so a later elimination step produces a direction that tries to decrease an already-zero coordinate. The step size is then $\alpha = 0$, $T$ does not change, $r$ is not reset, and the algorithm stalls.

The fix is to reset $r$ on either progress event.

## Run

```
pip install -r requirements.txt
python reproduce.py
```

Needs Python 3 and sympy. Exact rationals; no floating point.

## What you should see

`reproduce.py` has three steps.

1. **Published Algorithm 1** (reset $r$ only when $T$ changes): the walk hits $\alpha = 0$ and repeats forever.
2. **Not a bad circuit choice:** at that iterate, every circuit that could be used for the elimination step also has $\alpha = 0$.
3. **The fix** (also reset $r$ when a nonbasic coordinate is zeroed): the same instance reaches $`x^*`$, with at most $2m$ resets.

`algorithm1.py` is Algorithm 1. The proposed fix is the `elif reset_on_support_shrink ...` branch.

## Scope

This is a direct, line-by-line implementation of Algorithm 1 as published
(elementary vector enumeration, conformal decomposition, both the norm-
reduction and elimination steps), parameterized by `reset_on_support_shrink`
to switch between the published reset rule and the proposed fix. It works
on arbitrary instances (A, b, x0, $`x^*`$, B, N), not just the stalling example
in `reproduce.py`.

Note: `elementary_vectors()` enumerates circuits by brute-force combinatorial
search over column subsets — this is exponential in n and is meant for
small/illustrative instances, not a strongly-polynomial implementation of
the algorithm itself.

Developed with AI assistance; the reproduction itself is exact-arithmetic and independently verifiable by running the code.
