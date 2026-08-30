# Algorithm 1 stall (arXiv:2602.06958)

Reproduction of a non-termination bug in Algorithm 1 of

> B. Natura, *Circuit diameter of polyhedra is strongly polynomial*, arXiv:2602.06958

The **proofs are fine**. The gap is only in the pseudocode: the reference
index `r` is reset when the trapped set `T` changes, but not when a
non-basic coordinate is zeroed — the other progress event the
\(O(m^2\log m)\) analysis counts on.

## Run

```
pip install -r requirements.txt
python reproduce.py
```

Needs Python 3 and sympy. Exact rationals; no floating point.

## What you should see

1. **Published rule** (reset `r` only when `T` changes): the walk hits
   `alpha = 0` and repeats forever on a concrete instance (`m=3`, `n=7`).
2. **Decomposition-independent**: every elementary vector that can appear
   in a conformal decomposition of `z - x` and that makes progress on the
   intended coordinate `q` also decreases an already-zero coordinate, so
   `alpha = 0` is forced.
3. **One-line fix** (also reset `r` when `|supp(x_N)|` shrinks): the same
   instance converges, with at most `2m` resets.

`algorithm1.py` is Algorithm 1. The proposed fix is the
`elif reset_on_support_shrink ...` branch.
