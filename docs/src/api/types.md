# Types API Reference

This page documents the core data types in NetworkCore.jl.

## Abstract Types

### AbstractNetwork

```@docs
AbstractNetwork
```

## Network Types

### Network

The module is `NetworkCore`; the type it exports is `Network`. So
`using NetworkCore; Network(5)` works.

```@docs
Network
```

### BipartiteNetwork

```@docs
BipartiteNetwork
```

## Constructors

### network

```@docs
network
```

### network\_initialize

```@docs
network_initialize
```

## Missing-Dyad Mask

Dyads whose tie status is unobserved (statnet-style NA ties). A masked dyad
is *not* the same as an absent edge: structural queries keep reporting the
stored face value, while statistical packages consult the mask to exclude
the dyad from estimation. See the [Missing Data guide](../guide/missing_data.md).

```@docs
set_missing_dyad!
is_missing_dyad
delete_missing_dyad!
clear_missing_dyads!
missing_dyads
n_missing_dyads
```

## Shared Result Presentation

The presentation layer used by every model package in the ecosystem:
coefficient tables, the z → p-value step, p-value formatting, Monte-Carlo
p-values, and the goodness-of-fit containers behind the single `gof` generic.

[`z_pvalues`](@ref) is the ONE z → two-sided-normal-p helper: it uses
`erfc(|z|/√2)` so the tail stays accurate where the naive `2(1 − Φ(|z|))`
cancels to a false `0.0`, floors the result at `floatmin(Float64)` for finite
`z` (a finite statistic never has a p-value of exactly zero), and keeps a `NaN`
`z` — an undefined standard error — as a `NaN` p-value. The two-argument form
builds the z statistics from estimates and standard errors first (`NaN` where
the standard error is not strictly positive), so a `show` method is one call
away from a coefficient table:

```julia
using NetworkCore

est = [1.20, -0.30, 0.80]
se  = [0.40,  0.10, 0.00]          # the third SE is undefined
zp  = z_pvalues(est, se)           # (z = [3.0, -3.0, NaN], p = [0.0027, 0.0027, NaN])
print_coeftable(stdout, ["a", "b", "c"], est, se, zp.p; z_values=zp.z)

z_pvalues([50.0])[1] == floatmin(Float64)   # true: floored, never 0.0
```

`print_coeftable` also renders a table **without** a standard-error column
(`std_errors=nothing`, `z_values` mandatory, a three-title `header`) for a
statistic that has none — a permutation test — and takes the display floor
`p_floor=`, so a Monte-Carlo p-value with resolution `1/reps` prints as
`"<0.001"` rather than `0.0`.

[`CoefficientTable`](@ref) is the *inspectable* form of the same table: the
value every model package returns from [`coeftable`](@ref)`(fit)` — the
StatsAPI verb, re-exported here so that all packages extend the ONE binding.
It stores exactly the vectors `print_coeftable` consumes, fills `z`/`p` from
`z_pvalues` when they are omitted, and its `show` **is** `print_coeftable`, so
what a user inspects is what the fit prints:

```julia
using NetworkCore

tbl = CoefficientTable(["edges", "mutual"], [-2.1, 1.3], [0.2, 0.4])
tbl.p_values == z_pvalues([-2.1, 1.3], [0.2, 0.4]).p   # true
tbl["mutual"].estimate                               # 1.3
length(tbl)                                          # 2
[r.name for r in tbl]                                # ["edges", "mutual"]
tbl                                                  # prints the R-style table

# A permutation test: no Std.Error column, a t statistic, resolution 1/1000
CoefficientTable(["intercept", "x"], [0.12, 0.85], nothing;
                 z_values=[1.02, 6.7], p_values=[0.31, 0.0],
                 header=("Estimate", "t value", "Pr(>=|t|)"), p_floor=1/1000)
```

```@docs
print_coeftable
CoefficientTable
coeftable
coefnames
z_pvalues
format_pvalue
signif_code
SIGNIF_LEGEND
mc_pvalue
GOFStatistic
GOFResult
n_simulations
gof
```

For what a fit *did* — the objective it maximised, whether that objective is
exact, how the standard errors were computed — see the
[Result-Metadata Protocol](metadata.md).
