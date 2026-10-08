# Result-Metadata Protocol

Every model package in this ecosystem approximates *something*. ERGM's MPLE
maximises a pseudo-likelihood, not the likelihood; a panel TERGM's CMPLE is the
exact CMLE on a dyad-independent formula and an approximation with
anticonservative standard errors on a dyad-dependent one; a SAOM maximises a
moment condition; a relational-event model maximises a case-control partial
likelihood over a sampled risk set.

None of that should have to be read out of a `show` method or a warning string.
The **result-metadata protocol** makes it programmatically inspectable: seven
accessors that any fitted model may specialise, and one constructor,
[`fit_metadata`](@ref), that collects them into a [`ResultMetadata`](@ref).

```julia
using ERGM, NetworkCore

net = load_dataset(:sampson)
fit = ergm(net, [Edges(), Mutual()]; method=:mple)

meta = fit_metadata(fit)
meta.objective      # :pseudolikelihood
meta.is_exact       # false — the formula is dyad-dependent
meta.se_method      # :hessian  (anticonservative under dependence!)
meta.missing_method # :available_case
```

Each accessor has a conservative default, so a package that has not declared is
reported as `:unspecified` rather than silently assumed exact.

## The Container

```@docs
ResultMetadata
fit_metadata
```

## The Accessors

```@docs
estimand
objective
is_exact
se_method
missing_method
tie_method
approximations
```

`missing_method` and `tie_method` report the policy that was actually applied,
and connect this protocol to the [ecosystem contracts](contracts.md).

## The StatsAPI Surface

The metadata protocol says *how* a fit was obtained; the **StatsAPI surface**
says what can be read off it. Every fitted-model type in the ecosystem answers
the [StatsAPI](https://github.com/JuliaStats/StatsAPI.jl) verbs wherever the
quantity exists — `coef`, `stderror`, `vcov`, `confint`, `loglikelihood`,
`nobs`, `dof`, `aic`, `bic`, `coeftable` and `coefnames` — each defined by
the package on its own result type, never here. NetworkCore.jl contributes the
three shared pieces: it re-exports the [`coeftable`](@ref) and
[`coefnames`](@ref) generics so that every package extends the same bindings
(`NetworkCore.coeftable === StatsAPI.coeftable`, and `using NetworkCore, REM`
leaves the verbs defined), it provides
[`CoefficientTable`](@ref) as the value `coeftable(fit)` returns (the same
vectors the fit's `show` prints), and it provides [`check_statsapi`](@ref),
the one checker each package's testset runs (`check_statsapi(fit;
strict=true)`) so that "implements the surface" is a test, not prose. See
[Shared Inference](inference.md#The-StatsAPI-Surface) for the checker.
