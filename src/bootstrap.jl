# Shared resampling infrastructure.
#
# Four packages (ERGM.jl, ERGMCount.jl, ERGMRank.jl, ERGMMulti.jl) report
# inverse-pseudo-Hessian standard errors, which are anticonservative whenever
# the model is dyad-dependent, and a fifth (REM.jl) reports inverse-Hessian
# standard errors conditional on ONE draw of the case-control risk set. The
# honest alternative in every case has the same shape:
#
#   1. draw `n_boot` replicate datasets from the fitted model (a parametric
#      bootstrap) or from the sampling design (a repeated draw of the controls),
#   2. refit the SAME estimator on each replicate,
#   3. report the empirical covariance of the refitted coefficients.
#
# That loop is written ONCE, here, and the packages supply the two callbacks
# that differ. It lives in NetworkCore.jl rather than ERGM.jl because REM.jl
# deliberately does not depend on ERGM.jl — NetworkCore.jl is the one foundation
# all five share. The sibling shared optimizer `newton_fit` (src/newton.jl)
# moved here from ERGM.jl for the same reason.
#
# This file also hosts `check_se`, the ONE validator of the `se=` keyword that
# every fitter exposing `:hessian`/`:sandwich`/`:bootstrap` used to re-type
#.
#
# Only stdlib is used (Random, Statistics, LinearAlgebra): NetworkCore.jl stays
# dependency-light.

"""
    bootstrap_cov(refit, simulate, θ̂; n_boot=100, rng=Random.default_rng(),
                  threaded=true) -> NamedTuple

Resampling covariance of an estimator: draw `n_boot` replicate datasets, refit
on each, and return the empirical covariance of the refitted coefficients.

This is the ONE shared bootstrap loop of the ecosystem. It is used for the
parametric bootstrap of a pseudo-likelihood estimator (simulate from the model
at θ̂, refit, empirical covariance — `ERGM.mple`, `ERGMCount.count_mple`,
`ERGMRank.fit_ergm_rank`, `ERGMMulti.ergm_multi`, `TERGM.stergm` with
`se=:bootstrap`) and, with a different `simulate`, for repeated draws of a
*sampling design* rather than of the model: `REM.control_draw_cov` redraws the
case-control risk set and reports the between-draw covariance of the refits as
a control-draw sensitivity *diagnostic*, not as a standard error (REM's former
`se=:bootstrap`, which combined it with the within-draw covariance, was
withdrawn as double-counting — the sampled likelihood's own information is
already consistent under nested case-control sampling).

A refit that has no finite answer on its replicate (a boundary statistic in a
simulated network, a separated design) should return a `NaN` row: the loop
keeps every replicate as returned, and the caller decides what a non-finite
row means for its estimator — every ERGM-family `se=:bootstrap` excludes such
rows from the covariance, warns once with the count in its own words, and
records the exclusion in `approximations(fit)`.

Its sibling for the *point* estimate is [`newton_fit`](@ref), the shared
Newton–Raphson optimizer also hosted in NetworkCore.jl (a fitter's `se=:hessian`
is the `se` that `newton_fit` returns; `se=:bootstrap` replaces it with the
covariance computed here). Validate the `se=` keyword itself with
[`check_se`](@ref).

Replicates are opaque to this function: `simulate` may return anything `refit`
understands (a simulated network, an observations DataFrame, an `(index, data)`
pair — REM uses the last so each threaded refit can stash its within-draw
covariance).

**The point estimate is not touched.** Only the covariance is replaced: callers
keep θ̂ and swap in the returned `vcov`/`se`. The bootstrap centres on the
replicate mean, not on θ̂, so it estimates the sampling variability of the
estimator, not its bias.

# Arguments
- `refit`: `replicate -> AbstractVector{<:Real}` — refits the estimator on one
  replicate dataset and returns its coefficient vector (length `length(θ̂)`).
  Must be pure and thread-safe: it is called from a threaded loop.
- `simulate`: `(rng, n_boot) -> AbstractVector` — draws all `n_boot` replicate
  datasets in ONE call, so that a package with its own chain-splitting or
  seeding discipline (e.g. `ERGM.sample_networks`) keeps it. All randomness must
  flow through the `rng` it is handed, so a fixed `rng` reproduces the standard
  errors exactly.
- `θ̂::AbstractVector`: the point estimate — used for its length and returned
  unchanged by the caller; the covariance does not depend on it.

# Keyword Arguments
- `n_boot::Int=100`: number of replicates (at least 2; a covariance needs two)
- `rng::AbstractRNG=Random.default_rng()`: source of the replicate randomness
- `threaded::Bool=true`: refit the replicates on all available threads. The
  refits are deterministic given their replicate, so this is
  thread-count-independent (all the randomness lives in `simulate`).

# Returns
NamedTuple `(vcov, se, replicates)`:
- `vcov::Matrix{Float64}` — `n × n` empirical covariance of the refits
- `se::Vector{Float64}` — its square-rooted diagonal (negative-zero-guarded)
- `replicates::Matrix{Float64}` — the `n_boot × p` matrix of refitted
  coefficients, so callers can report percentile intervals or diagnose a refit
  that failed to converge

An exception thrown by `refit` propagates as itself in both modes: the
threaded loop unwraps the `CompositeException`/`TaskFailedException` the task
machinery adds and rethrows the first underlying exception.

# Example
```julia
# Bootstrap of a sample mean: the "model" is the fitted normal, the "refit" is
# the mean of a replicate.
using Random, Statistics
θ̂ = [3.0]
boot = bootstrap_cov(x -> [mean(x)],
                     (rng, B) -> [randn(rng, 50) .+ 3.0 for _ in 1:B],
                     θ̂; n_boot=200, rng=MersenneTwister(1))
boot.se[1]   # ≈ 1/sqrt(50)
```
"""
function bootstrap_cov(refit, simulate, θ̂::AbstractVector{<:Real};
                       n_boot::Int=100,
                       rng::Random.AbstractRNG=Random.default_rng(),
                       threaded::Bool=true)
    n_boot >= 2 ||
        throw(ArgumentError("n_boot must be at least 2 to form a covariance " *
                            "(got $n_boot)"))
    p = length(θ̂)
    p >= 1 || throw(ArgumentError("θ̂ must have at least one coefficient"))

    replicates = simulate(rng, n_boot)
    length(replicates) == n_boot ||
        throw(ArgumentError("`simulate` returned $(length(replicates)) " *
                            "replicates, expected n_boot = $n_boot"))

    coefs = Matrix{Float64}(undef, n_boot, p)
    if threaded
        try
            Threads.@threads for b in 1:n_boot
                coefs[b, :] .= refit(replicates[b])
            end
        catch e
            # A refit that throws arrives wrapped by the task machinery;
            # rethrow the underlying exception so that a caller catching, say,
            # `ArgumentError` sees it whether or not the loop is threaded.
            inner = unwrap_task_exception(e)
            inner === e ? rethrow() : throw(inner)
        end
    else
        for b in 1:n_boot
            coefs[b, :] .= refit(replicates[b])
        end
    end

    V = Matrix{Float64}(cov(coefs))
    return (vcov=V, se=sqrt.(max.(diag(V), 0.0)), replicates=coefs)
end

"""
    check_se(se::Symbol, allowed; context::AbstractString) -> Symbol

Validate the `se=` keyword of a fitter: return `se` unchanged when it is one
of `allowed` (a tuple or vector of `Symbol`s), otherwise throw an
`ArgumentError` whose message names `context` (the fitter), the allowed
vocabulary and the offending symbol:

    <context>: se must be one of (:hessian, :bootstrap) (got :sandwich)

This is the ONE `se=` validator of the ecosystem:
`ERGM.mple`, `ERGMCount.count_mple`, `ERGMRank.fit_ergm_rank`,
`ERGMMulti.ergm_multi`, `TERGM.stergm` and `REM.fit_rem` (both its
`EventSequence` and `DataFrame` methods take `(:hessian, :sandwich)`; REM
refuses `:bootstrap` first, with a pointer to `control_draw_cov`) each accept a
different subset of `(:hessian, :sandwich, :bootstrap)`, and each used to carry
its own copy of the `se in (...) || throw(...)` line with its own wording. One
helper, one message shape, one place to extend the vocabulary.

# Example
```julia
check_se(:hessian, (:hessian, :bootstrap); context="mple")   # :hessian
try check_se(:sandwich, (:hessian, :bootstrap); context="mple") catch err; err end
# ArgumentError("mple: se must be one of (:hessian, :bootstrap) (got :sandwich)")
```
"""
function check_se(se::Symbol, allowed; context::AbstractString)
    se in allowed && return se
    throw(ArgumentError("$context: se must be one of $(Tuple(allowed)) (got :$se)"))
end
