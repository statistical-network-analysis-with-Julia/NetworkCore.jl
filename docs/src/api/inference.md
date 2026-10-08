# Shared Inference Utilities

Machinery that more than one model package needs, defined once here so that it
cannot drift into several subtly different copies: the one Newton optimizer and
logistic-likelihood kernel, the one bootstrap loop and `se=` validator, the one
statistic protocol, and the golden-fixture harness that pins numerical
agreement with the R reference implementations.

## Shared Newton Optimizer

`newton_fit` is **the** Newton–Raphson-with-step-halving optimizer of the
ecosystem and `logistic_derivatives` **the** logistic log-likelihood kernel
behind every dyad-independent pseudo-likelihood (ERGM's MPLE, ERGMMulti,
TERGM's CMPLE, ERGMRank's swap MPLE). They live here rather than in ERGM.jl
because REM.jl and Revel.jl, which fit their models with the same
optimizer, deliberately do not depend on ERGM.jl; ERGM.jl
re-exports them, so `using ERGM` keeps working unchanged.

Both names are `public` but **not exported** from NetworkCore.jl: refer to them
as `NetworkCore.newton_fit` / `NetworkCore.logistic_derivatives`, or `import
NetworkCore: newton_fit, logistic_derivatives` (ERGM.jl's test suite calls
`newton_fit` unqualified after `using ERGM, NetworkCore`, and a second export of a
different binding would make that name ambiguous there until ERGM.jl aliases
its copy).

Two things distinguish this optimizer from the copies it replaces:

- **A step at which the objective is not finite is never taken.** It is halved
  like a decrease; if every halving is non-finite the iteration stops at the
  last finite iterate with `converged == false`, instead of stepping onto `NaN`.
- **An indefinite Hessian gives `NaN` standard errors and a warning**, never a
  finite number. Before Cholesky, information is scaled to unit diagonal;
  its smallest eigenvalue must exceed `p * eps(Float64)` times its largest.
  This rejects numerical non-identification without depending on individual
  coefficient units. Failed identification also sets `converged=false`.
  The covariance is unscaled back to the original coefficient units.

`logistic_derivatives` takes binary rows (`y`, optionally `weights=` for
collapsed dyad classes) or binomial rows (`n_tot`, `n_one`), through ONE
kernel whose workspaces are allocated once; an evaluation allocates exactly
the gradient and Hessian it returns.

```julia
using NetworkCore

# Edges-only MPLE of a 16-actor undirected network with 20 ties, in closed
# form: 120 dyads, so θ̂ = logit(20/120) and se = 1/sqrt(120 p (1 − p)).
d   = NetworkCore.logistic_derivatives(ones(1, 1), [120.0], [20.0])
fit = NetworkCore.newton_fit(d, [0.0])
fit.θ[1] ≈ log(20 / 100)                          # true
fit.se[1] ≈ sqrt(1 / (120 * (1 / 6) * (5 / 6)))    # true
fit.converged                                     # true

# A Poisson log-mean from k = 7 events: ll(θ) = kθ − e^θ, maximum at log k
k = 7.0
pfit = NetworkCore.newton_fit(θ -> (k * θ[1] - exp(θ[1]), [k - exp(θ[1])], hcat(-exp(θ[1]))),
                           [0.0])
pfit.θ[1] ≈ log(k) && pfit.se[1] ≈ 1 / sqrt(k)   # true
```

```@docs
NetworkCore.newton_fit
NetworkCore.logistic_derivatives
```

`converged == true` does not mean that a finite maximum exists. On separated
data the optimizer settles on an asymptote. Whether that happened is
answered by the [Separation Contract](contracts.md#Separation-Contract).

## Threaded Loops

Threaded samplers and refits run their tasks through `spawn_all`, which
waits for every task and rethrows the first failure as the task's own
exception, never wrapped in a `TaskFailedException`. A caller that catches an
`ArgumentError` therefore catches it on the threaded path as well as on the
serial one. `unwrap_task_exception` is the unwrapping step on its own, for
loops written with `Threads.@threads` (as `bootstrap_cov` is). Both are
`public`, not exported.

```@docs
NetworkCore.spawn_all
NetworkCore.unwrap_task_exception
```

## Bootstrap Covariance

`bootstrap_cov` is **the** resampling loop of the ecosystem. It backs the
parametric bootstrap of every pseudo-likelihood estimator (`ERGM.mple`,
`ERGMCount.count_mple`, `ERGMRank.fit_ergm_rank`, `ERGMMulti.ergm_multi` with
`se=:bootstrap`) and, with a different `simulate` callback, the repeated draw of
the case-control risk set in `REM.fit_rem`.

All randomness flows through the caller's `rng` and the refits are deterministic
given their replicate, so a fixed `rng` reproduces the standard errors exactly
and the threaded refits are thread-count-independent.

The `se=` keyword every such fitter exposes — a subset of `(:hessian,
:sandwich, :bootstrap)` — is validated by the one helper `check_se`, so a wrong
symbol produces the same message shape everywhere:

```julia
using NetworkCore

check_se(:hessian, (:hessian, :bootstrap); context="mple")    # :hessian
try
    check_se(:sandwich, (:hessian, :bootstrap); context="mple")
catch err
    println(err.msg)   # mple: se must be one of (:hessian, :bootstrap) (got :sandwich)
end
```

```@docs
bootstrap_cov
check_se
```

## The StatsAPI Surface

Every fitted-model type in the ecosystem is expected to implement the StatsAPI
verbs wherever the quantity exists — `coef`, `stderror`, `vcov`, `confint`,
`loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`, `coefnames` — on
its **own** type, with `coeftable` returning a [`CoefficientTable`](@ref) whose
names are `coefnames(fit)`. NetworkCore.jl defines no method on any of them; it
re-exports `coeftable` and `coefnames` (one binding each for the whole
ecosystem) and provides `check_statsapi`, which turns the claim into a testset
line. It reports, verb by verb, whether `hasmethod` finds an implementation,
checks that the ones found agree with each other (`coef`, `stderror`, `vcov`,
`coeftable` and `coefnames` have matching sizes, and `coefnames` equals the
table's names; a method that throws on its own fit counts as absent), and with
`strict=true` throws an error that names exactly what is missing. `coefnames`
is checked whenever a fit defines it, but is not yet required by default; list
it in `required=` to require it:

```julia
using NetworkCore
using LinearAlgebra: diag
import StatsAPI

struct TinyFit
    θ::Vector{Float64}
    V::Matrix{Float64}
end
StatsAPI.coef(f::TinyFit) = f.θ
StatsAPI.vcov(f::TinyFit) = f.V
StatsAPI.stderror(f::TinyFit) = sqrt.(diag(f.V))
StatsAPI.coeftable(f::TinyFit) =
    CoefficientTable(StatsAPI.coefnames(f), StatsAPI.coef(f), StatsAPI.stderror(f))
StatsAPI.coefnames(::TinyFit) = ["a", "b"]

fit = TinyFit([1.0, -0.5], [0.04 0.0; 0.0 0.09])
check_statsapi(fit; required=(:coef, :stderror, :vcov, :coeftable, :coefnames))
# (coef = true, stderror = true, vcov = true, coeftable = true, coefnames = true)

r = check_statsapi(fit)        # the ten verbs required by default
r.coef, r.aic                  # (true, false): no StatsAPI.aic(::TinyFit)
try
    check_statsapi(fit; strict=true)
catch err
    println(err.msg)   # "... missing or inconsistent: confint, loglikelihood, nobs, dof, aic, bic. ..."
end
```

A package's testset pins its surface with
`@test check_statsapi(fit; strict=true) !== nothing` (or a `required=` subset
naming the verbs the model legitimately has — a model without a likelihood has
no `aic`).

```@docs
check_statsapi
```

## Statistic Protocol

Three empty generics — [`compute`](@ref), [`name`](@ref) and
[`compute_all`](@ref). NetworkCore.jl defines no methods and knows nothing about
what a statistic *is*; every model package extends these with methods for its
own statistic types (`compute(term, net)` in ERGM.jl, `compute(term, net,
prev_net)` in TERGM.jl, `compute(stat, state, sender, receiver)` in REM.jl).
The signatures differ, so the methods coexist on one function.

This exists because ERGM.jl and REM.jl each used to define and export their
*own* `compute` and `name` — two different functions with one name — so
`using ERGM, REM` left the unqualified verbs undefined, which is Julia's rule
for conflicting exports. One generic cannot collide with itself.

The rule: **a verb that more than one package needs is defined once, here, and
imported by name (`import NetworkCore: compute, name`) wherever methods are
added** — never redefined locally.

```@docs
compute
name
compute_all
```

## Golden Fixtures

The harness behind the ecosystem's numerical-agreement tests: a fixture is a
stored set of reference values (typically produced by the corresponding R
statnet/RSiena routine) plus a per-key tolerance, and a test checks a computed
value against it.

A fixture only loads when it can be **audited and regenerated**: the
`[provenance]` block must name `r_version`, `seed` and `script`, and the
script must exist on disk. `load_golden` resolves `provenance.script` by
walking up the ancestors of the fixture's directory (or of `root=`, when
given) until `joinpath(dir, script)` is a file — so the conventional
`script = "test/fixtures/r/<name>.R"`, relative to the package root, resolves
from any working directory, and an ad-hoc fixture may keep `r/<name>.R`
beside itself. The resolved absolute path is the fixture's `script_path`
field, printed by `show` and by `golden_report`, so a maintainer facing a red
test sees not only which R and seed produced the reference but which script
to re-run. A fixture whose script cannot be found throws an `ArgumentError`
naming the fixture, the script string and every directory searched.

```julia
using NetworkCore

g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "florentine_sna.toml"))
g.provenance["sna_version"]              # "2.8"
endswith(g.script_path, "harness_selftest.R")   # true, an absolute path
check_golden(g, "density", 1 / 6)        # true
```

```@docs
GoldenFixture
load_golden
check_golden
golden_report
golden_tolerance
```
