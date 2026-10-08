# Shared Newton–Raphson optimizer and logistic-likelihood derivatives.
#
# Every likelihood-based fitter in the ecosystem — ERGM's MPLE, the
# pseudo-likelihoods of ERGMMulti/TERGM/ERGMRank/ERGMCount, REM's conditional
# logit and Revel's ordinal and interval models — maximises a smooth concave
# objective by Newton–Raphson with step halving, and most of them over the
# SAME logistic log-likelihood. Three packages carried the same optimizer and
# four the same logistic loop before 0.2.0. Both live here now, in the one foundation every model package already
# depends on (REM.jl and Revel.jl deliberately do not depend on ERGM.jl,
# which is also why `bootstrap_cov` is here).
#
# Only stdlib LinearAlgebra is used (`cholesky`, `Symmetric`, `mul!`, `norm`):
# NetworkCore.jl stays dependency-light.
#
# The names are `public`, not exported (see the `public` block in NetworkCore.jl):
# ERGM.jl re-exports them as `import NetworkCore: newton_fit, logistic_derivatives`
# + `export`, so `using ERGM` keeps working unchanged while there is exactly ONE
# definition.

using LinearAlgebra: eigvals

"""
    newton_fit(loglik_grad_hess, θ0::AbstractVector;
               maxiter=100, tol=1e-8, max_halvings=30,
               information_rtol=length(θ0)*eps(Float64)) -> NamedTuple

Maximize a smooth objective (typically a log-likelihood) by Newton–Raphson
with step halving.

`loglik_grad_hess(θ)` must return a tuple `(ll, grad, hess)`: the objective
value, its gradient vector, and its Hessian matrix at `θ` (the Hessian of the
objective itself, i.e. negative-definite near a maximum). Each Newton step
`−hess \\ grad` is halved up to `max_halvings` times until the objective is
finite and does not decrease (a decrease within a few ulps of `|ll|` counts as
"does not decrease" when the gradient norm shrinks with it: rounding noise at
the optimum is not a descent). Convergence is declared when the **Newton
decrement** `½ ∇ℓᵀ(−H)⁻¹∇ℓ` — the objective gain the full Newton step could
still deliver, scale-free in the parameters and in the same units as `tol` —
is below `tol` (the full step is then taken, as the polish that lands on the
maximum), or when the objective change is below `tol` and the gradient norm
below `sqrt(tol)`.

This is the ONE Newton optimizer of the ecosystem, hosted in NetworkCore.jl so
that packages which do not depend on ERGM.jl (REM.jl, Revel.jl) share it
too; ERGM.jl re-exports it. Pair it with [`logistic_derivatives`](@ref) for
any logistic (pseudo-)likelihood.

# Arguments
- `loglik_grad_hess`: Function `θ -> (ll, grad, hess)`
- `θ0::AbstractVector`: Starting values (copied, not mutated). The objective
  must be finite there (`ArgumentError` otherwise).

# Keywords
- `maxiter::Int=100`: Maximum Newton iterations
- `tol::Float64=1e-8`: Convergence tolerance on the objective change
- `max_halvings::Int=30`: Maximum step halvings per iteration, i.e. the smallest
  step tried is `2^-(max_halvings-1)` of the Newton step. Far from the maximum
  a likelihood with a rare indicator covariate proposes steps hundreds of times
  too long (a participation shift in a relational event model needed `2^-10`);
  the limit was 10 before 0.2.0, which stopped such fits one halving short.
- `information_rtol::Float64=length(θ0)*eps(Float64)`: Relative eigenvalue
  threshold for the final observed-information rank check, in `[0, 1)`.
  Objectives whose Hessians accumulate longer floating-point sums can supply
  a larger threshold reflecting their evaluation precision.

# Returns
NamedTuple `(θ, se, vcov, loglik, converged, iterations, stop)`:
- `vcov` is `inv(-hess)` at the final iterate, computed through a Cholesky
  factorization of diagonally equilibrated `-hess` — the usual observed-information covariance when
  the objective is a log-likelihood — and `se = sqrt.(diag(vcov))`.
- **If `-hess` is not positive definite** at the solution (a saddle, a flat
  direction, a non-identified coefficient), the covariance is *undefined*:
  `se` and `vcov` are filled with `NaN`, `converged` is false, and a warning is emitted.
  Numerical identification requires the smallest eigenvalue of the information
  scaled to unit diagonal to exceed `information_rtol` times its largest
  (by default, `p * eps(Float64)`).
  This criterion is invariant to changes of individual coefficient units;
  it detects singular information that roundoff can let Cholesky accept. The old
  `sqrt.(abs.(diag(pinv(-hess))))` would have returned finite, meaningless
  numbers here.
- `stop` says why the iteration ended, so a caller can name the reason
  instead of guessing it from `iterations`:
  - `:converged`;
  - `:maxiter` — `maxiter` iterations were taken;
  - `:singular_start` — there is no Newton direction at the starting value:
    the information (`-hess`) is singular there, so the step cannot be
    formed or is useless (typically collinear statistics). Raising `maxiter`
    does not help;
  - `:singular` — the same at a later iterate;
  - `:no_ascent` — no step halving found a finite, non-decreasing objective;
  - `:singular_information` — the iteration settled, but the information
    at the solution is singular, so the covariance is undefined.
- `converged == false` when `maxiter` is exhausted, when the Newton step
  cannot be formed (singular Hessian), or when no step halving finds a
  finite, non-decreasing objective. Only a failed factorisation of the
  Hessian is treated as "no Newton direction"; any other error raised while
  forming the step — a gradient of the wrong length (`DimensionMismatch`), a
  Hessian that is not a matrix (`MethodError`), an interrupt — is a bug in
  `loglik_grad_hess` and propagates instead of being reported as
  non-convergence. **A step at which the objective is not
  finite (`NaN`, `±Inf`) is never taken**: it is halved like a decrease, and
  if every halving is non-finite the iteration stops at the last finite
  iterate. When the halvings only find finite *decreasing* objectives, `θ`
  is left where it is and `converged` reports whether the gradient is
  already below `sqrt(tol)`. Before any halving, an iterate whose Newton
  decrement is below `tol` takes the full step — provided it is finite and
  loses no more than `max(tol, 8·eps(|ll|))`, which the quadratic model
  bounds — and is converged there: this is the polish R's `glm` reports as
  converged at the same point, and it is what makes the verdict
  scale-free.

The decrement criterion is what makes the verdict independent of the scale of the design:
on a pseudo-likelihood over thousands of rows with a statistic running to
~30 the gradient's own units put `sqrt(tol)` below the rounding floor of the
objective, so Newton stalled one ulp short of the maximum with
`‖∇ℓ‖ ≈ 1e-4` and reported a converged fit as unconverged — and BLAS
rounding then decided which bootstrap replicates counted as converged.

Callers should surface `converged == false` to the user (an unconverged fit
is a loud result in this ecosystem, never a silent one).

# `converged` does not certify a finite maximiser

On (quasi-)separated logistic data, a zero-count Poisson, or any objective
whose supremum is approached only as `θ → ±∞` along some direction, the
Newton decrement decays along that asymptote and the fit stops at a large,
finite `θ` (typically 20–40 for a logistic model) with enormous standard
errors and **`converged == true`**, exactly as R's `glm` does. `converged`
says the iteration settled, not that a finite maximum exists. Whether one
exists is a property of the data, decided by the separation verdict:
[`logistic_separation`](@ref), [`poisson_separation`](@ref) or
[`clogit_separation`](@ref). Every fitting function calls one of them and
applies the shared policy ([`warn_separation`](@ref)).

# Example
```julia
# Poisson mean via the log-likelihood of k events: ll(θ) = kθ − exp(θ)
import NetworkCore: newton_fit, logistic_derivatives    # public, not exported
k = 7.0
fit = newton_fit(θ -> (k*θ[1] - exp(θ[1]), [k - exp(θ[1])], hcat(-exp(θ[1]))),
                 [0.0])
fit.θ[1] ≈ log(k)         # true
fit.se[1] ≈ 1 / sqrt(k)   # true: observed information is k at the MLE
fit.converged             # true
fit.stop                  # :converged

# Collinear columns: the information is singular at the start, which is not
# an iteration limit
X = hcat(ones(6), [1.0, 2, 3, 1, 2, 3], [2.0, 4, 6, 2, 4, 6])
col = newton_fit(logistic_derivatives(X, Bool[1, 0, 1, 0, 1, 1]), zeros(3))
(col.converged, col.stop)  # (false, :singular_start)

# Complete separation: "converged", yet the estimate sits on an asymptote;
# the separation verdict says so
x = collect(-5.0:5.0)
sep = newton_fit(logistic_derivatives(hcat(ones(11), x), x .> 0), zeros(2))
sep.converged                                                        # true
NetworkCore.logistic_separation(hcat(ones(11), x), x .> 0).separated # true
```
"""
function newton_fit(loglik_grad_hess, θ0::AbstractVector{<:Real};
                    maxiter::Int=100, tol::Float64=1e-8, max_halvings::Int=30,
                    information_rtol::Float64=length(θ0)*eps(Float64))
    0 <= information_rtol < 1 ||
        throw(ArgumentError("newton_fit: information_rtol must be in [0, 1)"))
    θ = Vector{Float64}(θ0)
    p = length(θ)
    ll, grad, hess = loglik_grad_hess(θ)
    isfinite(ll) ||
        throw(ArgumentError("newton_fit: the objective is not finite at the " *
                            "starting values (got $ll); start inside the domain"))
    # A mis-shaped gradient or Hessian is a bug in the closure, and `\` is
    # lenient enough (1×1 and diagonal cases broadcast) that it would not
    # always surface it; say so here rather than iterate on nonsense.
    length(grad) == p ||
        throw(DimensionMismatch("newton_fit: the gradient has $(length(grad)) entries " *
                                "but θ has $p"))
    size(hess) == (p, p) ||
        throw(DimensionMismatch("newton_fit: the Hessian is $(size(hess)) but θ has $p entries"))
    converged = false
    iterations = 0
    stop = :maxiter
    # Is the information numerically singular (rank-deficient, not merely
    # indefinite) at the starting value? Decided once, scale-free, so that a
    # fit stopping at its first iterate can say WHY.
    singular_at_start = _information_singular(hess, p, information_rtol)

    for iter in 1:maxiter
        iterations = iter
        step = try
            -hess \ grad
        catch e
            # Only a failed factorisation means "no Newton direction". Any
            # other error — a wrong-length gradient (DimensionMismatch), a
            # non-matrix Hessian (MethodError), an interrupt — is the caller's
            # bug and must surface, never be reported as non-convergence.
            e isa Union{SingularException, LAPACKException, ZeroPivotException,
                        PosDefException} || rethrow()
            stop = iter == 1 ? :singular_start : :singular
            break            # singular Hessian: no Newton direction
        end
        if !all(isfinite, step)
            stop = iter == 1 ? :singular_start : :singular
            break
        end
        noise = 8 * eps(abs(ll))
        grad_norm = norm(grad)

        # Newton decrement: ½ gᵀ(−H)⁻¹g, the objective gain the full step
        # predicts — scale-free in θ, in the units of `tol` (negative only
        # when −H is not positive definite along g, and then it says
        # nothing). Below `tol` the iterate is converged in the objective's
        # own terms, and the FULL step is the polish that lands on the
        # maximum (R's glm stops at the same point): take it whenever it is
        # finite and does not lose more than the model itself could predict
        # (`tol`) or rounding can produce (`noise`) — a halved step would
        # leave the gradient half-way, one ulp of objective short.
        decrement = 0.5 * dot(grad, step)
        if 0 <= decrement < tol
            ll_full, grad_full, hess_full = loglik_grad_hess(θ .+ step)
            if isfinite(ll_full) && ll_full >= ll - max(tol, noise)
                θ .+= step
                ll, grad, hess = ll_full, grad_full, hess_full
                converged = true
                stop = :converged
                break
            end
            # (else: the quadratic model is wrong here; halve as usual)
        end

        # Step halving: shrink the Newton step until the objective is finite
        # and does not decrease. A non-finite objective (a step outside the
        # domain, an overflow) counts as a decrease. A decrease within a few
        # ulps of |ll| is rounding noise, not descent, PROVIDED the gradient
        # shrinks with it — at the optimum the last Newton step routinely
        # lowers the objective by one ulp while cutting the gradient by
        # orders of magnitude, and refusing it left the fit one step short.
        stepsize = 1.0
        ll_new, grad_new, hess_new = ll, grad, hess
        finite = false
        improved = false
        for _ in 1:max_halvings
            ll_new, grad_new, hess_new = loglik_grad_hess(θ .+ stepsize .* step)
            finite = isfinite(ll_new)
            if finite && (ll_new >= ll ||
                          (ll_new >= ll - noise && norm(grad_new) < grad_norm))
                improved = true
                break
            end
            stepsize /= 2
        end

        if !improved
            # No halving found a finite, non-decreasing objective: do not move.
            # A finite but decreasing objective at every halving means the
            # current point is (numerically) as good as Newton can do; it is
            # the optimum only if the gradient already says so (a decrement
            # below `tol` was handled above).
            converged = finite && grad_norm < sqrt(tol)
            stop = converged ? :converged :
                   (iter == 1 && singular_at_start) ? :singular_start : :no_ascent
            break
        end

        θ .+= stepsize .* step
        ll_change = abs(ll_new - ll)
        ll, grad, hess = ll_new, grad_new, hess_new

        if ll_change < tol && norm(grad) < sqrt(tol)
            converged = true
            stop = :converged
            break
        end
    end

    vcov, se = _observed_information_cov(hess, p; rtol=information_rtol)
    if converged && !all(isfinite, vcov)
        converged = false
        stop = :singular_information
    end
    return (θ=θ, se=se, vcov=vcov, loglik=ll, converged=converged,
            iterations=iterations, stop=stop)
end

# Is the information `-hess` numerically rank-deficient? A zero (or
# non-finite) diagonal is; otherwise, scaled to unit absolute diagonal, its
# smallest absolute eigenvalue is at most `rtol` times the largest. An
# indefinite but full-rank information (a minimum, a saddle) is not singular.
function _information_singular(hess, p::Int, rtol::Float64)
    p == 0 && return false
    H = Matrix{Float64}(-hess)
    all(isfinite, H) || return true
    d = abs.(diag(H))
    all(>(0), d) || return true
    scales = sqrt.(d)
    C = (H ./ scales) ./ scales'
    all(isfinite, C) || return true
    λ = abs.(eigvals(Symmetric(C)))
    return minimum(λ) <= rtol * maximum(λ)
end

# Covariance from the (negative) Hessian at the solution: `inv(-hess)` via a
# Cholesky factorization after a scale-invariant numerical rank check; NaN
# (with a warning) otherwise. Cholesky alone can accept a singular matrix
# when roundoff leaves an arbitrarily small positive final pivot.
function _observed_information_cov(hess, p::Int; rtol::Float64=p*eps(Float64))
    H = Matrix{Float64}(-hess)
    size(H) == (p, p) ||
        throw(ArgumentError("newton_fit: the Hessian is $(size(H)) for $p parameters"))
    p == 0 && return zeros(0, 0), Float64[]
    if all(isfinite, H) && all(>(0), diag(H))
        scales = sqrt.(diag(H))
        C = (H ./ scales) ./ scales'
        spectrum = all(isfinite, C) ? eigvals(Symmetric(C)) : [NaN]
        if first(spectrum) > rtol * last(spectrum)
            F = cholesky(Symmetric(C); check=false)
            if issuccess(F)
                vcov = (Matrix{Float64}(inv(F)) ./ scales) ./ scales'
                if all(isfinite, vcov) && all(>(0), diag(vcov))
                    return vcov, sqrt.(diag(vcov))
                end
            end
        end
    end
    @warn "newton_fit: the Hessian at the solution is not negative definite or numerically identifiable; " *
          "standard errors are undefined (returned as NaN)"
    return fill(NaN, p, p), fill(NaN, p)
end

"""
    logistic_derivatives(X, y::AbstractVector{Bool};
                         offset=nothing, weights=nothing) -> Function
    logistic_derivatives(X, n_tot::AbstractVector{<:Real},
                         n_one::AbstractVector{<:Real}; offset=nothing) -> Function

The `(ll, grad, hess)` closure of a logistic log-likelihood on design matrix
`X` (`n × p`), ready for [`newton_fit`](@ref).

**Binary rows** (`y`, optionally with non-negative `weights`, e.g. the
multiplicity of a collapsed dyad class):

    ℓ(β) = Σ_r w_r [ y_r η_r − log(1 + e^{η_r}) ],   η = Xβ + offset
    ∇ℓ   = X' w (y − p),   ∇²ℓ = −X' diag(w p(1−p)) X.

**Binomial rows** (`n_tot` trials of which `n_one` are successes, so that a
row stands for a whole class of identical dyads — the form ERGM's MPLE uses
after collapsing dyads with identical change statistics):

    ℓ(β) = Σ_r [ n_one,r η_r − n_tot,r log(1 + e^{η_r}) ]
    ∇ℓ   = X'(n_one − n_tot ⊙ p),   ∇²ℓ = −X' diag(n_tot ⊙ p(1−p)) X.

Both forms run through ONE kernel: the binary form is the binomial one with
`n_tot = weights` (or 1) and `n_one = n_tot ⋅ y`, so a weighted binary fit
and the equivalent duplicated-row fit give identical `ll`, `grad` and `hess`
to rounding. `log(1 + e^η)` is evaluated on its numerically stable branch on
both sides of zero.

Every ERGM-family pseudo-likelihood over dyad-independent rows is this
likelihood — `ERGM.mple`'s logistic regression on the change statistics,
`ERGMMulti`'s MPLE over the within-layer dyads, `TERGM`'s CMPLE over the free
dyads of the auxiliary networks, `ERGMRank`'s swap MPLE over the (ego,
alter-pair) comparisons (with `y ≡ true`: the observed order is always the
"success") — so it lives here, next to `newton_fit`, rather than pasted into
four packages. `offset` (a fixed per-row addition to the linear predictor) is
`ergm`'s offset mechanism: `ERGMMulti` uses it for terms whose coefficients
are held fixed.

**Allocation** : the workspaces are allocated ONCE, when
the closure is built. Each evaluation allocates only the length-`p` gradient
and `p×p` Hessian it returns — never the `n×p` weighted design or a per-row
`x * x'` outer product. The "Shared Newton optimizer" testset pins this with
`@allocated`, and the downstream packages mirror it.

Throws `ArgumentError` on a length mismatch, a negative or non-finite weight,
or a binomial row with `n_one` outside `[0, n_tot]`.

# Example
```julia
# Edges-only MPLE of the Florentine marriage network in closed form: 120
# dyads, 20 ties, so θ̂ = logit(20/120) and se = 1/sqrt(120 p (1 − p)).
import NetworkCore: newton_fit, logistic_derivatives    # public, not exported
d = logistic_derivatives(ones(1, 1), [120.0], [20.0])
fit = newton_fit(d, [0.0])
fit.θ[1] ≈ log(20 / 100)                              # true
fit.se[1] ≈ sqrt(1 / (120 * (1 / 6) * (5 / 6)))        # true

# The same fit from 120 binary rows, and from one weighted row
X = ones(120, 1); y = [i <= 20 for i in 1:120]
newton_fit(logistic_derivatives(X, y), [0.0]).θ[1] ≈ log(20 / 100)      # true
```
"""
function logistic_derivatives(X::AbstractMatrix{<:Real}, y::AbstractVector{Bool};
                              offset::Union{Nothing,AbstractVector{<:Real}}=nothing,
                              weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = size(X, 1)
    length(y) == n ||
        throw(ArgumentError("logistic_derivatives: y has $(length(y)) entries but X has $n rows"))
    if weights === nothing
        n_tot = ones(n)
    else
        length(weights) == n ||
            throw(ArgumentError("logistic_derivatives: weights has $(length(weights)) " *
                                "entries but X has $n rows"))
        n_tot = Vector{Float64}(weights)
        all(w -> isfinite(w) && w >= 0, n_tot) ||
            throw(ArgumentError("logistic_derivatives: weights must be finite and non-negative"))
    end
    n_one = Vector{Float64}(undef, n)
    @inbounds for r in 1:n
        n_one[r] = y[r] ? n_tot[r] : 0.0
    end
    return _logistic_closure(_design(X), n_tot, n_one, _offset_vector(offset, n))
end

function logistic_derivatives(X::AbstractMatrix{<:Real}, n_tot::AbstractVector{<:Real},
                              n_one::AbstractVector{<:Real};
                              offset::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = size(X, 1)
    length(n_tot) == n ||
        throw(ArgumentError("logistic_derivatives: n_tot has $(length(n_tot)) entries but X has $n rows"))
    length(n_one) == n ||
        throw(ArgumentError("logistic_derivatives: n_one has $(length(n_one)) entries but X has $n rows"))
    tot = Vector{Float64}(n_tot)
    one = Vector{Float64}(n_one)
    @inbounds for r in 1:n
        (isfinite(tot[r]) && tot[r] >= 0) ||
            throw(ArgumentError("logistic_derivatives: n_tot must be finite and non-negative " *
                                "(row $r has $(tot[r]))"))
        (isfinite(one[r]) && 0 <= one[r] <= tot[r]) ||
            throw(ArgumentError("logistic_derivatives: n_one must lie in [0, n_tot] " *
                                "(row $r has n_one = $(one[r]), n_tot = $(tot[r]))"))
    end
    return _logistic_closure(_design(X), tot, one, _offset_vector(offset, n))
end

# Keep a Float64 design as is (no copy); convert anything else once.
_design(X::AbstractMatrix{Float64}) = X
_design(X::AbstractMatrix{<:Real}) = Matrix{Float64}(X)

_offset_vector(::Nothing, n::Int) = nothing
function _offset_vector(offset::AbstractVector{<:Real}, n::Int)
    length(offset) == n ||
        throw(ArgumentError("logistic_derivatives: offset has $(length(offset)) entries but X has $n rows"))
    return Vector{Float64}(offset)
end

# The ONE kernel. `off` is `nothing` or a `Vector{Float64}`; this function is
# specialised on that type, so the closure captures it concretely and the
# per-row branch folds away.
function _logistic_closure(X::AbstractMatrix{Float64}, n_tot::Vector{Float64},
                           n_one::Vector{Float64}, off::Union{Nothing,Vector{Float64}})
    n, p = size(X)
    η = Vector{Float64}(undef, n)
    resid = Vector{Float64}(undef, n)
    WX = Matrix{Float64}(undef, n, p)

    return function (β)
        n == 0 && return (0.0, zeros(p), zeros(p, p))
        mul!(η, X, β)
        ll = 0.0
        @inbounds for r in 1:n
            ηr = off === nothing ? η[r] : η[r] + off[r]
            pr = 1.0 / (1.0 + exp(-ηr))
            # log(1 + e^η) on its stable branch
            l1pe = ηr > 0 ? ηr + log1p(exp(-ηr)) : log1p(exp(ηr))
            ll += n_one[r] * ηr - n_tot[r] * l1pe
            resid[r] = n_one[r] - n_tot[r] * pr
            w = n_tot[r] * pr * (1 - pr)
            for k in 1:p
                WX[r, k] = w * X[r, k]
            end
        end
        grad = Vector{Float64}(undef, p)
        mul!(grad, transpose(X), resid)
        hess = Matrix{Float64}(undef, p, p)
        mul!(hess, transpose(X), WX, -1.0, 0.0)
        return ll, grad, hess
    end
end
