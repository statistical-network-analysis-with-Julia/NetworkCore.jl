# The ecosystem's ONE separation verdict, and the ONE policy that goes with it.
#
# A log-likelihood that is concave in θ has no finite maximiser exactly when it
# has a *direction of recession*: a direction d along which it never decreases
# and somewhere strictly increases. For the three likelihood families the
# ecosystem fits, that condition is a set of linear inequalities on d — one
# "margin" row per observation, which d must not make negative:
#
#   logistic   a tie / success row x needs x'd ≥ 0, a non-tie row −x'd ≥ 0
#              (a binomial row with both successes and failures needs both,
#              i.e. x'd = 0);
#   Poisson    every row needs −x'd ≥ 0, and a row with a positive count
#              also x'd ≥ 0 (a zero count is what lets a rate run to 0);
#   clogit     in every stratum each chosen alternative c needs
#              (x_c − x_j)'d ≥ 0 against every alternative j.
#
# Separation is "some d makes every margin ≥ 0 and at least one > 0". A d that
# makes every margin exactly 0 is not separation but non-identification
# (collinear columns), which the optimizer's rank check reports instead.
#
# The decision is a linear programme, R ergm's `mple.existence` and the
# `detectseparation` package solve one as well. NetworkCore stays free of an
# LP dependency, so it carries a small revised simplex on the dual of
#
#     maximise 1'M u  subject to  M u ≥ 0,  −1 ≤ u ≤ 1,
#
# whose dual has only q (= number of coefficients) equality rows, so a basis
# is q × q whatever the number of observations. The optimum is positive iff
# the design is separated (Stiemke's lemma: otherwise some strictly positive
# weights make the margin rows sum to zero, which is the statement that the
# observed sufficient statistic lies in the relative interior of its convex
# support, i.e. the MLE exists). The optimal simplex multipliers are the
# direction. The direction is reduced to an inclusion-minimal set of
# coefficients in floating point, and only then CERTIFIED: the margins that
# are numerically zero are made exactly zero within the exact null space of
# those rows, and every margin's sign is proved (by a rigorous floating-point
# bound where it is far from zero, by exact integer arithmetic where it is
# not), so a certified verdict is a proof, not a tolerance. The certificate
# works on the reduced support, where most rows are zero, so its cost stays
# proportional to the floating-point programme.

using LinearAlgebra: I, norm

"""
    SeparationVerdict

The answer of [`logistic_separation`](@ref), [`poisson_separation`](@ref),
[`clogit_separation`](@ref) and [`separation_from_margins`](@ref): does the
likelihood have a finite maximiser, and if not, which coefficients run off
to infinity?

Fields:
- `separated::Bool`: `true` when the log-likelihood has a direction of
  recession, so no finite maximum-likelihood (or pseudo-likelihood) estimate
  exists.
- `terms::Vector{Int}`: the coefficients to flag: an inclusion-minimal set
  that carries a direction of recession on its own (the other coefficients
  held fixed), so no proper subset of it does. A design can have more than
  one such set (two covariates that each separate the data); the verdict
  names one, and a refit without it is checked again. Empty when not
  separated.
- `direction::Vector{Float64}`: a direction `d`, non-zero only on `terms`,
  along which the log-likelihood never decreases and somewhere increases,
  scaled to `maximum(abs, d) == 1`. Empty when not separated.
- `units::Vector{Int}`: the observations predicted perfectly along
  `direction`: design rows for the logistic and Poisson families, stratum
  labels for the conditional logit, margin rows for
  [`separation_from_margins`](@ref).
- `certified::Bool`: `true` when `direction` was verified in exact
  arithmetic, so the verdict is a proof. A separated verdict that could not
  be certified (it rests on the linear programme's floating-point optimum
  alone) has `certified == false`. That happens when the design is separated
  only up to the rounding of its entries, or when proving it would need more
  than 200,000 rows evaluated in big-integer arithmetic (a cap that keeps
  the verdict's cost bounded); neither has been observed on a fitter's
  design.
- `family::Symbol`: `:logistic`, `:poisson`, `:clogit` or `:margins`.

`verdict.separated` is the whole policy trigger; see
[`warn_separation`](@ref) for what a fitting function does with it.

# Example
```julia
using NetworkCore
x = collect(-3.0:3.0)
v = NetworkCore.logistic_separation(hcat(ones(7), x), x .> 0)
v.separated, v.certified      # (true, true)
v.terms                       # [2]: the slope alone separates the data
```
"""
struct SeparationVerdict
    separated::Bool
    terms::Vector{Int}
    direction::Vector{Float64}
    units::Vector{Int}
    certified::Bool
    family::Symbol
end

function Base.show(io::IO, v::SeparationVerdict)
    if v.separated
        print(io, "SeparationVerdict(separated, family=:", v.family, ", terms=", v.terms,
              ", ", length(v.units), " unit", length(v.units) == 1 ? "" : "s",
              " predicted perfectly", v.certified ? ", certified)" : ", not certified)")
    else
        print(io, "SeparationVerdict(not separated, family=:", v.family, ")")
    end
end

_not_separated(family::Symbol) =
    SeparationVerdict(false, Int[], Float64[], Int[], false, family)

# ----------------------------------------------------------------------------
# The core criterion
# ----------------------------------------------------------------------------

"""
    separation_from_margins(M::AbstractMatrix{<:Real}; family=:margins) -> SeparationVerdict

The core criterion every family entry point reduces to. Each row `m_k` of
`M` is a linear form the direction must keep non-negative: the design is
**separated** iff some `d` has `M d ≥ 0` with at least one strictly
positive entry. A `d` with `M d == 0` only (a collinear design) is not
separation.

The decision is the linear programme `max 1'M u` subject to `M u ≥ 0` and
`-1 ≤ u ≤ 1`, solved by a revised simplex on its q-row dual (q =
`size(M, 2)`), after the rows and columns are equilibrated. Its optimum is
positive iff the design is separated. The direction it returns is reduced to
an inclusion-minimal set of columns and then certified exactly: the margins
that are numerically zero are set exactly to zero inside their exact null
space, and the sign of every margin is proved, by a rigorous floating-point
bound when it is far from zero and in exact integer arithmetic otherwise.
The cost is that of a few floating-point programmes: a REM design of 640,000
margin rows takes about a second whether or not it is separated. `units`
lists the rows of `M` with a strictly positive margin. A row of zeros
constrains nothing and is ignored. Non-finite entries are an
`ArgumentError`.

Use it directly for a likelihood that is none of the three families, after
writing down its margin rows; otherwise call the family entry point, which
builds them.

# Example
```julia
using NetworkCore
# Two forms a direction must keep non-negative: d₁ ≥ 0 and d₂ − d₁ ≥ 0.
v = NetworkCore.separation_from_margins([1.0 0.0; -1.0 1.0])
v.separated                        # true: d = (0, 1) makes the second one positive
# d₁ ≥ 0 and −d₁ ≥ 0 force d₁ = 0: an equality, nothing can grow
NetworkCore.separation_from_margins([1.0; -1.0;;]).separated     # false
```
"""
function separation_from_margins(M::AbstractMatrix{<:Real}; family::Symbol=:margins)
    A = M isa Matrix{Float64} ? M : Matrix{Float64}(M)
    all(isfinite, A) || throw(ArgumentError(
        "separation_from_margins: the margin rows must be finite " *
        "(drop rows fixed by an infinite offset before testing)"))
    m, q = size(A)
    q == 0 && return _not_separated(family)
    keep = [k for k in 1:m if _nonzero_row(A, k)]
    isempty(keep) && return _not_separated(family)
    Ak = length(keep) == m ? A : A[keep, :]
    allcols = collect(1:q)

    # 1. The linear programme on the full design, in floating point.
    Sbuf = Matrix{Float64}(undef, size(Ak))
    lp = _lp_direction(Ak, allcols, Sbuf)
    (lp === nothing || lp.value <= 1e-10) && return _not_separated(family)
    full_cert = nothing
    if lp.value <= _LP_SEPARATION_VALUE
        # An optimum at the edge of rounding: only an exact certificate on the
        # full design can establish separation.
        full_cert = _certify_ray(Ak, allcols, lp.u)
        full_cert === nothing && return _not_separated(family)
    end

    # 2. Reduce the direction to an inclusion-minimal support, in floating
    # point: a coefficient is dropped (its direction component fixed at 0,
    # i.e. its column deleted) whenever the remaining ones still carry a
    # direction of recession. The least involved go first. Without this, a
    # design separated by one statistic reports every other statistic too,
    # since the programme's optimum also moves coefficients that merely may
    # move.
    support = allcols
    best = lp
    for j in sortperm(abs.(lp.u .* lp.colscale))
        length(support) == 1 && break
        trial = filter(!=(j), support)
        sub = _lp_direction(Ak, trial, Sbuf)
        (sub === nothing || sub.value <= _LP_SEPARATION_VALUE) && continue
        support = trial
        best = sub
    end

    # 3. Certify the reduced direction exactly. The reduced support is small
    # and most rows are zero on it, so this costs little; the full-support
    # certificate is the fallback when the reduced one cannot be found.
    cert = nothing
    if length(support) < q
        cert = _certify_ray(Ak, support, best.u)
        if cert === nothing
            full_cert === nothing && (full_cert = _certify_ray(Ak, allcols, lp.u))
            if full_cert !== nothing
                support, best, cert = allcols, lp, full_cert
            end
        end
    else
        cert = full_cert === nothing ? _certify_ray(Ak, allcols, lp.u) : full_cert
    end

    u = zeros(q)
    u[support] = cert === nothing ? best.u : Float64.(cert.ray)
    u ./= maximum(abs, u)
    positive = cert === nothing ?
        findall(k -> ((mf, mag) = _float_margin(Ak, k, support, best.u); mf > 1e-9 * mag),
                axes(Ak, 1)) :
        cert.positive
    return SeparationVerdict(true, support, u, keep[positive], cert !== nothing, family)
end

function _nonzero_row(A::Matrix{Float64}, k::Int)
    @inbounds for j in axes(A, 2)
        iszero(A[k, j]) || return true
    end
    return false
end

# The margin of row k on the columns `cols` along `u`, and the sum of the
# magnitudes of its terms, in floating point.
@inline function _float_margin(A::Matrix{Float64}, k::Int, cols::Vector{Int}, u::Vector{Float64})
    mf = 0.0; mag = 0.0
    @inbounds for jj in eachindex(cols)
        t = A[k, cols[jj]] * u[jj]
        mf += t; mag += abs(t)
    end
    return mf, mag
end

# The linear programme on the columns `cols` of the non-zero rows `A`:
# `nothing` when no row is non-zero on them, otherwise the floating-point
# optimum `u` (original coordinates of those columns), the column scales and
# the programme's value. Positive column scales are a change of variables and
# positive row scales leave each inequality unchanged, so the cone is the
# same and the tolerances become relative. The scaled matrix is written into
# `Sbuf` (size(A) or larger), one buffer for every programme of a verdict, so
# the programmes on subsets of columns allocate O(m), not a copy of the design.
function _lp_direction(A::Matrix{Float64}, cols::Vector{Int}, Sbuf::Matrix{Float64})
    q = length(cols)
    q == 0 && return nothing
    m = size(A, 1)
    colscale = ones(q)
    @inbounds for (jj, j) in enumerate(cols)
        s = 0.0
        for k in 1:m
            s = max(s, abs(A[k, j]))
        end
        s > 0 && (colscale[jj] = s)
    end
    rowscale = zeros(m)
    @inbounds for (jj, j) in enumerate(cols), k in 1:m
        rowscale[k] = max(rowscale[k], abs(A[k, j]) / colscale[jj])
    end
    keep = findall(>(0), rowscale)
    isempty(keep) && return nothing
    S = view(Sbuf, 1:length(keep), 1:q)
    @inbounds for (jj, j) in enumerate(cols), (i, k) in enumerate(keep)
        S[i, jj] = A[k, j] / colscale[jj] / rowscale[k]
    end
    ũ, value = _recession_lp(S)
    return (u=ũ ./ colscale, value=value, colscale=colscale)
end

# A floating-point LP optimum above this (in equilibrated units: every row has
# max-abs 1 and |u| ≤ 1, so a genuine separation has an optimum of order the
# smallest positive margin) counts as separation even when the direction
# cannot be certified exactly.
const _LP_SEPARATION_VALUE = 1e-6

# Revised simplex on the dual of  max 1'S u  s.t.  S u ≥ 0, −1 ≤ u ≤ 1.
# Dual in standard form: variables (y ∈ ℝᵐ, v⁺, v⁻ ∈ ℝ^q) ≥ 0,
#     −S'y + v⁺ − v⁻ = c := S'1,   minimise 1'v⁺ + 1'v⁻.
# Only q equality rows, so the basis is q × q. At the optimum the simplex
# multipliers π satisfy S π ≥ 0 and |π| ≤ 1 (the reduced costs), i.e. π is
# the primal optimum, and c'π its value. Dantzig pricing, switching to
# Bland's rule after a run of degenerate pivots (the programme is degenerate
# at its apex, u = 0), which cannot cycle.
function _recession_lp(S::AbstractMatrix{Float64})
    m, q = size(S)
    c = vec(sum(S; dims=1))
    basis = [c[j] >= 0 ? m + j : m + q + j for j in 1:q]
    Binv = Matrix{Float64}(I, q, q)
    for j in 1:q
        c[j] < 0 && (Binv[j, j] = -1.0)
    end
    xB = abs.(c)
    cost(k) = k <= m ? 0.0 : 1.0
    function column!(a, k)
        if k <= m
            @inbounds for j in 1:q
                a[j] = -S[k, j]
            end
        else
            fill!(a, 0.0)
            k <= m + q ? (a[k - m] = 1.0) : (a[k - m - q] = -1.0)
        end
        return a
    end
    π = zeros(q); r = zeros(m); a = zeros(q); d = zeros(q); cB = zeros(q)
    # reduced cost of variable k (rows, then v⁺, then v⁻: ONE index order,
    # used by both halves of Bland's rule)
    reduced(k) = k <= m ? r[k] : (k <= m + q ? 1.0 - π[k - m] : 1.0 + π[k - m - q])
    degenerate_run = 0
    tol_opt = 1e-10
    tol_piv = 1e-11
    maxiter = 50_000 + 50 * (m + q)
    optimal = false
    for iter in 1:maxiter
        for i in 1:q
            cB[i] = cost(basis[i])
        end
        mul!(π, transpose(Binv), cB)
        mul!(r, S, π)
        bland = degenerate_run >= 20
        entering = 0
        best = -tol_opt
        for k in 1:(m + 2q)
            rc = reduced(k)
            rc < best || continue
            k in basis && continue
            entering = k
            bland && break                       # Bland: the smallest index
            best = rc                            # Dantzig: the most negative
        end
        if entering == 0
            optimal = true
            break
        end
        column!(a, entering)
        mul!(d, Binv, a)
        leave = 0
        θ = Inf
        for i in 1:q
            d[i] > tol_piv || continue
            t = max(xB[i], 0.0) / d[i]
            if leave == 0 || t < θ - 1e-14
                θ = t; leave = i
            elseif abs(t - θ) <= 1e-14
                # ties: Bland takes the smallest variable index; Dantzig the
                # largest pivot, for stability
                if bland ? basis[i] < basis[leave] : d[i] > d[leave]
                    θ = t; leave = i
                end
            end
        end
        leave == 0 && break     # unbounded: impossible, the objective is ≥ 0
        piv = d[leave]
        for j in 1:q
            Binv[leave, j] /= piv
        end
        xB[leave] /= piv
        for i in 1:q
            i == leave && continue
            f = d[i]
            f == 0 && continue
            for j in 1:q
                Binv[i, j] -= f * Binv[leave, j]
            end
            xB[i] -= f * xB[leave]
        end
        basis[leave] = entering
        degenerate_run = θ <= 1e-14 ? degenerate_run + 1 : 0
        if iter % 32 == 0                        # refactorise against drift
            B = zeros(q, q)
            for i in 1:q
                B[:, i] = column!(a, basis[i])
            end
            # in place: the closures above capture these arrays, and
            # rebinding a captured variable would box it
            Binv .= inv(B)
            mul!(xB, Binv, c)
            xB .= max.(xB, 0.0)
        end
    end
    optimal || error("separation_from_margins: the simplex stopped without an optimum " *
                     "(m = $m rows, q = $q columns); please report this design")
    for i in 1:q
        cB[i] = cost(basis[i])
    end
    mul!(π, transpose(Binv), cB)
    clamp!(π, -1.0, 1.0)
    return π, max(dot(c, π), 0.0)
end

# ----------------------------------------------------------------------------
# The exact certificate
# ----------------------------------------------------------------------------
#
# A certificate is a rational ray d on the columns `cols` with every margin
# A[k, cols]·d ≥ 0 and at least one > 0, decided exactly. The floating-point
# direction is corrected first: margins within a relative threshold of zero
# are made exactly zero inside the exact null space of those rows (forcing an
# extra equality is always safe, since the result is verified). Thresholds are
# tried from tight to loose; the first ray that verifies wins.
#
# The cost is kept proportional to the floating-point programme:
#   - a row that is zero on `cols` has margin exactly 0 and is never touched;
#   - the null space is that of a BASIS of the near-zero rows, at most
#     length(cols) of them, picked in floating point (Gram–Schmidt on the
#     equilibrated rows) and eliminated exactly; a row the floating-point
#     test misjudged shows up in the verification and is added;
#   - a margin is decided in floating point when it is far from zero relative
#     to its terms (a rigorous bound: the rounding error of a dot product of
#     ≤ 10⁵ terms with the rounded ray is below 1e-10 of their magnitude), and
#     exactly otherwise, by an integer dot product of the Float64 entries'
#     significands with the ray scaled to integers: Int128 with overflow
#     checks, BigInt only when that overflows;
#   - the free components of the ray are first tried as small rationals near
#     the floating-point ones (an LP vertex of an integer design is a ratio of
#     small integers), which keeps the scaled ray in Int64.
# BigInt evaluations are capped at `_EXACT_BUDGET` rows per certificate; past
# the cap the certificate is abandoned and the verdict says `certified=false`.

const _EXACT_BUDGET = 200_000

# `nothing`, or `(ray=Vector{Rational{BigInt}}, positive=Vector{Int})` with the
# rows of A whose margin is strictly positive.
function _certify_ray(A::Matrix{Float64}, cols::Vector{Int}, u::Vector{Float64})
    all(isfinite, u) || return nothing
    um = maximum(abs, u; init=0.0)
    um > 0 || return nothing
    u = u ./ um
    s = length(cols)
    m = size(A, 1)
    # rows non-zero on the support, with their floating-point margins
    rows = Int[]
    margins = Float64[]
    mags = Float64[]
    for k in 1:m
        nz = false
        @inbounds for jj in 1:s
            iszero(A[k, cols[jj]]) || (nz = true; break)
        end
        nz || continue
        mf, mag = _float_margin(A, k, cols, u)
        push!(rows, k); push!(margins, mf); push!(mags, mag)
    end
    isempty(rows) && return nothing
    colscale = zeros(s)
    for (jj, j) in enumerate(cols), k in rows
        colscale[jj] = max(colscale[jj], abs(A[k, j]))
    end
    replace!(colscale, 0.0 => 1.0)
    budget = Ref(_EXACT_BUDGET)
    tried = Set{Vector{Int}}()
    for rel in (64 * eps(Float64), 1e-12, 1e-9, 1e-6)
        # a clearly negative margin beyond the threshold rules this level out
        any(i -> margins[i] < -rel * mags[i], eachindex(rows)) && continue
        boundary = [rows[i] for i in eachindex(rows) if abs(margins[i]) <= rel * mags[i]]
        boundary in tried && continue          # the same rows as a tighter level
        push!(tried, boundary)
        for snap in (true, false)
            cert = _certify_with(A, cols, u, rows, boundary, colscale, snap, budget)
            cert === nothing || return cert
            budget[] <= 0 && return nothing
        end
    end
    return nothing
end

# One attempt: the null correction for `boundary`, verified on `rows`.
function _certify_with(A, cols, u, rows, boundary, colscale, snap::Bool, budget)
    Q = Rational{BigInt}
    s = length(cols)
    # the free components: the floating-point ones, or small rationals near them
    base = snap ? Q[_snap_rational(x) for x in u] : Q.(u)
    basis = _BoundaryBasis(s)
    for k in boundary
        _try_add_row!(basis, A, k, cols, colscale) || continue
        basis.rank == s && return nothing          # full rank: only the zero ray
    end
    for _ in 0:s
        ray = _back_substitute(basis, base)
        all(iszero, ray) && return nothing
        N, N64 = _integer_ray(ray)
        uf = Float64.(ray ./ maximum(abs, ray))
        filter_ok = all(i -> iszero(ray[i]) || abs(uf[i]) >= floatmin(Float64), eachindex(ray))
        positive = Int[]
        violated = 0
        for k in rows
            sg = _margin_sign(A, k, cols, uf, filter_ok, N, N64, budget)
            sg === nothing && return nothing       # past the exact-arithmetic budget
            if sg < 0
                violated = k
                break
            end
            sg > 0 && push!(positive, k)
        end
        if violated == 0
            isempty(positive) && return nothing
            return (ray=ray, positive=positive)
        end
        # A boundary row the floating-point basis judged dependent: add it
        # exactly and correct again. Any other negative row ends the attempt.
        violated in boundary || return nothing
        _add_row_exact!(basis, A, violated, cols) || return nothing
        basis.rank == s && return nothing
    end
    return nothing
end

# Small rational near x (an LP vertex of an integer design is a ratio of small
# integers); x itself when none is within 1e-12.
function _snap_rational(x::Float64)
    r = rationalize(Int64, x; tol=1e-12)
    return Rational{BigInt}(denominator(r) <= 1 << 20 ? r : x)
end

# Row echelon form of the rows forced to margin zero: exact rows and pivots,
# plus an orthonormal floating-point basis of the same (equilibrated) span,
# used only to skip the boundary rows already in it. (A row added after the
# verification found it violated needs no floating-point column: every
# boundary row has been offered by then.)
mutable struct _BoundaryBasis
    rows::Vector{Vector{Rational{BigInt}}}
    pivots::Vector{Int}
    F::Matrix{Float64}
    rank::Int
    v::Vector{Float64}
end
_BoundaryBasis(s::Int) = _BoundaryBasis(Vector{Rational{BigInt}}[], Int[], zeros(s, s), 0, zeros(s))

# Add row k when the floating-point test finds it outside the span (and exact
# elimination confirms it); `true` if added.
function _try_add_row!(b::_BoundaryBasis, A, k, cols, colscale)
    s = length(cols)
    v = b.v
    for jj in 1:s
        v[jj] = A[k, cols[jj]] / colscale[jj]
    end
    nv = norm(v)
    nv == 0 && return false
    v ./= nv
    for _ in 1:2, i in 1:b.rank
        f = view(b.F, :, i)
        v .-= dot(f, v) .* f
    end
    res = norm(v)
    res > 1e-8 || return false
    _add_row_exact!(b, A, k, cols) || return false
    b.F[:, b.rank] .= v ./ res
    return true
end

# Exact elimination of row k against the basis; `true` if it adds a pivot.
function _add_row_exact!(b::_BoundaryBasis, A, k, cols)
    Q = Rational{BigInt}
    s = length(cols)
    row = Q[A[k, cols[jj]] for jj in 1:s]
    for (bi, pivot) in zip(b.rows, b.pivots)
        iszero(row[pivot]) && continue
        factor = row[pivot] / bi[pivot]
        for j in 1:s
            iszero(bi[j]) || (row[j] -= factor * bi[j])
        end
    end
    pivot = findfirst(!iszero, row)
    pivot === nothing && return false
    push!(b.rows, row)
    push!(b.pivots, pivot)
    b.rank += 1
    return true
end

# The ray: the free components from `base`, the pivot components solved so
# that every basis row has margin exactly zero.
function _back_substitute(b::_BoundaryBasis, base::Vector{Rational{BigInt}})
    ray = copy(base)
    k = length(ray)
    for i in reverse(eachindex(b.rows))
        row, pivot = b.rows[i], b.pivots[i]
        acc = zero(Rational{BigInt})
        for j in 1:k
            (j == pivot || iszero(row[j]) || iszero(ray[j])) && continue
            acc += row[j] * ray[j]
        end
        ray[pivot] = -acc / row[pivot]
    end
    return ray
end

# The ray as coprime integers (same direction), and as Int64 when they fit.
function _integer_ray(ray::Vector{Rational{BigInt}})
    den = foldl(lcm, (denominator(x) for x in ray); init=one(BigInt))
    N = [numerator(x) * (den ÷ denominator(x)) for x in ray]
    g = foldl(gcd, N; init=zero(BigInt))
    g > 1 && (N .= N .÷ g)
    N64 = all(x -> typemin(Int64) < x <= typemax(Int64), N) ? Int64.(N) : nothing
    return N, N64
end

# A non-zero finite Float64 as an odd integer times a power of two, exactly.
@inline function _dyadic(a::Float64)
    e = exponent(a)
    mant = Int64(ldexp(a, 52 - e))
    tz = trailing_zeros(mant)
    return mant >> tz, e - 52 + tz
end

# The sign of row k's margin along the ray: in floating point when it is
# clearly away from zero, otherwise exactly. `nothing` past the budget.
function _margin_sign(A, k, cols, uf, filter_ok, N, N64, budget)
    if filter_ok
        mf, mag = _float_margin(A, k, cols, uf)
        if mag > 1e-290 && abs(mf) > 1e-10 * mag
            return mf > 0 ? 1 : -1
        end
    end
    if N64 !== nothing
        sg = _exact_sign_int128(A, k, cols, N64)
        sg === nothing || return sg
    end
    budget[] -= 1
    budget[] < 0 && return nothing
    return _exact_sign_bigint(A, k, cols, N)
end

# Σ_jj A[k, cols[jj]] * N[jj] exactly in Int128; `nothing` on overflow.
function _exact_sign_int128(A, k, cols, N64::Vector{Int64})
    emin = typemax(Int)
    @inbounds for jj in eachindex(cols)
        a = A[k, cols[jj]]
        (iszero(a) || iszero(N64[jj])) && continue
        emin = min(emin, last(_dyadic(a)))
    end
    emin == typemax(Int) && return 0
    acc = Int128(0)
    @inbounds for jj in eachindex(cols)
        a = A[k, cols[jj]]
        (iszero(a) || iszero(N64[jj])) && continue
        mant, ex = _dyadic(a)
        t = Int128(mant) * Int128(N64[jj])      # |t| < 2^116: no overflow
        sh = ex - emin
        sh >= 0 || return nothing
        (sh > 0 && leading_zeros(abs(t)) <= sh + 1) && return nothing
        t <<= sh
        acc, ovf = Base.Checked.add_with_overflow(acc, t)
        ovf && return nothing
    end
    return Int(sign(acc))
end

function _exact_sign_bigint(A, k, cols, N::Vector{BigInt})
    emin = typemax(Int)
    for jj in eachindex(cols)
        a = A[k, cols[jj]]
        (iszero(a) || iszero(N[jj])) && continue
        emin = min(emin, last(_dyadic(a)))
    end
    emin == typemax(Int) && return 0
    acc = zero(BigInt)
    for jj in eachindex(cols)
        a = A[k, cols[jj]]
        (iszero(a) || iszero(N[jj])) && continue
        mant, ex = _dyadic(a)
        acc += (BigInt(mant) * N[jj]) << (ex - emin)
    end
    return Int(sign(acc))
end

# Map a margin-level verdict back to the family's observation units:
# `owner[k]` is the unit of margin row k.
function _family_verdict(M::Matrix{Float64}, owner::Vector{Int}, family::Symbol)
    v = separation_from_margins(M; family=family)
    v.separated || return v
    units = sort!(unique(owner[v.units]))
    return SeparationVerdict(true, v.terms, v.direction, units, v.certified, family)
end

_check_design(X, n, what) = size(X, 1) == n || throw(ArgumentError(
    "$what has $n entries but the design has $(size(X, 1)) rows"))

# ----------------------------------------------------------------------------
# The three likelihood families
# ----------------------------------------------------------------------------

"""
    logistic_separation(X, y::AbstractVector{Bool}; weights=nothing) -> SeparationVerdict
    logistic_separation(X, n_tot::AbstractVector{<:Real}, n_one::AbstractVector{<:Real}) -> SeparationVerdict

Whether the logistic (pseudo-)likelihood of the binary responses `y` (or of
binomial rows, `n_one[r]` successes in `n_tot[r]` trials: ERGM's compressed
MPLE design) on design `X` has a finite maximiser. Separation means a
direction `d` with `x_r'd ≥ 0` on every success row and `x_r'd ≤ 0` on every
failure row, strictly somewhere: complete or quasi-complete separation in the
sense of Albert & Anderson (1984), the condition R ergm's `mple.existence`
tests. A binomial row with both successes and failures forces `x_r'd = 0`.
Rows with zero weight or `n_tot[r] == 0` carry no information and are
skipped. An offset does not change the answer (a finite shift of the linear
predictor); drop rows fixed by an infinite offset first. `units` are the
design rows predicted perfectly.

The verdict depends on the data alone, not on where an optimizer stopped:
call it before or after the fit. See [`separation_from_margins`](@ref) for
the criterion and [`warn_separation`](@ref) for the policy.

# Example
```julia
using NetworkCore
X = hcat(ones(6), [-2.0, -1, 0, 0, 1, 2])
NetworkCore.logistic_separation(X, Bool[0, 0, 0, 1, 1, 1]).separated   # true: quasi-complete
NetworkCore.logistic_separation(X, Bool[0, 1, 0, 1, 0, 1]).separated   # false
# the binomial-row form: a nodematch-like column with no ties in its cell
Xb = [1.0 0.0; 1.0 1.0]
NetworkCore.logistic_separation(Xb, [40.0, 10.0], [8.0, 0.0]).terms    # [2]
```
"""
function logistic_separation(X::AbstractMatrix{<:Real}, y::AbstractVector{Bool};
                             weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = length(y)
    _check_design(X, n, "y")
    weights === nothing || _check_design(X, length(weights), "weights")
    rows = [r for r in 1:n if weights === nothing || weights[r] > 0]
    M = Matrix{Float64}(undef, length(rows), size(X, 2))
    for (k, r) in enumerate(rows)
        s = y[r] ? 1.0 : -1.0
        for j in axes(X, 2)
            M[k, j] = s * X[r, j]
        end
    end
    return _family_verdict(M, rows, :logistic)
end

function logistic_separation(X::AbstractMatrix{<:Real}, n_tot::AbstractVector{<:Real},
                             n_one::AbstractVector{<:Real})
    n = length(n_tot)
    _check_design(X, n, "n_tot")
    _check_design(X, length(n_one), "n_one")
    all(r -> 0 <= n_one[r] <= n_tot[r], 1:n) || throw(ArgumentError(
        "logistic_separation: n_one must lie in [0, n_tot] on every row"))
    owner = Int[]
    signs = Float64[]
    for r in 1:n
        n_tot[r] > 0 || continue
        n_one[r] > 0 && (push!(owner, r); push!(signs, 1.0))
        n_one[r] < n_tot[r] && (push!(owner, r); push!(signs, -1.0))
    end
    M = Matrix{Float64}(undef, length(owner), size(X, 2))
    for k in eachindex(owner), j in axes(X, 2)
        M[k, j] = signs[k] * X[owner[k], j]
    end
    return _family_verdict(M, owner, :logistic)
end

"""
    poisson_separation(X, y::AbstractVector{<:Real}; weights=nothing) -> SeparationVerdict

Whether the Poisson log-linear likelihood `Σ_r y_r η_r − exp(η_r)`, `η = Xθ`
(plus any finite offset), has a finite maximiser. Along `d` the term
`exp(η_r)` grows without bound unless `x_r'd ≤ 0`, and a positive count
`y_r` stops `η_r` from falling, so the condition is: `x_r'd ≤ 0` on every
row, `x_r'd = 0` on every row with `y_r > 0`, and `x_r'd < 0` somewhere.
That is the zero-count boundary (Haberman 1974): a rate whose every
observation is zero has its maximum at a log-rate of `-Inf`. Counts must be
non-negative; zero-weight rows are skipped. `units` are the design rows
whose fitted rate goes to zero.

The ERGM family's count models condition on each dyad's support instead of
fitting a Poisson GLM: their likelihood is a conditional logit over the
support values, so they call [`clogit_separation`](@ref).

# Example
```julia
using NetworkCore
X = [1.0 0.0; 1.0 0.0; 1.0 1.0; 1.0 1.0]
NetworkCore.poisson_separation(X, [3, 1, 0, 0]).terms      # [2]: the x = 1 rate is 0
NetworkCore.poisson_separation(X, [3, 1, 0, 2]).separated  # false
```
"""
function poisson_separation(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real};
                            weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = length(y)
    _check_design(X, n, "y")
    weights === nothing || _check_design(X, length(weights), "weights")
    all(v -> isfinite(v) && v >= 0, y) || throw(ArgumentError(
        "poisson_separation: counts must be finite and non-negative"))
    owner = Int[]
    signs = Float64[]
    for r in 1:n
        (weights === nothing || weights[r] > 0) || continue
        push!(owner, r); push!(signs, -1.0)
        y[r] > 0 && (push!(owner, r); push!(signs, 1.0))
    end
    M = Matrix{Float64}(undef, length(owner), size(X, 2))
    for k in eachindex(owner), j in axes(X, 2)
        M[k, j] = signs[k] * X[owner[k], j]
    end
    return _family_verdict(M, owner, :poisson)
end

"""
    clogit_separation(X, chosen::AbstractVector{Bool}, strata::AbstractVector{<:Integer};
                      weights=nothing) -> SeparationVerdict

Whether a conditional-logit (multinomial-choice, Cox partial) likelihood has a
finite maximiser. Row `i` of `X` is one alternative of stratum `strata[i]`;
`chosen[i]` marks the alternatives that were observed (the case: one per
stratum for a relational event, any number for a stratum with case weights,
such as the support values observed in one dyad class of a count ERGM).
Rows with `weights[i] == 0` are not in the stratum's choice set (an
alternative with zero weight in the denominator).

Along `d` a stratum's contribution stays bounded below only if every chosen
alternative scores at least as high as every alternative,
`(x_c − x_j)'d ≥ 0`; separation is that for every stratum, with a strict
inequality somewhere. Two chosen alternatives in one stratum must therefore
tie on `d`. This is what `survival::coxph` diagnoses as "coefficient may be
infinite", decided here exactly rather than from the Newton increment.
Strata are those of the risk sets actually used: with sampled controls the
question is about the sampled likelihood, which is the one being fitted.
A stratum without a chosen alternative is an `ArgumentError`. `units` are
the stratum labels in which the case beats some alternative strictly.

# Example
```julia
using NetworkCore
# Three events, each against two controls; one statistic.
X = reshape([1.0, 0.0, 0.0,   2.0, 1.0, 1.0,   1.0, 0.0, 1.0], :, 1)
chosen = Bool[1, 0, 0, 1, 0, 0, 1, 0, 0]
strata = [1, 1, 1, 2, 2, 2, 3, 3, 3]
v = NetworkCore.clogit_separation(X, chosen, strata)
v.separated, v.units          # (true, [1, 2, 3]): the case always scores highest
X[9] = 2.0                    # in event 3 a control now beats the case
NetworkCore.clogit_separation(X, chosen, strata).separated   # false
```
"""
function clogit_separation(X::AbstractMatrix{<:Real}, chosen::AbstractVector{Bool},
                           strata::AbstractVector{<:Integer};
                           weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = length(chosen)
    _check_design(X, n, "chosen")
    _check_design(X, length(strata), "strata")
    weights === nothing || _check_design(X, length(weights), "weights")
    groups = Dict{eltype(strata),Vector{Int}}()
    order = eltype(strata)[]
    for i in 1:n
        (weights === nothing || weights[i] > 0) || continue
        s = strata[i]
        haskey(groups, s) || (groups[s] = Int[]; push!(order, s))
        push!(groups[s], i)
    end
    # Per stratum: ties among the chosen (pairs with the first case, both
    # signs), and the first case against every unchosen alternative.
    nrows = 0
    for s in order
        g = groups[s]
        nc = count(i -> chosen[i], g)
        nc == 0 && throw(ArgumentError(
            "clogit_separation: stratum $s has no chosen alternative"))
        nrows += 2 * (nc - 1) + (length(g) - nc)
    end
    p = size(X, 2)
    M = Matrix{Float64}(undef, nrows, p)
    owner = Vector{Int}(undef, nrows)
    labels = Int[]
    k = 0
    for (si, s) in enumerate(order)
        g = groups[s]
        push!(labels, Int(s))
        c0 = g[findfirst(i -> chosen[i], g)]
        for i in g
            i == c0 && continue
            if chosen[i]
                for sgn in (1.0, -1.0)
                    k += 1
                    for j in 1:p
                        M[k, j] = sgn * (X[c0, j] - X[i, j])
                    end
                    owner[k] = si
                end
            else
                k += 1
                for j in 1:p
                    M[k, j] = X[c0, j] - X[i, j]
                end
                owner[k] = si
            end
        end
    end
    v = _family_verdict(M, owner, :clogit)
    v.separated || return v
    return SeparationVerdict(true, v.terms, v.direction, labels[v.units], v.certified, :clogit)
end

# ----------------------------------------------------------------------------
# The policy: one warning, one caveat, inference withheld
# ----------------------------------------------------------------------------

_unit_noun(family::Symbol) =
    family === :clogit ? ("stratum", "strata") :
    family === :margins ? ("constraint", "constraints") : ("observation", "observations")

function _separated_term_list(v::SeparationVerdict, names)
    isempty(v.terms) && return "a combination of the coefficients"
    labels = names === nothing ? ["coefficient $j" for j in v.terms] :
                                 [string('`', names[j], '`') for j in v.terms]
    return join(labels, ", ")
end

# A pseudo-likelihood estimate (MPLE, CMPLE, swap-MPLE) is not a likelihood:
# name the objective the estimate maximises.
_objective_noun(estimate::AbstractString) =
    occursin(r"MPLE|pseudo"i, estimate) ? "log-pseudo-likelihood" : "log-likelihood"

function _units_phrase(v::SeparationVerdict)
    one, many = _unit_noun(v.family)
    k = length(v.units)
    return "$k $(k == 1 ? one : many)"
end

"""
    separation_message(context, verdict::SeparationVerdict, names=nothing;
                       estimate="maximum-likelihood estimate", note="") -> String

The ONE warning sentence for a separated fit, prefixed with `context` (the
fitting function's name). It names the separated terms (`names[verdict.terms]`,
or their indices when `names` is `nothing`), says how many observations are
predicted perfectly, and states the policy: the coefficients are where the
optimizer stopped, the fit is reported with `converged == false`, and z
values, p-values and confidence intervals are withheld. `estimate` names
what does not exist (`"MPLE"` for a pseudo-likelihood); `note` is appended
as is, for the R function's own wording.

# Example
```julia
using NetworkCore
v = NetworkCore.logistic_separation([1.0 0.0; 1.0 1.0], [40.0, 10.0], [8.0, 0.0])
msg = NetworkCore.separation_message("mple", v, ["edges", "nodematch.g"]; estimate="MPLE")
occursin("`nodematch.g`", msg)     # true
```
"""
function separation_message(context::AbstractString, v::SeparationVerdict, names=nothing;
                            estimate::AbstractString="maximum-likelihood estimate",
                            note::AbstractString="")
    terms = _separated_term_list(v, names)
    return "$context: the $estimate does not exist (separation). The $(_objective_noun(estimate)) " *
           "keeps increasing as the coefficient(s) on $terms run to ±Inf along a " *
           "direction that predicts $(_units_phrase(v)) perfectly, so the returned " *
           "values are where the optimizer stopped, not estimates. The fit is reported " *
           "with `converged == false`, the separated terms are flagged, and z values, " *
           "p-values and confidence intervals are withheld (NaN). Remove, merge or " *
           "coarsen the separating term(s), or collect more data." *
           (isempty(note) ? "" : " " * note)
end

"""
    warn_separation(context, verdict::SeparationVerdict, names=nothing; kwargs...) -> Bool

Apply the ecosystem's separation policy's warning: when `verdict.separated`,
emit [`separation_message`](@ref) as a warning and return `true`; otherwise
do nothing and return `false`. The policy, defined once here and followed by
every fitting function (R's `glm`, `coxph` and ergm's `mple.existence`
behave the same way):

1. **warn** with this message;
2. return the fit with **`converged == false`**;
3. **flag** the separated terms (`verdict.terms`) on the result;
4. **withhold inference**: z values and p-values are `NaN`
   ([`z_pvalues`](@ref) with the verdict does this), and so are confidence
   intervals; [`separation_caveat`](@ref) is the entry for `approximations`.

Bootstrap refits and other internal fits pass on the warning but keep the
verdict.

# Example
```julia
using NetworkCore
x = collect(-3.0:3.0)
v = NetworkCore.logistic_separation(hcat(ones(7), x), x .> 0)
NetworkCore.warn_separation("my_fit", v, ["(Intercept)", "x"])   # warns; true
```
"""
function warn_separation(context::AbstractString, v::SeparationVerdict, names=nothing; kwargs...)
    v.separated || return false
    @warn separation_message(context, v, names; kwargs...)
    return true
end

"""
    separation_caveat(verdict::SeparationVerdict, names=nothing) -> Union{String,Nothing}
    separation_caveat(separated_terms::AbstractVector{<:AbstractString}) -> Union{String,Nothing}

The entry a separated fit lists in `approximations(fit)` (and prints under
its coefficient table): which terms are infinite at the maximum, and that
inference is withheld. `nothing` when the verdict is not separated, so a
caller can write `c = separation_caveat(v, names); c === nothing || push!(out, c)`.
The second form serves a result that stores only the names of its separated
terms (`fit.separated`); it gives the same sentence without the count of
observations, and `nothing` for an empty list.

# Example
```julia
using NetworkCore
v = NetworkCore.poisson_separation([1.0 0.0; 1.0 1.0], [4, 0])
NetworkCore.separation_caveat(v, ["(Intercept)", "x"])
NetworkCore.separation_caveat(NetworkCore.poisson_separation([1.0 0.0; 1.0 1.0], [4, 2])) === nothing  # true
NetworkCore.separation_caveat(["x"])          # the same sentence, from stored names
NetworkCore.separation_caveat(String[]) === nothing                                                    # true
```
"""
function separation_caveat(v::SeparationVerdict, names=nothing)
    v.separated || return nothing
    return "separation: the coefficient(s) on $(_separated_term_list(v, names)) are " *
           "infinite: the objective has no finite maximum ($(_units_phrase(v)) predicted " *
           "perfectly); the reported values are where the optimizer stopped, and z " *
           "values, p-values and confidence intervals are withheld"
end

function separation_caveat(terms::AbstractVector{<:AbstractString})
    isempty(terms) && return nothing
    return "separation: the coefficient(s) on " * join(("`" * t * "`" for t in terms), ", ") *
           " are infinite: the objective has no finite maximum; the reported values are " *
           "where the optimizer stopped, and z values, p-values and confidence " *
           "intervals are withheld"
end

"""
    z_pvalues(estimates, std_errors, verdict::SeparationVerdict) -> NamedTuple{(:z, :p)}

[`z_pvalues`](@ref) under the separation policy: the usual Wald z values and
two-sided p-values when the verdict is not separated, and `NaN` for every
coefficient when it is (a fit with no finite maximum supports no Wald
inference on any coefficient, since the others are estimated conditionally
on the separated ones being infinite).

# Example
```julia
using NetworkCore
x = collect(-3.0:3.0)
v = NetworkCore.logistic_separation(hcat(ones(7), x), x .> 0)
z_pvalues([1.0, 20.0], [0.5, 900.0], v).p          # [NaN, NaN]
```
"""
function z_pvalues(estimates::AbstractVector{<:Real}, std_errors::AbstractVector{<:Real},
                   v::SeparationVerdict)
    v.separated || return z_pvalues(estimates, std_errors)
    length(estimates) == length(std_errors) || throw(DimensionMismatch(
        "z_pvalues: $(length(estimates)) estimates but $(length(std_errors)) standard errors"))
    return (z=fill(NaN, length(estimates)), p=fill(NaN, length(estimates)))
end
