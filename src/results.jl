# Shared result-presentation infrastructure.
#
# Every statistical package built on NetworkCore.jl (ERGM.jl, TERGM.jl, Siena.jl,
# REM.jl, ...) prints fitted-model results and goodness-of-fit assessments.
# This file factors the presentation layer out of the model packages so that
# all of them render the same way:
#
# - `print_coeftable` — an R-style coefficient table (Estimate / Std.Error /
#   z value / Pr(>|z|) with significance codes) that model packages call from
#   their `Base.show` methods, passing names and vectors.
# - `format_pvalue` / `signif_code` — p-value display with a `"<1e-16"`-style
#   floor (never `0.0`) and the R significance codes.
# - `z_pvalues` — the ONE z → two-sided-normal-p helper (floored at
#   `floatmin`, NaN-aware), so no model package rolls its own `2(1 - Φ(|z|))`
#   (which underflows to a false 0 beyond |z| ≈ 8.3).
# - `GOFStatistic` / `GOFResult` + `Base.show` — a common container for
#   observed-vs-simulated goodness-of-fit comparisons.
# - `gof` — an empty generic function. Model packages add methods for their
#   own result types to this ONE generic, so `gof(fit)` works uniformly and
#   the packages' exports never collide.
# - `mc_pvalue` — the shared two-sided Monte-Carlo p-value, computed with the
#   `(1 + k)/(N + 1)` estimator so it is never exactly zero.
# - `CoefficientTable` — the generic, inspectable coefficient table that every
#   package's `StatsAPI.coeftable(fit)` returns, holding exactly the vectors
#   `print_coeftable` consumes so the inspected table IS the printed one.
#
# Only stdlib (Printf, Statistics) plus SpecialFunctions (`erfc`, for
# `z_pvalues`) and the dependency-free StatsAPI (the `coeftable` generic) is
# used: NetworkCore.jl stays dependency-light. (StatsBase's `CoefTable` would pull
# in StatsBase and its dependency tree, which is why this file rolls its own
# Printf-based table and its own `CoefficientTable` — named so that `using
# NetworkCore, StatsBase` never makes `CoefTable` ambiguous; Distributions would
# pull in a far larger tree for one `ccdf(Normal(), z)` call, and a hand-rolled
# `erfc` is a correctness liability, so the small, stdlib-like SpecialFunctions
# is the one non-stdlib numerical dependency of the foundation.)

"""
    format_pvalue(p::Real; digits::Int=4, floor::Real=1e-16) -> String

Format a p-value for display. Never prints `0.0`: values below `floor`
render as `"<1e-16"` (R's `format.pval` convention), values too small for
`digits` decimal places switch to scientific notation, and everything else
is printed with `digits` decimal places. `NaN` renders as `"NaN"` (e.g. a
coefficient with an undefined standard error).

Throws `ArgumentError` if `p` is outside `[0, 1]` (and not `NaN`), or if
`digits < 1`.

See also [`signif_code`](@ref), [`print_coeftable`](@ref).

# Example
```julia
format_pvalue(0.03511)    # "0.0351"
format_pvalue(3.2e-7)     # "3.2e-07"
format_pvalue(0.0)        # "<1e-16"
format_pvalue(1e-300)     # "<1e-16"
```
"""
function format_pvalue(p::Real; digits::Int=4, floor::Real=1e-16)
    digits >= 1 || throw(ArgumentError("digits must be at least 1 (got $digits)"))
    isnan(p) && return "NaN"
    0 <= p <= 1 || throw(ArgumentError("p-value must be in [0, 1] (got $p)"))
    if p < floor
        return "<" * Printf.format(Printf.Format("%g"), floor)
    elseif p < 10.0^(-digits)
        return Printf.format(Printf.Format("%.1e"), p)
    else
        return Printf.format(Printf.Format("%.$(digits)f"), p)
    end
end

"""
    z_pvalues(z::AbstractVector{<:Real}) -> Vector{Float64}
    z_pvalues(estimates::AbstractVector{<:Real},
              std_errors::AbstractVector{<:Real}) -> NamedTuple{(:z, :p)}

Two-sided normal (Wald) p-values `2 · P(Z ≥ |z|)` for a vector of z
statistics, computed as `erfc(|z| / √2)` so that the tail is accurate far
beyond where the naive `2 · (1 − Φ(|z|))` cancels to a false `0.0` (around
`|z| ≈ 8.3`). The result is **floored at `floatmin(Float64)`** for finite `z`
(`max(p, floatmin)`, so it is also monotone through the subnormal range) — a
finite statistic never has a p-value of exactly zero — while `±Inf` gives `0.0`
and a `NaN` `z` (an undefined standard error) stays `NaN`.

The two-argument form first forms `z = estimate / std_error`, with
`z = NaN` wherever the standard error is not strictly positive (zero, negative
or `NaN` — the coefficient's uncertainty is undefined, so its test is too), and
returns the NamedTuple `(z = ..., p = ...)` ready for [`print_coeftable`](@ref).

This is the ONE z → p helper of the ecosystem: every model package's
`Base.show` and `coeftable` builds its `Pr(>|z|)` column through it, so a
p-value cannot mean two things in two packages (floored in one, underflowing
in another). Throws `ArgumentError` when the two vectors differ in length.

See also [`format_pvalue`](@ref), [`signif_code`](@ref),
[`print_coeftable`](@ref).

# Examples
```julia
z_pvalues([0.0, 1.96, -2.5758])            # [1.0, 0.04999579029644087, 0.009999713929143034]
z_pvalues([50.0])[1] == floatmin(Float64)  # true — floored, never 0.0

zp = z_pvalues([1.2, -0.3, 0.8], [0.4, 0.1, 0.0])
zp.z                                       # [3.0, -3.0, NaN]
isnan(zp.p[3])                             # true — an undefined SE gives an undefined p-value
```
"""
function z_pvalues(z::AbstractVector{<:Real})
    p = Vector{Float64}(undef, length(z))
    @inbounds for k in eachindex(z, p)
        zk = Float64(z[k])
        if isnan(zk)
            p[k] = NaN
        else
            pk = erfc(abs(zk) / sqrt(2.0))
            # A finite statistic never has p exactly 0 (or a subnormal that
            # would print as one): floor at floatmin, which also keeps the map
            # monotone through the underflow region.
            p[k] = isfinite(zk) ? max(pk, floatmin(Float64)) : pk
        end
    end
    return p
end

function z_pvalues(estimates::AbstractVector{<:Real},
                   std_errors::AbstractVector{<:Real})
    length(estimates) == length(std_errors) ||
        throw(ArgumentError("z_pvalues: $(length(estimates)) estimates but " *
                            "$(length(std_errors)) standard errors"))
    z = Vector{Float64}(undef, length(estimates))
    @inbounds for k in eachindex(estimates, std_errors, z)
        se = Float64(std_errors[k])
        z[k] = se > 0 ? Float64(estimates[k]) / se : NaN
    end
    return (z=z, p=z_pvalues(z))
end

"""
    signif_code(p::Real) -> String

R-style significance code for a p-value: `"***"` for `p < 0.001`, `"**"` for
`p < 0.01`, `"*"` for `p < 0.05`, `"."` for `p < 0.1`, and `""` otherwise
(including `NaN`).

The legend, printed by [`print_coeftable`](@ref):

    Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example
```julia
signif_code(0.0004)    # "***"
signif_code(0.03)      # "*"
signif_code(0.2)       # ""
```
"""
function signif_code(p::Real)
    isnan(p) && return ""
    return p < 0.001 ? "***" :
           p < 0.01  ? "**"  :
           p < 0.05  ? "*"   :
           p < 0.1   ? "."   : ""
end

"""
    SIGNIF_LEGEND

The significance-code legend line printed beneath every coefficient table in the
ecosystem, matching R's convention. Paired with [`signif_code`](@ref), which
maps a p-value to one of the codes.

# Example
```julia
SIGNIF_LEGEND    # "Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"
println(SIGNIF_LEGEND)
```
"""
const SIGNIF_LEGEND = "Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"

"""
    print_coeftable(io::IO, names, estimates, std_errors, p_values;
                    z_values=nothing, digits::Int=4,
                    signif_codes::Bool=true, legend::Bool=signif_codes,
                    header=nothing, p_floor::Real=1e-16)

Print an R-style coefficient table with aligned columns and significance
codes. This is the shared presentation layer for every model package built
on NetworkCore.jl — call it from a result type's `Base.show` method, passing the
coefficient names and vectors (or `show` a [`CoefficientTable`](@ref), which
delegates here).

All vectors must have equal length. Column widths adapt to the widest entry,
so long term names and large coefficients stay aligned. p-values are rendered
with [`format_pvalue`](@ref) (so an underflowed p-value prints as `"<1e-16"`,
never `0.0`).

# Arguments
- `names`: coefficient names (strings or symbols), left-aligned first column
- `estimates`, `p_values`: equal-length numeric vectors
- `std_errors`: standard errors, or `nothing` for a table **without** a
  `Std.Error` column (a permutation test such as `SNA.netlm` has a
  test statistic and a p-value but no standard error). With `nothing`,
  `z_values` must be given explicitly (there is nothing to divide by) and the
  table has three numeric columns.

# Keyword Arguments
- `z_values`: test statistics; defaults to `estimates ./ std_errors`
  (required when `std_errors === nothing`)
- `digits::Int=4`: decimal places for the numeric columns
- `signif_codes::Bool=true`: append `***`/`**`/`*`/`.` per row
- `legend::Bool=signif_codes`: print the significance-code legend below the
  table (pass `legend=false` when printing several tables in one report)
- `header`: the column titles — a 4-tuple `("Estimate", "Std.Error",
  "z value", "Pr(>|z|)")` (the default) with standard errors, a 3-tuple
  `("Estimate", "z value", "Pr(>|z|)")` without; override e.g. to label the
  statistic `"t value"` / `"Pr(>|t|)"` for a t-based model. A tuple of the
  wrong length is an `ArgumentError`.
- `p_floor::Real=1e-16`: the display floor forwarded to `format_pvalue`. A
  Monte-Carlo or permutation p-value has a resolution of `1/reps` and is never
  exactly zero in truth, so pass `p_floor=1/reps` and a zero count prints as
  `"<0.001"` (for 1000 replicates), never `0.0`.

A method without `io` prints to `stdout`.

# Examples
```julia
print_coeftable(stdout,
                ["edges", "gwesp.fixed.0.5"],
                [-3.2145, 0.8113],
                [0.1936, 0.2854],
                [1.2e-62, 0.0045])

# A permutation test: no Std.Error column, a t statistic, resolution 1/1000
print_coeftable(stdout, ["intercept", "x"], [0.12, 0.85], nothing, [0.31, 0.0];
                z_values=[1.02, 6.7],
                header=("Estimate", "t value", "Pr(>=|t|)"), p_floor=1/1000)
```
"""
function print_coeftable(io::IO,
                         names::AbstractVector,
                         estimates::AbstractVector{<:Real},
                         std_errors::Union{Nothing, AbstractVector{<:Real}},
                         p_values::AbstractVector{<:Real};
                         z_values::Union{Nothing, AbstractVector{<:Real}}=nothing,
                         digits::Int=4,
                         signif_codes::Bool=true,
                         legend::Bool=signif_codes,
                         header::Union{Nothing, Tuple{Vararg{AbstractString}}}=nothing,
                         p_floor::Real=1e-16)
    n = length(names)
    has_se = std_errors !== nothing
    if z_values === nothing
        has_se || throw(ArgumentError("print_coeftable: z_values= is required when " *
                                      "std_errors is nothing (there is no standard " *
                                      "error to divide by)"))
        z_values = estimates ./ std_errors
    end
    ncol = has_se ? 4 : 3
    # Bound once, never reassigned: the `ntuple` closure below captures it,
    # and a reassigned captured variable is boxed (JET).
    hdr = header === nothing ?
          (has_se ? ("Estimate", "Std.Error", "z value", "Pr(>|z|)") :
                    ("Estimate", "z value", "Pr(>|z|)")) : header
    length(hdr) == ncol ||
        throw(ArgumentError("print_coeftable: header must have $ncol titles for a " *
                            "table $(has_se ? "with" : "without") a Std.Error " *
                            "column (got $(length(hdr)))"))
    checks = has_se ? (("estimates", estimates), ("std_errors", std_errors),
                       ("z_values", z_values), ("p_values", p_values)) :
                      (("estimates", estimates), ("z_values", z_values),
                       ("p_values", p_values))
    for (label, v) in checks
        length(v) == n ||
            throw(ArgumentError("length mismatch: $n names but " *
                                "$(length(v)) $label"))
    end

    fmt = Printf.Format("%.$(digits)f")
    fnum(x) = Printf.format(fmt, x)
    name_strs = string.(names)
    pcol = [format_pvalue(p; digits=digits, floor=p_floor) for p in p_values]
    cols = has_se ? ([fnum(e) for e in estimates],
                     [fnum(s) for s in std_errors],
                     [fnum(z) for z in z_values],
                     pcol) :
                    ([fnum(e) for e in estimates],
                     [fnum(z) for z in z_values],
                     pcol)

    name_w = maximum(length, name_strs; init=0)
    widths = ntuple(k -> max(length(hdr[k]), maximum(length, cols[k]; init=0)), ncol)

    # Header row (the name column is left blank, as in R's printCoefmat)
    print(io, " "^name_w)
    for k in 1:ncol
        print(io, "  ", lpad(hdr[k], widths[k]))
    end
    println(io)

    for i in 1:n
        print(io, rpad(name_strs[i], name_w))
        for k in 1:ncol
            print(io, "  ", lpad(cols[k][i], widths[k]))
        end
        if signif_codes
            code = signif_code(p_values[i])
            isempty(code) || print(io, " ", code)
        end
        println(io)
    end

    if signif_codes && legend
        println(io, "---")
        println(io, SIGNIF_LEGEND)
    end
    return nothing
end

print_coeftable(names::AbstractVector, estimates::AbstractVector{<:Real},
                std_errors::Union{Nothing, AbstractVector{<:Real}},
                p_values::AbstractVector{<:Real}; kwargs...) =
    print_coeftable(stdout, names, estimates, std_errors, p_values; kwargs...)

#==============================================================================#
# CoefficientTable — the generic, inspectable coefficient table
#==============================================================================#

"""
    CoefficientTable(names, estimates, std_errors;
                     z_values=nothing, p_values=nothing, header=nothing,
                     p_floor=1e-16)

The ecosystem's generic coefficient table: the value every model package
returns from `StatsAPI.coeftable(fit)`. It holds exactly the vectors
[`print_coeftable`](@ref) consumes — names, estimates, standard errors (or
`nothing` for a statistic without one, e.g. a permutation test), test
statistics, p-values — plus the column titles and the p-value display floor,
so `show` renders it through `print_coeftable` and prints identically to the
fit's own `show` method.

When `z_values`/`p_values` are omitted they are filled from
[`z_pvalues`](@ref)`(estimates, std_errors)` — `z = NaN` where the standard
error is not strictly positive, p floored at `floatmin`, never `0.0`. With
`std_errors === nothing`, `z_values` is required and `p_values` defaults to
`z_pvalues(z_values)`. `header` defaults to the four R-style titles (three
without standard errors) and must have the matching length. All vectors must
have the same length (`ArgumentError` otherwise).

The name is `CoefficientTable`, not `CoefTable`: `StatsBase.CoefTable` exists,
and `using NetworkCore, StatsBase` must never leave the name ambiguous. This one
carries no StatsBase dependency and is a plain struct any package can build.

# Fields
- `names::Vector{String}`, `estimates::Vector{Float64}`,
  `std_errors::Union{Nothing, Vector{Float64}}`, `z_values::Vector{Float64}`,
  `p_values::Vector{Float64}`
- `header::Tuple{Vararg{String}}`: column titles (length 4, or 3 without SEs)
- `p_floor::Float64`: the floor `show` passes to `format_pvalue`

# Interface
- `length(tbl)`: number of coefficients
- `tbl[i]`: row `i` as a NamedTuple `(name, estimate, std_error, z_value,
  p_value)` (`std_error === nothing` in a table without SEs)
- `tbl["edges"]` / `tbl[:edges]`: the row with that name (`KeyError` if absent)
- iteration yields the rows in order
- `show(io, tbl)`: [`print_coeftable`](@ref) with the stored header and floor

# Examples
```julia
using NetworkCore

tbl = CoefficientTable(["edges", "mutual"], [-2.1, 1.3], [0.2, 0.4])
tbl.p_values == z_pvalues([-2.1, 1.3], [0.2, 0.4]).p   # true
tbl["mutual"].estimate                               # 1.3
length(tbl)                                          # 2
tbl                                                  # prints the R-style table

# A permutation test: no standard errors, a t statistic, resolution 1/1000
CoefficientTable(["intercept", "x"], [0.12, 0.85], nothing;
                 z_values=[1.02, 6.7], p_values=[0.31, 0.0],
                 header=("Estimate", "t value", "Pr(>=|t|)"), p_floor=1/1000)
```

See also [`coeftable`](@ref), [`check_statsapi`](@ref).
"""
struct CoefficientTable
    names::Vector{String}
    estimates::Vector{Float64}
    std_errors::Union{Nothing, Vector{Float64}}
    z_values::Vector{Float64}
    p_values::Vector{Float64}
    header::Tuple{Vararg{String}}
    p_floor::Float64

    function CoefficientTable(names::AbstractVector,
                              estimates::AbstractVector{<:Real},
                              std_errors::Union{Nothing, AbstractVector{<:Real}};
                              z_values::Union{Nothing, AbstractVector{<:Real}}=nothing,
                              p_values::Union{Nothing, AbstractVector{<:Real}}=nothing,
                              header::Union{Nothing, Tuple{Vararg{AbstractString}}}=nothing,
                              p_floor::Real=1e-16)
        n = length(names)
        has_se = std_errors !== nothing
        length(estimates) == n ||
            throw(ArgumentError("length mismatch: $n names but " *
                                "$(length(estimates)) estimates"))
        if has_se
            length(std_errors) == n ||
                throw(ArgumentError("length mismatch: $n names but " *
                                    "$(length(std_errors)) std_errors"))
            if z_values === nothing
                zp = z_pvalues(estimates, std_errors)
                z_values = zp.z
                p_values === nothing && (p_values = zp.p)
            end
        else
            z_values === nothing &&
                throw(ArgumentError("CoefficientTable: z_values= is required when " *
                                    "std_errors is nothing (there is no standard " *
                                    "error to divide by)"))
        end
        length(z_values) == n ||
            throw(ArgumentError("length mismatch: $n names but " *
                                "$(length(z_values)) z_values"))
        p_values === nothing && (p_values = z_pvalues(z_values))
        length(p_values) == n ||
            throw(ArgumentError("length mismatch: $n names but " *
                                "$(length(p_values)) p_values"))
        ncol = has_se ? 4 : 3
        if header === nothing
            header = has_se ? ("Estimate", "Std.Error", "z value", "Pr(>|z|)") :
                              ("Estimate", "z value", "Pr(>|z|)")
        end
        length(header) == ncol ||
            throw(ArgumentError("CoefficientTable: header must have $ncol titles " *
                                "for a table $(has_se ? "with" : "without") a " *
                                "Std.Error column (got $(length(header)))"))
        0 <= p_floor <= 1 ||
            throw(ArgumentError("CoefficientTable: p_floor must be in [0, 1] " *
                                "(got $p_floor)"))
        return new(string.(names),
                   Vector{Float64}(estimates),
                   has_se ? Vector{Float64}(std_errors) : nothing,
                   Vector{Float64}(z_values),
                   Vector{Float64}(p_values),
                   Tuple(String.(collect(header))),
                   Float64(p_floor))
    end
end

Base.length(tbl::CoefficientTable) = length(tbl.names)
Base.firstindex(tbl::CoefficientTable) = 1
Base.lastindex(tbl::CoefficientTable) = length(tbl)

function Base.getindex(tbl::CoefficientTable, i::Integer)
    1 <= i <= length(tbl) || throw(BoundsError(tbl, i))
    return (name=tbl.names[i],
            estimate=tbl.estimates[i],
            std_error=tbl.std_errors === nothing ? nothing : tbl.std_errors[i],
            z_value=tbl.z_values[i],
            p_value=tbl.p_values[i])
end

function Base.getindex(tbl::CoefficientTable, name::Union{Symbol, AbstractString})
    i = findfirst(==(string(name)), tbl.names)
    i === nothing && throw(KeyError(name))
    return tbl[i]
end

Base.iterate(tbl::CoefficientTable, i::Int=1) =
    i > length(tbl) ? nothing : (tbl[i], i + 1)

Base.show(io::IO, tbl::CoefficientTable) =
    print_coeftable(io, tbl.names, tbl.estimates, tbl.std_errors, tbl.p_values;
                    z_values=tbl.z_values, header=tbl.header, p_floor=tbl.p_floor)

"""
    coeftable(fit) -> CoefficientTable

The `StatsAPI.coeftable` generic, re-exported by NetworkCore.jl so that every
model package in the ecosystem extends the ONE binding (`NetworkCore.coeftable ===
StatsAPI.coeftable`) and `using NetworkCore, REM` — or any other pair — leaves
the verb defined. NetworkCore.jl adds **no method**: each package defines
`StatsAPI.coeftable(fit::ItsOwnResult) = CoefficientTable(...)` on its own
result type, feeding it the same vectors its `show` method passes to
[`print_coeftable`](@ref), so the table a user inspects is the table they see
printed.

This is one verb of the ecosystem's **StatsAPI surface**: `coef`, `stderror`,
`vcov`, `confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`
and [`coefnames`](@ref), all defined on the fitted-model type wherever the quantity
exists. [`check_statsapi`](@ref) reports which of them a result type
implements.

# Example
```julia
using NetworkCore
import StatsAPI

struct MyFit
    est::Vector{Float64}
    se::Vector{Float64}
end
StatsAPI.coeftable(f::MyFit) = CoefficientTable(["a", "b"], f.est, f.se)

tbl = coeftable(MyFit([1.2, -0.3], [0.4, 0.1]))
tbl["a"].p_value ≈ z_pvalues([3.0])[1]   # true
```
"""
coeftable

"""
    coefnames(fit) -> Vector{String}

The `StatsAPI.coefnames` generic: the coefficient labels of a fitted model, in
`coef(fit)` order, as R's `names(coef(fit))`. Re-exported by NetworkCore.jl for
the same reason as [`coeftable`](@ref): every model package extends the ONE
binding (`NetworkCore.coefnames === StatsAPI.coefnames`, also the binding
StatsBase and GLM export), so `using NetworkCore, ERGM, REM` leaves the verb
defined. NetworkCore.jl adds **no method**. A package defines
`StatsAPI.coefnames(fit::ItsOwnResult)` returning the same labels its
`coeftable(fit).names` holds; [`check_statsapi`](@ref) checks that they agree.

# Example
```julia
using NetworkCore
import StatsAPI

struct MyFit
    est::Vector{Float64}
    se::Vector{Float64}
end
StatsAPI.coefnames(::MyFit) = ["edges", "mutual"]
StatsAPI.coef(f::MyFit) = f.est
StatsAPI.coeftable(f::MyFit) = CoefficientTable(coefnames(f), f.est, f.se)

fit = MyFit([-1.2, 0.8], [0.3, 0.2])
coefnames(fit)                                  # ["edges", "mutual"]
coefnames(fit) == coeftable(fit).names          # true
check_statsapi(fit; required=(:coef, :coeftable, :coefnames))
```
"""
coefnames

#==============================================================================#
# Goodness of fit
#==============================================================================#

"""
    mc_pvalue(sim::AbstractVector{<:Real}, obs::Real) -> Float64

Two-sided Monte-Carlo p-value of an observed statistic relative to a vector
of simulated values, using the `(1 + k)/(N + 1)` estimator (the observed
value counts as one draw from its null distribution), so the result is never
exactly zero. This is the shared p-value used by goodness-of-fit assessments
across the ecosystem.

**Undefined statistics give `NaN`, never a small p-value.** A `NaN` observed
value (a statistic undefined on the observed network) returns `NaN`.
Simulated values that are `NaN` (undefined on that draw) are dropped and `N` is
the number of defined draws; if none is defined the result is `NaN`. (Comparing
`NaN` counts no draw as extreme, so before 0.2.0 both cases returned the
"significant" `2/(N+1)`.)

Throws `ArgumentError` if `sim` is empty.

# Example
```julia
mc_pvalue([1.0, 2.0, 3.0, 4.0], 5.0)        # 0.4 — 2 * (1 + 0)/(4 + 1)
mc_pvalue([1.0, NaN, 3.0, 4.0], 5.0)        # 0.5 — 2 * (1 + 0)/(3 + 1), the NaN draw dropped
mc_pvalue([1.0, 2.0], NaN)                  # NaN — undefined observed statistic
```
"""
function mc_pvalue(sim::AbstractVector{<:Real}, obs::Real)
    length(sim) > 0 || throw(ArgumentError("mc_pvalue needs at least one simulated value"))
    isnan(obs) && return NaN
    N = count(!isnan, sim)
    N > 0 || return NaN
    # `>=`/`<=` are false for a NaN draw, so the counts already skip them.
    p_ge = (1 + count(>=(obs), sim)) / (N + 1)
    p_le = (1 + count(<=(obs), sim)) / (N + 1)
    return min(1.0, 2.0 * min(p_ge, p_le))
end

"""
    GOFStatistic(name, labels, observed, simulated; p_values=nothing)

Goodness-of-fit comparison for ONE statistic category (e.g. the degree
distribution, the edgewise-shared-partner distribution, a triad census):
observed values versus a matrix of simulated values, one column per level of
the statistic.

# Arguments
- `name`: statistic name shown as the table heading (e.g. `"degree"`)
- `labels`: one label per level (e.g. `["0", "1", "2", ...]`)
- `observed`: observed value per level (length `k`)
- `simulated`: `n_sim × k` matrix of the statistic in simulated networks

# Keyword Arguments
- `p_values`: per-level p-values; when omitted they are computed with
  [`mc_pvalue`](@ref) (two-sided, `(1 + k)/(N + 1)` estimator)

Collect one or more `GOFStatistic`s in a [`GOFResult`](@ref) for formatted
display.

# Example
```julia
sims = [rand(0:5) for _ in 1:100, _ in 1:3]        # 100 simulations × 3 levels
stat = GOFStatistic("degree", ["0", "1", "2"], [2, 3, 1], sims)
stat.p_values                                      # Monte-Carlo p-value per level
stat                                               # prints the observed / envelope / p table
```
"""
struct GOFStatistic
    name::String
    labels::Vector{String}
    observed::Vector{Float64}
    simulated::Matrix{Float64}
    p_values::Vector{Float64}

    function GOFStatistic(name::AbstractString, labels::AbstractVector,
                          observed::AbstractVector{<:Real},
                          simulated::AbstractMatrix{<:Real};
                          p_values::Union{Nothing, AbstractVector{<:Real}}=nothing)
        k = length(observed)
        length(labels) == k ||
            throw(ArgumentError("$(length(labels)) labels for $k observed values"))
        size(simulated, 2) == k ||
            throw(ArgumentError("simulated must have one column per level: " *
                                "got $(size(simulated, 2)) columns for $k levels"))
        size(simulated, 1) > 0 ||
            throw(ArgumentError("simulated must contain at least one simulation"))
        pv = p_values === nothing ?
             [mc_pvalue(view(simulated, :, j), observed[j]) for j in 1:k] :
             convert(Vector{Float64}, p_values)
        length(pv) == k ||
            throw(ArgumentError("$(length(pv)) p-values for $k levels"))
        return new(String(name), string.(labels),
                   convert(Vector{Float64}, observed),
                   convert(Matrix{Float64}, simulated), pv)
    end
end

"""
    GOFResult(statistics; model="", p_overall=nothing)
    GOFResult(statistic::GOFStatistic; kwargs...)

Goodness-of-fit assessment of a fitted model: one [`GOFStatistic`](@ref) per
statistic category, displayed as observed value, simulation envelope
(2.5% / mean / 97.5% quantiles of the simulated distribution), and
Monte-Carlo p-value per level.

This is the common return type for [`gof`](@ref) methods across the model
packages, so every model's goodness of fit prints the same way.

# Keyword Arguments
- `model`: description shown in the heading (e.g. `"ERGM"`, `"SAOM"`)
- `p_overall`: optional overall p-value (e.g. the Monte-Carlo p-value of a
  Mahalanobis distance, as in RSiena's `sienaGOF`)

# Example
```julia
sims = [rand(0:5) for _ in 1:100, _ in 1:3]
stat = GOFStatistic("degree", ["0", "1", "2"], [2, 3, 1], sims)
GOFResult([stat]; model="ERGM")
```
"""
struct GOFResult
    model::String
    statistics::Vector{GOFStatistic}
    p_overall::Union{Nothing, Float64}

    function GOFResult(statistics::AbstractVector{GOFStatistic};
                       model::AbstractString="",
                       p_overall::Union{Nothing, Real}=nothing)
        isempty(statistics) &&
            throw(ArgumentError("GOFResult needs at least one GOFStatistic"))
        return new(String(model), collect(statistics),
                   p_overall === nothing ? nothing : Float64(p_overall))
    end
end

GOFResult(statistic::GOFStatistic; kwargs...) = GOFResult([statistic]; kwargs...)

"""
    n_simulations(result::GOFResult) -> Int

Number of simulated networks the goodness-of-fit assessment is based on
(the number of rows of the first statistic's `simulated` matrix).

# Example
```julia
sims = [rand(0:5) for _ in 1:100, _ in 1:3]
res = GOFResult(GOFStatistic("degree", ["0", "1", "2"], [2, 3, 1], sims))
n_simulations(res)    # 100
```
"""
n_simulations(result::GOFResult) = size(first(result.statistics).simulated, 1)

# One observed-vs-simulated table (rows = levels of the statistic).
function _print_gof_table(io::IO, stat::GOFStatistic)
    println(io, "Goodness-of-fit for ", stat.name)
    println(io)

    fnum(x) = Printf.format(Printf.Format("%.2f"), x)
    fobs(x) = Printf.format(Printf.Format("%g"), x)
    header = ("Obs", "2.5%", "Mean", "97.5%", "MC p-value")
    k = length(stat.labels)
    # The envelope is taken over the draws on which the statistic is defined
    # (the same draws `mc_pvalue` counts); `quantile` refuses NaN outright.
    defined(j) = filter(!isnan, view(stat.simulated, :, j))
    envelope(j, f) = (d = defined(j); isempty(d) ? "NaN" : fnum(f(d)))
    cols = (map(fobs, stat.observed),
            [envelope(j, d -> quantile(d, 0.025)) for j in 1:k],
            [envelope(j, mean) for j in 1:k],
            [envelope(j, d -> quantile(d, 0.975)) for j in 1:k],
            [format_pvalue(p) for p in stat.p_values])

    label_w = maximum(length, stat.labels; init=0)
    widths = ntuple(c -> max(length(header[c]), maximum(length, cols[c]; init=0)), 5)

    print(io, " "^label_w)
    for c in 1:5
        print(io, "  ", lpad(header[c], widths[c]))
    end
    println(io)
    for i in 1:k
        print(io, rpad(stat.labels[i], label_w))
        for c in 1:5
            print(io, "  ", lpad(cols[c][i], widths[c]))
        end
        println(io)
    end
    return nothing
end

function Base.show(io::IO, result::GOFResult)
    heading = isempty(result.model) ? "Goodness-of-fit assessment" :
              "Goodness-of-fit assessment: " * result.model
    println(io, heading)
    println(io, "=" ^ length(heading))
    println(io, "Based on ", n_simulations(result), " simulated networks")
    if result.p_overall !== nothing
        println(io, "Overall p-value: ", format_pvalue(result.p_overall))
    end
    for stat in result.statistics
        println(io)
        _print_gof_table(io, stat)
    end
    return nothing
end

Base.show(io::IO, stat::GOFStatistic) = _print_gof_table(io, stat)

"""
    gof(result; kwargs...) -> GOFResult

Goodness-of-fit assessment of a fitted statistical network model.

This is an empty generic function: NetworkCore.jl defines no methods. Each model
package (ERGM.jl, TERGM.jl, Siena.jl, REM.jl, ...) adds methods for its own
result types to this ONE shared generic, so `gof(fit)` works uniformly across
the ecosystem and loading several model packages together never produces
colliding `gof` exports. Methods should return a [`GOFResult`](@ref).

# Example
```julia
using Random
import NetworkCore: gof                      # extend the ONE generic, never define a new `gof`
struct ToyFit end
function gof(::ToyFit; n_sim=100, rng=Random.default_rng())
    sims = [rand(rng, 0:5) for _ in 1:n_sim, _ in 1:3]
    GOFResult(GOFStatistic("degree", ["0", "1", "2"], [2, 3, 1], sims); model="toy")
end
n_simulations(gof(ToyFit(); n_sim=50))    # 50
```
"""
function gof end

# ============================================================================
# The shared result-metadata protocol
# ============================================================================
#
# Six issues across the ecosystem (TERGM#1, ERGMCount#2, ERGMRank#1,
# ERGMMulti#1, ERGMEgo#1, REM#2) independently ask the same question in
# different words: *what did this fit actually do?* Is the number an exact
# likelihood or a pseudo-likelihood? Are the standard errors an inverse Hessian
# or a bootstrap? Were unobserved dyads dropped, or conditioned on at face
# value? Prose caveats in a `show` method cannot be inspected, compared, or
# tabulated — and they are exactly what a reader skips.
#
# So the answer is a PROTOCOL, not a struct field. Each model package adds
# methods for its own result type; nothing needs a new field, no constructor
# breaks, and packages that have not opted in are honestly reported as
# `:unspecified` rather than silently assumed to be exact.
#
# `fit_metadata(fit)` collects the seven accessors into one `ResultMetadata`,
# which is what the ecosystem capability matrix is generated from.

"""
    estimand(fit) -> Symbol

The population quantity the fit targets (e.g. `:ergm`, `:stergm_formation`,
`:saom`, `:relational_event`). Defaults to `:unspecified`.

Part of the shared result-metadata protocol; see [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: estimand
struct MyFit end
estimand(::MyFit) = :ergm
estimand(MyFit())          # :ergm
estimand("anything else")  # :unspecified — the default for an undeclared type
```
"""
estimand(::Any) = :unspecified

"""
    objective(fit) -> Symbol

The objective function that was actually maximised — **not** the one the model
is named after. Conventional values:

- `:likelihood` — the exact likelihood of the model
- `:pseudolikelihood` — a dyadwise (or swapwise) pseudo-likelihood
- `:conditional_pseudolikelihood` — CMPLE, as in panel TERGM
- `:mc_likelihood` — a Monte-Carlo approximation to the likelihood (MCMLE)
- `:moment` — method of moments / Robbins-Monro (SAOM)
- `:partial_likelihood` — conditional-logit / case-control partial likelihood

Defaults to `:unspecified`. See [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: objective
struct MPLEFit end
objective(::MPLEFit) = :pseudolikelihood
objective(MPLEFit())    # :pseudolikelihood
```
"""
objective(::Any) = :unspecified

"""
    is_exact(fit) -> Bool

Whether [`objective`](@ref) coincides with the exact likelihood **for this
particular fitted model**, not in general.

This is deliberately a property of the fit rather than of the estimator: a
CMPLE fit of a dyad-independent formula *is* the CMLE, while the same estimator
on a dyad-dependent formula is an approximation with anticonservative standard
errors. The distinction is precisely what a user needs and cannot see.

Defaults to `false` — the conservative answer for a package that has not
declared. See [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: is_exact
struct MPLEFit; dyad_independent::Bool; end
is_exact(fit::MPLEFit) = fit.dyad_independent    # MPLE is the MLE only for a dyad-independent model
is_exact(MPLEFit(true)), is_exact(MPLEFit(false))    # (true, false)
```
"""
is_exact(::Any) = false

"""
    se_method(fit) -> Symbol

How the reported standard errors were obtained. Conventional values:

- `:hessian` — inverse negative Hessian of the objective (for a
  pseudo-likelihood this is generally **anticonservative** under dependence)
- `:fisher` — Fisher information from an MC sample
- `:bootstrap` — parametric or block bootstrap
- `:sandwich` — Godambe / robust sandwich
- `:none` — no standard errors reported

Defaults to `:unspecified`. See [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: se_method
struct BootFit end
se_method(::BootFit) = :bootstrap
se_method(BootFit())    # :bootstrap
```
"""
se_method(::Any) = :unspecified

"""
    missing_method(fit) -> Symbol

How unobserved (masked) dyads were treated. Conventional values:

- `:none` — the data had no masked dyads
- `:available_case` — masked dyads excluded from the objective (e.g. ERGM MPLE)
- `:condition_on_face` — masked dyads frozen at their stored face value, an
  explicitly different estimand from missing-data maximum likelihood
- `:rejected` — the estimator refuses masked data

Defaults to `:unspecified`. See the ecosystem missing-data contract in
`missing.jl` ([`supports_missing`](@ref), [`require_observed`](@ref)).

# Example
```julia
import NetworkCore: missing_method
struct MPLEFit end
missing_method(::MPLEFit) = :available_case    # masked dyads left out of the design
missing_method(MPLEFit())    # :available_case
```
"""
missing_method(::Any) = :unspecified

"""
    tie_method(fit) -> Symbol

How tied event times were handled, for the relational-event models. `:none`
when ties cannot arise or none were present; `:not_applicable` for models
without event times. Otherwise one of [`TIE_POLICIES`](@ref) — the policy the
fit **actually ran**, which is why `:error` can never appear here: a fit that
would have had to break a tie under `ties=:error` threw instead of returning.

Defaults to `:not_applicable`. See [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: tie_method
struct REMFit; ties::Symbol; end
tie_method(fit::REMFit) = fit.ties
tie_method(REMFit(:efron))    # :efron
tie_method("an ERGM fit")     # :not_applicable — no event clock
```
"""
tie_method(::Any) = :not_applicable

"""
    TIE_POLICIES

The ecosystem vocabulary for **tied event times** — one keyword (`ties=`), one
set of symbols, one meaning per symbol, across every model with an event clock
(`REM.fit_rem`, `Revel.fit_revel`).

A tie is not a nuisance to be sorted away: an exact-time or exact-order
likelihood evaluated on tied data is evaluated on information the data does not
contain. So the vocabulary starts from refusal and every alternative is an
explicit, recorded claim about what a simultaneous pair *means*:

- `:error` — refuse to fit when two events share a timestamp. **The default
  wherever the model makes an exact-time or exact-order claim**, which is
  everywhere in this ecosystem. A user with tied data is told, not handed a
  number that depends on an arbitrary sort.
- `:ordered` — break ties in sequence order and apply no correction. The
  legacy behaviour; the tied events are treated as if one had genuinely
  preceded the other and had been able to influence it. Invents information.
- `:breslow` — the Breslow tie correction of the Cox partial likelihood: the
  tied events share one risk set (statistics frozen at the pre-tie state) and
  each contributes the same denominator. Defined only for partial-likelihood /
  ordinal models.
- `:efron` — the Efron tie correction: as Breslow, but the tied cases enter the
  denominator of the *j*-th term with weight `1 − (j−1)/d`, averaging over the
  `d!` orderings instead of ignoring them. Strictly the better approximation and
  what `survival::coxph` defaults to. Defined only for partial-likelihood /
  ordinal models.
- `:batch` — the tied events are one simultaneous batch: none of them could
  have influenced any other (statistics frozen at the pre-tie state) and the
  batch consumes one exposure interval. The coarsened-observation reading, and
  the only one available to an exact-time (hazard) likelihood, where Breslow and
  Efron are undefined because there is no partial likelihood to correct.

Not every symbol makes sense for every model, and one that does not must be
**refused**, never silently accepted — see [`check_tie_policy`](@ref).

# Example
```julia
TIE_POLICIES                       # (:error, :ordered, :breslow, :efron, :batch)
:efron in TIE_POLICIES             # true
check_tie_policy(:efron, (:error, :ordered, :breslow, :efron); model="fit_rem")   # :efron
```
"""
const TIE_POLICIES = (:error, :ordered, :breslow, :efron, :batch)

"""
    check_tie_policy(ties::Symbol, supported; model, reasons=Dict()) -> Symbol

Validate a `ties=` keyword against the [`TIE_POLICIES`](@ref) vocabulary and
against the subset a particular model can actually honour. Returns `ties` when
it is supported; otherwise throws an `ArgumentError` that names `model`, says
**why** the policy does not apply to it (from `reasons`, keyed by symbol) and
lists what is supported instead.

This is the guard that keeps the ecosystem's rule: an unimplemented option
fails loudly, never no-ops.

```julia
ties = :efron
check_tie_policy(ties, (:error, :ordered, :breslow, :efron);
                 model = "fit_rem (conditional-logit partial likelihood)",
                 reasons = Dict(:batch => "a batch of simultaneous events with " *
                                          "the risk set held fixed IS the Breslow " *
                                          "correction here; pass ties=:breslow"))
```
"""
function check_tie_policy(ties::Symbol, supported::Tuple{Vararg{Symbol}};
                          model::AbstractString,
                          reasons::AbstractDict{Symbol,<:AbstractString}=
                              Dict{Symbol,String}())
    ties in supported && return ties
    sup = join(("`:$s`" for s in supported), ", ")
    if ties in TIE_POLICIES
        why = get(reasons, ties, "")
        throw(ArgumentError(
            "tie policy `:$ties` is not defined for $model" *
            (isempty(why) ? "" : ": $why") *
            ". Supported here: $sup."))
    end
    voc = join(("`:$s`" for s in TIE_POLICIES), ", ")
    throw(ArgumentError(
        "unknown tie policy `:$ties`. The ecosystem vocabulary is $voc " *
        "(`NetworkCore.TIE_POLICIES`); $model supports $sup."))
end

"""
    ResultMetadata

What a fitted model actually did, in a form that can be inspected and tabulated
rather than read out of a warning. Built by [`fit_metadata`](@ref) from the seven
accessors of the shared result-metadata protocol.

The `approximations` field carries any free-text caveats the package wants to
attach (e.g. "support truncated at 20; boundary mass 3e-6").

# Example
```julia
md = fit_metadata("a value whose type declares nothing")
md isa ResultMetadata           # true
md.objective, md.is_exact       # (:unspecified, false) — the conservative defaults
md                              # prints one line per field
```
"""
struct ResultMetadata
    estimand::Symbol
    objective::Symbol
    is_exact::Bool
    se_method::Symbol
    missing_method::Symbol
    tie_method::Symbol
    approximations::Vector{String}
end

"""
    approximations(fit) -> Vector{String}

Free-text caveats attached to a fit — the approximations a reader would need to
know about and that do not fit the fixed vocabulary of the other accessors.
Defaults to empty. See [`fit_metadata`](@ref).

# Example
```julia
import NetworkCore: approximations
struct TruncFit end
approximations(::TruncFit) = ["support truncated at 20; boundary mass 3e-6"]
approximations(TruncFit())    # 1-element Vector{String}
approximations(1.0)           # String[] — the default
```
"""
approximations(::Any) = String[]

"""
    fit_metadata(fit) -> ResultMetadata

Collect the shared result-metadata protocol for a fitted model: what it
estimates, what objective it maximised, whether that objective is exact for
*this* model, how the standard errors were computed, how missing dyads and tied
events were treated, and any further caveats.

This is what makes an approximation programmatically inspectable instead of a
sentence in a `show` method, and it is the source the ecosystem capability
matrix is generated from.

```julia
import NetworkCore: estimand, objective, is_exact, se_method, missing_method
struct CMPLEFit end                       # a package's result type
estimand(::CMPLEFit) = :stergm_formation
objective(::CMPLEFit) = :conditional_pseudolikelihood
is_exact(::CMPLEFit) = false
se_method(::CMPLEFit) = :hessian
missing_method(::CMPLEFit) = :available_case

md = fit_metadata(CMPLEFit())
md.is_exact          # false — CMPLE of a dyad-dependent formula
md.se_method         # :hessian, therefore anticonservative here
md.missing_method    # :available_case
```

A package opts in by adding methods for [`estimand`](@ref), [`objective`](@ref),
[`is_exact`](@ref), [`se_method`](@ref), [`missing_method`](@ref),
[`tie_method`](@ref) and [`approximations`](@ref) on its result type. Anything
undeclared reports `:unspecified` (and `is_exact` reports `false`) rather than
being optimistically assumed.
"""
fit_metadata(fit) = ResultMetadata(estimand(fit), objective(fit), is_exact(fit),
                                   se_method(fit), missing_method(fit),
                                   tie_method(fit), approximations(fit))

function Base.show(io::IO, md::ResultMetadata)
    println(io, "ResultMetadata")
    println(io, "  estimand:       ", md.estimand)
    println(io, "  objective:      ", md.objective,
                md.is_exact ? " (exact for this model)" : " (approximate)")
    println(io, "  standard errors: ", md.se_method)
    md.missing_method === :none ||
        println(io, "  missing dyads:  ", md.missing_method)
    md.tie_method === :not_applicable ||
        println(io, "  tied events:    ", md.tie_method)
    isempty(md.approximations) && return
    println(io, "  approximations:")
    for a in md.approximations
        println(io, "    • ", a)
    end
end
