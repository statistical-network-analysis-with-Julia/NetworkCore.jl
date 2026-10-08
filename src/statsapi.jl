# The ecosystem's StatsAPI surface, and the ONE checker every model package's
# testset calls.
#
# Every fitted-model type in the ecosystem (ERGMResult, SienaResult, REMResult,
# ...) is expected to answer the StatsAPI verbs wherever the quantity exists:
#
#   coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic, coeftable,
#   coefnames
#
# NetworkCore.jl defines NO method on any of them (it has no fitted-model type and
# never touches a foreign type); it re-exports `coeftable` and `coefnames` so that
# every package extends the same bindings, provides `CoefficientTable`
# (results.jl) as the common return value of `coeftable`, and provides
# `check_statsapi` so that the claim "this fit implements the surface" is a
# testset line, not prose: a migration guide once promised six verbs "on every
# fitted model" while `hasmethod` found them on a third of the types.

# The verbs `check_statsapi` requires by default.
const STATSAPI_VERBS = (:coef, :stderror, :vcov, :confint, :loglikelihood,
                        :nobs, :dof, :aic, :bic, :coeftable)
# Verbs of the surface that are not (yet) required by default: `check_statsapi`
# checks them when they are listed in `required=`, and checks the consistency of
# any that a fit implements even when they are not listed. `coefnames` joined the
# surface after the other ten; it becomes a default requirement once every model
# package defines it.
const STATSAPI_OPTIONAL_VERBS = (:coefnames,)

"""
    check_statsapi(fit; required=(:coef, :stderror, :vcov, :confint,
                                  :loglikelihood, :nobs, :dof, :aic, :bic,
                                  :coeftable),
                   strict::Bool=false) -> NamedTuple

Report which verbs of the ecosystem's **StatsAPI surface** a fitted-model
object implements: for each name in `required`, whether
`hasmethod(StatsAPI.<verb>, Tuple{typeof(fit)})` holds. When the methods exist
their results are also checked for consistency — `length(coef(fit)) ==
length(stderror(fit)) == size(vcov(fit), 1) == size(vcov(fit), 2)`,
`coeftable(fit)` must be a [`CoefficientTable`](@ref) with one row per
coefficient, and `coefnames(fit)` must be a vector of strings, one per
coefficient, equal to `coeftable(fit).names` — and a verb whose method exists
but throws or disagrees is reported `false` too.

The surface has an eleventh verb, [`coefnames`](@ref), which is not in the
default `required` list yet: list it (`required=(..., :coefnames)`) to require
it. A fit that defines `coefnames` has it checked for consistency even when it
is not listed, and with `strict=true` an inconsistent `coefnames` throws.

Returns a `NamedTuple` keyed by `required` with `Bool` values, so a package's
testset can read `@test all(check_statsapi(fit))` or inspect the gaps. With
`strict=true` an incomplete surface throws an `ArgumentError` naming every
missing or inconsistent verb, so the idiom
`@test check_statsapi(fit; strict=true) !== nothing` pins the full surface in
one line and the failure message says exactly what is absent.

`required` may be any subset of the eleven verbs (a model with no likelihood
legitimately omits `loglikelihood`/`aic`/`bic`; pass the verbs it *should*
have). A name that is not a `StatsAPI` function is an `ArgumentError`.

# Example
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
check_statsapi(fit).aic     # false — no method StatsAPI.aic(::TinyFit)
```

See also [`CoefficientTable`](@ref), [`coeftable`](@ref), [`coefnames`](@ref).
"""
function check_statsapi(fit; required=STATSAPI_VERBS, strict::Bool=false)
    verbs = Tuple(Symbol.(collect(required)))
    isempty(verbs) && throw(ArgumentError("check_statsapi: required must name at " *
                                          "least one StatsAPI verb"))
    T = typeof(fit)
    tname = string(nameof(T))
    problems = String[]
    status = Dict{Symbol, Bool}()
    surface = (STATSAPI_VERBS..., STATSAPI_OPTIONAL_VERBS...)
    for f in verbs
        f in surface && isdefined(StatsAPI, f) && getproperty(StatsAPI, f) isa Function ||
            throw(ArgumentError("check_statsapi: :$f is not a verb of the StatsAPI " *
                                "surface $(surface)"))
        ok = hasmethod(getproperty(StatsAPI, f), Tuple{T})
        ok || push!(problems, "no method StatsAPI.$f(::$tname)")
        status[f] = ok
    end
    # An optional verb the fit implements is checked even when not required.
    checked = copy(status)
    for f in STATSAPI_OPTIONAL_VERBS
        haskey(checked, f) || (checked[f] = hasmethod(getproperty(StatsAPI, f), Tuple{T}))
    end

    # Consistency of what exists. A method that throws on its own fit is not
    # an implementation, so it is reported false with the error attached.
    function _call(f)
        try
            return getproperty(StatsAPI, f)(fit)
        catch err
            status[f] = false
            push!(problems, "StatsAPI.$f(::$tname) threw " *
                            sprint(showerror, err; context=:compact => true))
            return nothing
        end
    end
    p = get(status, :coef, false) ? _call(:coef) : nothing
    n = p === nothing ? nothing : length(p)
    if get(status, :stderror, false)
        se = _call(:stderror)
        if se !== nothing && n !== nothing && length(se) != n
            status[:stderror] = false
            push!(problems, "stderror(::$tname) has length $(length(se)) but " *
                            "coef(::$tname) has length $n")
        end
    end
    if get(status, :vcov, false)
        V = _call(:vcov)
        if V !== nothing
            if ndims(V) != 2 || size(V, 1) != size(V, 2)
                status[:vcov] = false
                push!(problems, "vcov(::$tname) is not a square matrix " *
                                "(size $(size(V)))")
            elseif n !== nothing && size(V, 1) != n
                status[:vcov] = false
                push!(problems, "vcov(::$tname) is $(size(V, 1))×$(size(V, 2)) " *
                                "but coef(::$tname) has length $n")
            end
        end
    end
    tbl = nothing
    if get(status, :coeftable, false)
        tbl = _call(:coeftable)
        if tbl !== nothing
            if !(tbl isa CoefficientTable)
                status[:coeftable] = false
                push!(problems, "coeftable(::$tname) returns a " *
                                "$(nameof(typeof(tbl))), not a CoefficientTable")
                tbl = nothing
            elseif n !== nothing && length(tbl) != n
                status[:coeftable] = false
                push!(problems, "coeftable(::$tname) has $(length(tbl)) rows but " *
                                "coef(::$tname) has length $n")
            end
        end
    end
    if checked[:coefnames]
        # `coefnames` needs the coefficients to be checked against: call `coef`
        # here when it was not among the required verbs.
        if n === nothing && !haskey(status, :coef) && hasmethod(StatsAPI.coef, Tuple{T})
            p = try StatsAPI.coef(fit) catch; nothing end
            n = p === nothing ? nothing : length(p)
        end
        if tbl === nothing && !haskey(status, :coeftable) &&
           hasmethod(StatsAPI.coeftable, Tuple{T})
            tbl = try StatsAPI.coeftable(fit) catch; nothing end
            tbl isa CoefficientTable || (tbl = nothing)
        end
        cn = _call(:coefnames)
        bad = nothing
        if cn !== nothing
            if !(cn isa AbstractVector{<:AbstractString})
                bad = "coefnames(::$tname) returns a $(nameof(typeof(cn))), not a " *
                      "vector of strings"
            elseif n !== nothing && length(cn) != n
                bad = "coefnames(::$tname) has length $(length(cn)) but " *
                      "coef(::$tname) has length $n"
            elseif tbl !== nothing && tbl.names != cn
                bad = "coefnames(::$tname) differs from coeftable(::$tname).names"
            end
        end
        if bad !== nothing
            status[:coefnames] = false
            push!(problems, bad)
        end
    end

    result = NamedTuple{verbs}(ntuple(k -> status[verbs[k]], length(verbs)))
    if strict && !isempty(problems)
        missing_verbs = [f for f in (verbs..., STATSAPI_OPTIONAL_VERBS...)
                         if !get(status, f, true)]
        unique!(missing_verbs)
        throw(ArgumentError("check_statsapi: $tname does not implement the " *
                            "StatsAPI surface — missing or inconsistent: " *
                            "$(join(string.(missing_verbs), ", ")). " *
                            join(problems, "; ")))
    end
    return result
end
