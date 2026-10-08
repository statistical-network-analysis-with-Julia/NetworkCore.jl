# The ecosystem conversion contract.
#
# The network stack has adapters running in every direction — `Network` ↔
# `DynamicNetwork`, `DynamicNetwork` → snapshots / contact sequences /
# relational-event sequences, `Network` → Siena panel matrices. Every one of
# them maps a rich source object onto a target that can hold *less*, and the
# dangerous failure mode is the implicit one: returning a plausible-looking
# object from which a field quietly vanished. A snapshot of a partially
# observed network that reports zero unobserved dyads is not a snapshot, it is
# a fabrication.
#
# So the ecosystem states one rule for conversions, mirroring the missing-data
# contract in `missing.jl`:
#
#   1. Preserve everything the target can represent. A field that *can* survive
#      the trip must survive it — silently dropping representable metadata is
#      never acceptable.
#   2. What the target cannot represent is either REJECTED (the conversion
#      refuses the input) or DROPPED UNDER AN EXPLICIT POLICY (the caller asked
#      for it in writing). Missing-dyad masks are rejected by default via
#      `require_observed`, because "unobserved" silently becoming "absent" is
#      the exact failure the missing-data contract exists to prevent.
#   3. What is dropped is REPORTED. Every lossy conversion accepts
#      `report=true` and then returns `(result, ::ConversionReport)` naming
#      each dropped field and why.
#
# The per-path invariant table lives in `docs/src/guide/conversion_invariants.md`.

"""
    ConversionReport

Metadata describing what a conversion between network representations dropped.

Every lossy conversion in the ecosystem accepts a `report::Bool=false` keyword;
with `report=true` it returns `(result, ::ConversionReport)` instead of just
`result`. The report names the source and target representations and lists each
dropped field with the reason it could not be carried across.

# Fields
- `source::Symbol` — the source representation (e.g. `:DynamicNetwork`)
- `target::Symbol` — the target representation (e.g. `:Network`)
- `dropped::Vector{Pair{Symbol,String}}` — dropped field ⇒ why

See also [`is_lossless`](@ref), [`dropped_fields`](@ref), [`record_drop!`](@ref).

# Example
```julia
# What an adapter builds (DynamicNetworks' `network_extract(dnet, t; report=true)`
# returns one of these beside the snapshot):
rep = ConversionReport(:DynamicNetwork, :Network)
record_drop!(rep, :spells, "a snapshot has no time axis")
record_drop!(rep, :observation_period, "a snapshot has no observation window")
is_lossless(rep)      # false — a snapshot cannot carry spells
dropped_fields(rep)   # [:spells, :observation_period]
```
"""
struct ConversionReport
    source::Symbol
    target::Symbol
    dropped::Vector{Pair{Symbol, String}}
end

ConversionReport(source::Symbol, target::Symbol) =
    ConversionReport(source, target, Pair{Symbol, String}[])

"""
    record_drop!(report::ConversionReport, field::Symbol, why::AbstractString) -> report

Record that `field` could not be carried across a conversion, and why. Used by
the conversion adapters to build their [`ConversionReport`](@ref); calling it
twice for the same field is not deduplicated.

# Example
```julia
rep = ConversionReport(:DynamicNetwork, :Network)
record_drop!(rep, :spells, "a snapshot has no time axis")
is_lossless(rep)      # false
dropped_fields(rep)   # [:spells]
```
"""
function record_drop!(report::ConversionReport, field::Symbol, why::AbstractString)
    push!(report.dropped, field => String(why))
    return report
end

"""
    is_lossless(report::ConversionReport) -> Bool

Whether the conversion dropped nothing: every field of the source object is
recoverable from the target.

# Example
```julia
rep = ConversionReport(:Network, :Matrix)
is_lossless(rep)                          # true — nothing recorded yet
record_drop!(rep, :vertex_attrs, "a matrix has no vertex attributes")
is_lossless(rep)                          # false
```
"""
is_lossless(report::ConversionReport) = isempty(report.dropped)

"""
    dropped_fields(report::ConversionReport) -> Vector{Symbol}

The names of the fields the conversion dropped, in the order recorded.

# Example
```julia
rep = ConversionReport(:Network, :Matrix)
record_drop!(rep, :vertex_attrs, "a matrix has no vertex attributes")
record_drop!(rep, :network_attrs, "a matrix has no network attributes")
dropped_fields(rep)    # [:vertex_attrs, :network_attrs]
```
"""
dropped_fields(report::ConversionReport) = [f for (f, _) in report.dropped]

function Base.show(io::IO, report::ConversionReport)
    print(io, "ConversionReport($(report.source) → $(report.target): ")
    if is_lossless(report)
        print(io, "lossless)")
    else
        print(io, "$(length(report.dropped)) dropped)")
    end
end

function Base.show(io::IO, ::MIME"text/plain", report::ConversionReport)
    println(io, "ConversionReport: $(report.source) → $(report.target)")
    if is_lossless(report)
        print(io, "  lossless — nothing dropped")
        return
    end
    println(io, "  dropped $(length(report.dropped)) field(s):")
    for (i, (field, why)) in enumerate(report.dropped)
        print(io, "    • ", field, ": ", why)
        i < length(report.dropped) && println(io)
    end
end
