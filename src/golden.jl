# Golden-fixture harness: frozen R output, with provenance.
#
# Issue #8 across the ecosystem ("Estimator validation is uneven"): several
# packages carry `# golden value from R` comments with no record of WHICH R,
# which package version, which seed, or what tolerance was thought acceptable.
# A bare number in a comment cannot be regenerated, cannot be audited, and
# silently rots when the reference implementation changes. For long-lived
# academic software that is not evidence — it is folklore.
#
# So a fixture here is a FILE, not a number, and it must carry:
#
#   [provenance]  the R version, the exact package versions, the seed, the
#                 script that produced it, and the date. Enough to regenerate.
#   [tolerance]   what agreement was deemed acceptable, and WHY — Monte-Carlo
#                 estimators cannot be compared to machine precision, and the
#                 tolerance is a scientific claim that must be stated and checkable.
#   [values]      the frozen output itself.
#
# The R script that generates a fixture lives beside it (`r/<name>.R`), in the
# repository, so the whole chain is reproducible by anyone with R.
#
# TOML is stdlib (read-only parse), so this adds no dependency to NetworkCore.jl.

using TOML

"""
    GoldenFixture

Frozen output from a reference implementation (R statnet / RSiena / relevent),
with the provenance needed to regenerate and audit it.

# Fields
- `name` — the fixture's identifier
- `values::Dict{String,Any}` — the frozen numbers (coefficients, SEs, targets…)
- `tolerance::Dict{String,Any}` — the agreed tolerances, per value or global
- `provenance::Dict{String,Any}` — R version, package versions, seed, script, date
- `script_path::String` — the absolute path of the generating script named by
  `provenance["script"]`, resolved on disk by [`load_golden`](@ref) (a fixture
  whose script cannot be found does not load)

Load with [`load_golden`](@ref); compare with [`check_golden`](@ref). The
struct has no external constructor: a `GoldenFixture` only ever comes from a
provenanced file on disk.

# Example
```julia
g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "network_density.toml"))
g isa GoldenFixture             # true
g.provenance["network_version"] # the R `network` version the values came from
isfile(g.script_path)           # true — the generating R script, resolved on disk
```
"""
struct GoldenFixture
    name::String
    values::Dict{String,Any}
    tolerance::Dict{String,Any}
    provenance::Dict{String,Any}
    script_path::String
end

"""
    load_golden(path::AbstractString; root=nothing) -> GoldenFixture

Load a golden fixture from a TOML file and **refuse to load one that cannot be
audited**. A fixture missing its `[provenance]` block — the R version, the
package versions, the seed, and the generating script — throws, because an
unprovenanced golden value is not a check on anything: it cannot be regenerated
when the reference implementation changes, and a future maintainer cannot tell
whether a mismatch means their code broke or R moved.

The named `provenance.script` must also **exist**. It is resolved by walking up
the ancestors of the fixture's directory — the directory itself, its parent,
and so on to the filesystem root — until `joinpath(dir, script)` is a file, so
the conventional `script = "test/fixtures/r/<name>.R"` (relative to the
package root) resolves from any working directory, and an ad-hoc fixture can
keep its script beside it as `r/<name>.R`. Pass `root=` to start the walk from
a different directory (a fixture copied out of its package, say). An absolute
`script` is used as given. Relative paths are normalized to the platform's
path separators. The resolved absolute path is stored in the
returned fixture's `script_path` field and printed by `show` and
[`golden_report`](@ref). A fixture whose script cannot be found throws an
`ArgumentError` naming the fixture, the script string and every directory
searched, for the same reason an unprovenanced one does.

# Example
```julia
using NetworkCore

g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "florentine_sna.toml"))
g.provenance["sna_version"]      # "2.8"
g.values["density"]              # the frozen sna::gden value
isfile(g.script_path)            # true — test/fixtures/r/harness_selftest.R
```
"""
function load_golden(path::AbstractString; root::Union{Nothing,AbstractString}=nothing)
    isfile(path) || throw(ArgumentError("no golden fixture at $path"))
    raw = TOML.parsefile(path)

    prov = get(raw, "provenance", nothing)
    prov isa Dict ||
        throw(ArgumentError(
            "golden fixture $path has no [provenance] block. A frozen reference " *
            "value without provenance cannot be regenerated or audited — when it " *
            "disagrees with the Julia code, nobody can tell whether the code broke " *
            "or the reference implementation moved. Record at least `r_version`, " *
            "the reference package's version, `seed`, and `script`."))

    for key in ("r_version", "seed", "script")
        haskey(prov, key) ||
            throw(ArgumentError("golden fixture $path is missing " *
                                "provenance.$key (required: r_version, seed, script)"))
    end

    # The generating script must EXIST: a provenance block naming a script
    # nobody can find is folklore with extra steps. `script` is recorded
    # relative to the package root (`test/fixtures/r/<name>.R`), so walk up
    # from the fixture's directory (or from `root`) until it is found; the
    # first ancestor tried is the fixture's own directory, which serves ad-hoc
    # fixtures that keep their script beside them.
    script = String(prov["script"])
    start = isnothing(root) ? dirname(abspath(path)) : abspath(root)
    script_path, searched = _resolve_fixture_script(start, script)
    isnothing(script_path) &&
        throw(ArgumentError(
            "golden fixture $path names provenance.script = $(repr(script)), but " *
            "no such file exists in any ancestor of $(start) (searched: " *
            join(searched, ", ") * "). A fixture whose generating script is " *
            "gone cannot be regenerated; check the script in beside the fixture " *
            "(test/fixtures/r/), fix the path, or pass `root=` to name the " *
            "directory it is relative to."))

    values = get(raw, "values", Dict{String,Any}())
    isempty(values) &&
        throw(ArgumentError("golden fixture $path has no [values] block"))

    return GoldenFixture(get(raw, "name", basename(path)),
                         values,
                         get(raw, "tolerance", Dict{String,Any}()),
                         prov,
                         script_path)
end

# Locate `script` (as recorded in a fixture's provenance) on disk: an absolute
# path as given, otherwise `joinpath(dir, script)` for `dir` = `start` and each
# of its ancestors up to the filesystem root. Returns `(path, searched)`:
# the first existing candidate (or `nothing`) and the directories tried.
function _resolve_fixture_script(start::AbstractString, script::AbstractString)
    if isabspath(script)
        return (isfile(script) ? String(script) : nothing), [dirname(script)]
    end
    searched = String[]
    dir = String(abspath(start))
    while true
        push!(searched, dir)
        candidate = normpath(joinpath(dir, script))
        isfile(candidate) && return String(candidate), searched
        parent = dirname(dir)
        parent == dir && return nothing, searched
        dir = parent
    end
end

"""
    golden_tolerance(g::GoldenFixture, key::AbstractString) -> Float64

The tolerance agreed for `key`, falling back to the fixture's `default` and then
to `0.0`. Comparing a Monte-Carlo estimator to a reference at machine precision
is meaningless, so the tolerance is stated in the fixture rather than invented at
the call site — which makes it auditable.

# Example
```julia
g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "network_density.toml"))
golden_tolerance(g, "directed_noloops")    # 1.0e-12 — the fixture's `default`
```
"""
function golden_tolerance(g::GoldenFixture, key::AbstractString)
    haskey(g.tolerance, key) && return Float64(g.tolerance[key])
    haskey(g.tolerance, "default") && return Float64(g.tolerance["default"])
    return 0.0
end

"""
    check_golden(g::GoldenFixture, key::AbstractString, actual) -> Bool

Whether `actual` agrees with the frozen value `g.values[key]` within the
tolerance the fixture declares for it (see [`golden_tolerance`](@ref)).
Element-wise for vectors; lengths must match.

Returns a `Bool` so it drops straight into `@test`. Use
[`golden_report`](@ref) to see *why* a comparison failed.

# Example
```julia
g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "network_density.toml"))
check_golden(g, "directed_noloops", 0.25)           # true
check_golden(g, "directed_noloops", 0.25 + 1e-6)    # false — outside the declared tolerance
```
"""
function check_golden(g::GoldenFixture, key::AbstractString, actual)
    haskey(g.values, key) ||
        throw(ArgumentError("golden fixture $(g.name) has no value named $(repr(key)); " *
                            "it has: $(join(sort(collect(keys(g.values))), ", "))"))
    expected = g.values[key]
    atol = golden_tolerance(g, key)

    if expected isa AbstractVector && actual isa AbstractVector
        length(expected) == length(actual) || return false
        return all(isapprox(Float64(e), Float64(a); atol=atol)
                   for (e, a) in zip(expected, actual))
    end
    return isapprox(Float64(expected), Float64(actual); atol=atol)
end

"""
    golden_report(g::GoldenFixture, key::AbstractString, actual) -> String

A human-readable comparison of `actual` against the frozen value: the reference,
the observed value, the discrepancy, and the tolerance — plus the provenance, so
a maintainer facing a red test knows which R version and seed produced the
number they are being measured against.

# Example
```julia
g = load_golden(joinpath(pkgdir(NetworkCore), "test", "fixtures", "network_density.toml"))
println(golden_report(g, "directed_noloops", 0.3))
# golden fixture network_density [directed_noloops]  (atol = 1.0e-12)
#   reference (R): 0.25
#   observed (jl): 0.3
#   |difference| : 0.04999999999999999
#   provenance   : date=2026-09-09, ..., r_version=4.6.1, ...
#   script       : /.../test/fixtures/r/network_density.R
```
"""
function golden_report(g::GoldenFixture, key::AbstractString, actual)
    expected = g.values[key]
    atol = golden_tolerance(g, key)
    io = IOBuffer()
    println(io, "golden fixture $(g.name) [$key]  (atol = $atol)")
    println(io, "  reference (R): ", expected)
    println(io, "  observed (jl): ", actual)
    if expected isa AbstractVector && actual isa AbstractVector &&
       length(expected) == length(actual)
        println(io, "  |difference| : ",
                [abs(Float64(e) - Float64(a)) for (e, a) in zip(expected, actual)])
    elseif !(expected isa AbstractVector) && !(actual isa AbstractVector)
        println(io, "  |difference| : ", abs(Float64(expected) - Float64(actual)))
    end
    println(io, "  provenance   : ", join(["$k=$v" for (k, v) in sort(collect(g.provenance); by=first)], ", "))
    print(io, "  script       : ", g.script_path)
    return String(take!(io))
end

function Base.show(io::IO, g::GoldenFixture)
    println(io, "GoldenFixture(", g.name, ")")
    println(io, "  values:     ", join(sort(collect(keys(g.values))), ", "))
    println(io, "  provenance: ",
            join(["$k=$v" for (k, v) in sort(collect(g.provenance); by=first)], ", "))
    print(io, "  script:     ", g.script_path)
end
