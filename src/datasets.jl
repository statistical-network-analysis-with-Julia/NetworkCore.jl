# Bundled classic teaching datasets.
#
# Data files live in the package's data/ directory as tab-separated text and
# are parsed with a minimal hand-rolled reader (no extra dependencies). See
# data/README.md for provenance details (source URL, upstream package and
# version, file hashes, licence, citation, variable coding); every value was
# extracted programmatically from the upstream R objects by the checked-in
# script data/r/extract_datasets.R (base R + `network` only) and
# cross-checked against the golden fixtures of SNA.jl, ERGM.jl, ERGMEgo.jl
# and Siena.jl.

# Resolved at RUN time, not baked in at precompile time: a `const` built from
# `@__DIR__` froze the build machine's path into the precompile cache, which
# breaks under a relocated depot or a PackageCompiler sysimage. `pkgdir` asks
# the loaded module where it lives now.
_data_dir() = pkgdir(@__MODULE__, "data")

# The names `load_dataset` answers to, in the order the docstring lists them.
const _DATASET_NAMES = (:florentine_marriage, :florentine_business, :sampson,
                        :faux_mesa_high, :s50, :wtc_police_calls)

# Read a tab-separated file with a header line. Returns
# (header::Vector{String}, rows::Vector{Vector{String}}).
function _read_tsv(file::AbstractString)
    lines = readlines(joinpath(_data_dir(), file))
    isempty(lines) && error("empty data file: $file")
    header = split(lines[1], '\t')
    rows = [split(line, '\t') for line in lines[2:end] if !isempty(line)]
    for row in rows
        length(row) == length(header) ||
            error("malformed row in $file: expected $(length(header)) fields")
    end
    return String.(header), [String.(row) for row in rows]
end

# Read a header-less tab-separated integer matrix (one row per line). Every
# row must have the same number of fields; returns a dense `Matrix{Int}`.
function _read_int_matrix_tsv(file::AbstractString)
    lines = filter(!isempty, readlines(joinpath(_data_dir(), file)))
    isempty(lines) && error("empty data file: $file")
    rows = [parse.(Int, split(line, '\t')) for line in lines]
    ncols = length(rows[1])
    for (i, row) in enumerate(rows)
        length(row) == ncols ||
            error("malformed row $i in $file: expected $ncols fields, got $(length(row))")
    end
    M = Matrix{Int}(undef, length(rows), ncols)
    for (i, row) in enumerate(rows)
        M[i, :] = row
    end
    return M
end

function _load_florentine(edge_file::AbstractString)
    net = network(16; directed=false)
    _, vrows = _read_tsv("florentine_vertices.tsv")
    for row in vrows
        v = parse(Int, row[1])
        set_vertex_attribute!(net, :name, v, row[2])
        set_vertex_attribute!(net, :wealth, v, parse(Int, row[3]))
        set_vertex_attribute!(net, :priorates, v, parse(Int, row[4]))
        set_vertex_attribute!(net, :totalties, v, parse(Int, row[5]))
    end
    _, erows = _read_tsv(edge_file)
    for row in erows
        add_edge!(net, parse(Int, row[1]), parse(Int, row[2]))
    end
    return net
end

function _load_sampson()
    net = network(18; directed=true)
    _, vrows = _read_tsv("sampson_vertices.tsv")
    for row in vrows
        v = parse(Int, row[1])
        set_vertex_attribute!(net, :name, v, row[2])
        set_vertex_attribute!(net, :group, v, row[3])
        set_vertex_attribute!(net, :cloisterville, v, parse(Bool, row[4]))
    end
    _, erows = _read_tsv("sampson_edges.tsv")
    for row in erows
        i, j = parse(Int, row[1]), parse(Int, row[2])
        add_edge!(net, i, j)
        set_edge_attribute!(net, :nominations, i, j, parse(Int, row[3]))
    end
    return net
end

# statnet's faux.mesa.high: attribute names are kept EXACTLY as statnet spells
# them (`Grade`, `Race`, `Sex`, capitalised) so that R code and its golden
# fixtures (`nodematch("Grade")`) translate without renaming.
function _load_faux_mesa_high()
    _, vrows = _read_tsv("faux_mesa_high_vertices.tsv")
    net = network(length(vrows); directed=false)
    for row in vrows
        v = parse(Int, row[1])
        set_vertex_attribute!(net, :Grade, v, parse(Int, row[2]))
        set_vertex_attribute!(net, :Race, v, row[3])
        set_vertex_attribute!(net, :Sex, v, row[4])
    end
    _, erows = _read_tsv("faux_mesa_high_edges.tsv")
    for row in erows
        add_edge!(net, parse(Int, row[1]), parse(Int, row[2]))
    end
    return net
end

function _load_s50()
    friendship = Network{Int,true}[]
    for w in 1:3
        A = _read_int_matrix_tsv(joinpath("s50", "friendship_wave$w.tsv"))
        # 0/1 matrices: no edge attribute is stored (store_values=false makes
        # that explicit rather than relying on the all-binary heuristic).
        push!(friendship, network_from_matrix(A; directed=true, store_values=false))
    end
    alcohol = _read_int_matrix_tsv(joinpath("s50", "alcohol.tsv"))
    smoke = _read_int_matrix_tsv(joinpath("s50", "smoke.tsv"))
    return (friendship=friendship, alcohol=alcohol, smoke=smoke)
end

function _load_wtc_police_calls()
    _, erows = _read_tsv("wtc_police_calls_events.tsv")
    events = Matrix{Int}(undef, length(erows), 3)
    for (k, row) in enumerate(erows)
        events[k, 1] = parse(Int, row[1])
        events[k, 2] = parse(Int, row[2])
        events[k, 3] = parse(Int, row[3])
    end
    _, arows = _read_tsv("wtc_police_calls_actors.tsv")
    is_icr = Vector{Bool}(undef, length(arows))
    for row in arows
        is_icr[parse(Int, row[1])] = parse(Int, row[2]) == 1
    end
    return (events=events, n_actors=length(is_icr), is_icr=is_icr)
end

"""
    load_dataset(name::Symbol) -> Network | NamedTuple

Load one of the classic teaching datasets bundled with NetworkCore.jl. Every
call returns a fresh, independent object (nothing is cached), so mutating
what you get back never affects a later load.

Three datasets are single cross-sectional networks and return a `Network`;
three are longitudinal or event data and return a documented `NamedTuple`,
because a panel of waves or a stream of events is not one `Network`.

# Cross-sectional networks (return a `Network`)
- `:florentine_marriage` — Padgett Florentine families marriage network
  (statnet `flomarriage`): 16 families, 20 undirected marriage ties, with
  vertex attributes `:name`, `:wealth` (net wealth in thousands of lira,
  1427), `:priorates` (seats on the civic council, 1282–1344), and
  `:totalties` (total ties in the full 116-family dataset).
- `:florentine_business` — Padgett Florentine families business network
  (statnet `flobusiness`): the same 16 families with the same vertex
  attributes and 15 undirected business ties.
- `:sampson` — Sampson monastery "liking" network (statnet `samplike`):
  18 monks, 88 directed ties (i → j if monk i named monk j among his
  top-three liked monks at any of the three time points), with vertex
  attributes `:name`, `:group` (faction: `"Loyal"`, `"Outcasts"`, or
  `"Turks"`), and `:cloisterville` (`Bool`), and edge attribute
  `:nominations` (number of time points, 1–3, at which the tie was named).
- `:faux_mesa_high` — statnet's `faux.mesa.high` (`ergm` package): a
  simulated in-school friendship network modelled on an Add Health school,
  205 students, 203 undirected friendship ties, `Network{Int,false}`. Vertex
  attributes keep statnet's spelling so R model formulas translate directly:
  `:Grade` (`Int`, 7–12), `:Race` (`String`: `"Black"`, `"Hisp"`, `"NatAm"`,
  `"Other"`, `"White"`) and `:Sex` (`String`: `"F"`, `"M"`). 163 of the 203
  ties join students of the same grade (the `nodematch("Grade")` statistic
  ERGMEgo.jl's golden fixture pins).

# Longitudinal and event data (return a `NamedTuple`)
- `:s50` — the RSiena `s50` excerpt of the Teenage Friends and Lifestyle
  Study (50 girls, three yearly waves). Returns
  `(friendship, alcohol, smoke)`:
  `friendship::Vector{Network{Int,true}}` — three directed 50-vertex
  friendship networks (113, 116 and 122 ties), `alcohol::Matrix{Int}` (50 × 3,
  per-wave alcohol use coded 1 = none … 5 = more than once a week) and
  `smoke::Matrix{Int}` (50 × 3, per-wave smoking coded 1 = none, 2 =
  occasional, 3 = regular). Feeds Siena.jl directly:
  `DependentNetwork(:friendship, s50.friendship)` (with NetworkCore.jl loaded)
  and `DependentBehavior(:alcohol, [s50.alcohol[:, w] for w in 1:3])`.
- `:wtc_police_calls` — the World Trade Center police radio-call event
  stream of Butts, Petrescu-Prahova & Cross (2007), the running example of
  R's `relevent`. Returns `(events, n_actors, is_icr)`:
  `events::Matrix{Int}` (481 × 3: event `number` 1…481, `sender`, `receiver`,
  in observed order — ordinal time, no clock), `n_actors = 37` (35 of whom
  ever call or are called; the actor universe is the unit, not the
  participants) and `is_icr::Vector{Bool}` (length 37; `true` for the three
  institutionalised coordinator roles, actors 13, 14 and 32). Feeds REM.jl
  and Revel.jl:
  `[Event(events[k, 2], events[k, 3], Float64(events[k, 1])) for k in
  1:size(events, 1)]` is the event vector, wrapped as
  `EventSequence(evs; actors=ActorSet(1:37))` for `fit_rem` or passed with
  `n_actors` to `fit_revel`. Use `Float64` times: REM.jl accepts `Int` or
  `Float64`, but Revel.jl's statistics run on a
  `Float64` clock and an `Int`-time event vector fails with a `MethodError`.

# Provenance and licenses
The Florentine data are from Padgett & Ansell (1993), "Robust action and
the rise of the Medici, 1400–1434", *AJS* 98(6); the Sampson data are from
Sampson (1968), *A novitiate in a period of change*, PhD thesis, Cornell
University; `faux.mesa.high` is a simulated network of Hunter, Goodreau &
Handcock (2008) fitted to Add Health data (Resnick et al. 1997); the WTC
calls are from Butts, Petrescu-Prahova & Cross (2007), *J. Math. Sociol.*
31(2); s50 is from West & Sweeting (1996) / Michell & Amos (1997) as
distributed with RSiena. All are redistributed as extracted from the
upstream R objects (statnet `ergm`, GPL-3; statnet's `relevent` workshop,
GPL-2+; RSiena, GPL-2+) by `data/r/extract_datasets.R`; `data/README.md`
records URLs, file hashes, versions and citations. The package code is MIT;
the data keep the provenance stated here, and the Datasets guide tabulates the
upstream licence of each dataset.

# Example
```julia
using NetworkCore

flo = load_dataset(:florentine_marriage)
nv(flo)                                   # 16
ne(flo)                                   # 20
get_vertex_attribute(flo, :wealth, 9)     # 103 (Medici)

fmh = load_dataset(:faux_mesa_high)
count(==(7), vertex_attribute_vector(fmh, :Grade, Int))   # 62 seventh-graders

s50 = load_dataset(:s50)
ne.(s50.friendship)                       # [113, 116, 122]
size(s50.alcohol)                         # (50, 3)

wtc = load_dataset(:wtc_police_calls)
size(wtc.events)                          # (481, 3)
findall(wtc.is_icr)                       # [13, 14, 32]
```
"""
function load_dataset(name::Symbol)
    if name === :florentine_marriage
        return _load_florentine("florentine_marriage_edges.tsv")
    elseif name === :florentine_business
        return _load_florentine("florentine_business_edges.tsv")
    elseif name === :sampson
        return _load_sampson()
    elseif name === :faux_mesa_high
        return _load_faux_mesa_high()
    elseif name === :s50
        return _load_s50()
    elseif name === :wtc_police_calls
        return _load_wtc_police_calls()
    else
        throw(ArgumentError(
            "unknown dataset :$name; available datasets are " *
            join((":$n" for n in _DATASET_NAMES), ", ", " and ")))
    end
end
