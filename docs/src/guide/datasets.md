# Datasets

NetworkCore.jl bundles a few classic teaching datasets, loaded with
[`load_dataset`](@ref). They are the same networks used throughout the
social network analysis literature and in the R statnet documentation, so
examples translate directly.

```julia
using NetworkCore

flo  = load_dataset(:florentine_marriage)   # Network
flob = load_dataset(:florentine_business)   # Network
samp = load_dataset(:sampson)               # Network
fmh  = load_dataset(:faux_mesa_high)        # Network
s50  = load_dataset(:s50)                   # NamedTuple: (friendship, alcohol, smoke)
wtc  = load_dataset(:wtc_police_calls)      # NamedTuple: (events, n_actors, is_icr)
```

Three of the datasets are single cross-sectional networks and return a
`Network`; the other three are a longitudinal panel or an event stream and
return a documented `NamedTuple`, because a set of waves or a sequence of
events is not one `Network`. Every call returns a fresh object.

## Available datasets

### `:florentine_marriage`

Padgett's Florentine families marriage network (statnet's `flomarriage`):
16 Renaissance Florentine families and 20 undirected marriage ties. The
Pucci family (vertex 12) is an isolate. Vertex attributes:

- `:name` — family name (`"Medici"`, `"Strozzi"`, ...)
- `:wealth` — net wealth in thousands of lira (1427)
- `:priorates` — seats on the civic council, 1282–1344
- `:totalties` — total ties in the full 116-family Padgett dataset

```julia
flo = load_dataset(:florentine_marriage)
nv(flo), ne(flo)                        # (16, 20)
get_vertex_attribute(flo, :wealth, 9)   # 103 (Medici)
wealth = vertex_attribute_vector(flo, :wealth, Float64)
```

### `:florentine_business`

The business (loan/partnership) ties among the same 16 families (statnet's
`flobusiness`): 15 undirected edges, same vertex attributes as
`:florentine_marriage`.

### `:sampson`

Sampson's monastery "liking" network (statnet's `samplike`): 18 trainee
monks, 88 directed ties. A tie i → j means monk i named monk j among his
top-three liked monks at any of the three time points of the study. Vertex
attributes:

- `:name` — monk's name (`"John Bosco"`, ...)
- `:group` — faction: `"Loyal"` (Loyal Opposition), `"Outcasts"`, or
  `"Turks"` (Young Turks)
- `:cloisterville` — `Bool`, attended the "Cloisterville" minor seminary

Edge attribute `:nominations` gives the number of time points (1–3) at
which each tie was named.

```julia
samp = load_dataset(:sampson)
groups = vertex_attribute_vector(samp, :group, String)
noms   = get_edge_attribute(samp, :nominations, Int)
```

### `:faux_mesa_high`

statnet's `faux.mesa.high` (from the `ergm` package): a *simulated*
in-school friendship network fitted to an Add Health school, 205 students
and 203 undirected friendship ties. It is the running example of the ergm
and ergm.ego tutorials, so R formulas such as
`faux.mesa.high ~ edges + nodematch("Grade")` translate directly. Vertex
attributes keep statnet's capitalised names:

- `:Grade` — `Int`, 7–12
- `:Race` — `String`: `"Black"`, `"Hisp"`, `"NatAm"`, `"Other"`, `"White"`
- `:Sex` — `String`: `"F"`, `"M"`

```julia
using Graphs: src, dst                  # edge endpoints (Graphs.jl accessors)
fmh = load_dataset(:faux_mesa_high)
nv(fmh), ne(fmh)                        # (205, 203)
grade = vertex_attribute_vector(fmh, :Grade, Int)
count(==(7), grade)                     # 62 seventh-graders
# the nodematch("Grade") statistic: 163 of the 203 ties are within-grade
count(grade[src(e)] == grade[dst(e)] for e in edges(fmh))
```

Being simulated, it carries no real student's data; the counts above (and
`[203, 163]` as the edges/nodematch targets) are what ERGMEgo.jl's
provenanced `fauxmesa_ego_census` fixture pins against R's `ergm.ego`.

### `:s50`

The RSiena `s50` excerpt of the Teenage Friends and Lifestyle Study
(West & Sweeting 1996; Michell & Amos 1997): 50 girls in a Glasgow school,
three yearly waves of a directed friendship network plus two behaviour
variables. Returns a `NamedTuple`:

- `friendship::Vector{Network{Int,true}}` — three directed 50-vertex
  networks (RSiena `s501`, `s502`, `s503`) with 113, 116 and 122 ties
- `alcohol::Matrix{Int}` — 50 × 3 (RSiena `s50a`), per-wave alcohol use:
  1 = none, 2 = once or twice a year, 3 = once a month, 4 = once a week,
  5 = more than once a week
- `smoke::Matrix{Int}` — 50 × 3 (RSiena `s50s`), per-wave smoking:
  1 = none, 2 = occasional, 3 = regular

```julia
s50 = load_dataset(:s50)
ne.(s50.friendship)                     # [113, 116, 122]
size(s50.alcohol), size(s50.smoke)      # ((50, 3), (50, 3))
# wave-to-wave stability: ties present in both wave 1 and wave 2
A1, A2 = as_matrix(s50.friendship[1]), as_matrix(s50.friendship[2])
sum(A1 .* A2)
```

This is the data every Siena.jl example uses. With `using Siena` the waves
go straight into a dependent variable and the behaviours into a behaviour
variable or covariate:

<!-- skip-check -->
```julia
using NetworkCore, Siena
s50 = load_dataset(:s50)
data = siena_data()
add_nodeset!(data, NodeSet(50))
add_dependent!(data, DependentNetwork(:friendship, s50.friendship))
add_dependent!(data, DependentBehavior(:alcohol, [s50.alcohol[:, w] for w in 1:3]))
add_covariate!(data, ConstantCovariate(:smoke1, s50.smoke[:, 1]))
```

### `:wtc_police_calls`

The World Trade Center police radio-communication event stream of Butts,
Petrescu-Prahova & Cross (2007), the running example of R's `relevent`
package (`WTCPoliceCalls`, `WTCPoliceIsICR`): 481 radio calls among the 37
named communicants of one police unit on the morning of 11 September 2001.
The data are **ordinal** — the calls are known in order, not by clock time.
Returns a `NamedTuple`:

- `events::Matrix{Int}` — 481 × 3: event `number` (1…481), `sender`,
  `receiver`, in observed order; no self-calls
- `n_actors` — 37 (the actor universe; 35 of them ever call or are called,
  and a relational-event likelihood is conditional on the *universe*, so the
  two silent non-participants matter)
- `is_icr::Vector{Bool}` — length 37, `true` for the three actors (13, 14
  and 32) in *institutionalised coordinator roles*, the covariate the
  original analysis is built around

```julia
wtc = load_dataset(:wtc_police_calls)
size(wtc.events)                        # (481, 3)
findall(wtc.is_icr)                     # [13, 14, 32]
# who calls most?
senders = wtc.events[:, 2]
argmax([count(==(a), senders) for a in 1:wtc.n_actors])   # actor 32, an ICR
```

With `using REM`, or Revel.jl (relational-event effects, goodness of fit
and interval timing, using REM's `Event` type),
the matrix becomes an event sequence — the event number serves as the ordinal
time; Revel's README fits it with `fit_revel(calls, stats, wtc.n_actors)`.
Build the events with `Float64` times: REM.jl accepts either `Int` or
`Float64` event times, but Revel.jl's statistics are defined on a
`Float64` clock, so an `Event` vector with `Int` times fits `fit_rem` and
fails in Revel.jl with a `MethodError`:

<!-- skip-check -->
```julia
using REM
wtc = load_dataset(:wtc_police_calls)
# Float64 time: accepted by REM.jl, required by Revel.jl
evs = [Event(wtc.events[k, 2], wtc.events[k, 3], Float64(wtc.events[k, 1]))
       for k in 1:size(wtc.events, 1)]
seq = EventSequence(evs; actors=ActorSet(1:wtc.n_actors))   # REM.fit_rem
# Revel.fit_revel(evs, stats, wtc.n_actors)   — see Revel.jl's README
```

## Provenance

The Florentine data are from Padgett & Ansell (1993), "Robust action and
the rise of the Medici, 1400–1434", *American Journal of Sociology* 98(6);
the Sampson data are from Sampson (1968), *A novitiate in a period of
change*, PhD thesis, Cornell University. Both are redistributed as parsed
from the statnet `ergm` R package's data files (GPL-3). `faux.mesa.high`
(Hunter, Goodreau & Handcock 2008, simulated from Add Health data, Resnick
et al. 1997) is extracted from the same package's `data/faux.mesa.high.RData`;
the WTC calls (Butts, Petrescu-Prahova & Cross 2007) from statnet's
`relevent` workshop data (the CRAN tarball ships none; GPL-2+); s50 (West &
Sweeting 1996; Michell & Amos 1997) is RSiena's bundled copy (GPL-2+), the
matrices Siena.jl's provenanced `s50_siena07` fixture was generated from
with RSiena 1.6.6. The three are produced by the checked-in script
`data/r/extract_datasets.R` (base R + `network` only), and
`data/README.md` in the package repository records the source URL and
SHA-256 of every upstream file, the written-file hashes, licences, citations
and the variable codings.

## Licences of the bundled data

NetworkCore.jl's code is MIT-licensed. The data files under `data/` were
extracted from R packages distributed under the GNU GPL and are redistributed
here with that provenance stated; whether the extracted values carry the
upstream licence is a question this package does not settle, and anyone
redistributing them should consult the upstream terms.

| Dataset | Original source (citation) | Taken from | Upstream licence |
|---|---|---|---|
| `:florentine_marriage`, `:florentine_business` | Padgett & Ansell (1993); Padgett (1994) | statnet `ergm`, `data/florentine.RData` | GPL-3 |
| `:sampson` | Sampson (1968) | statnet `ergm`, `data/sampson.RData` | GPL-3 |
| `:faux_mesa_high` | Hunter, Goodreau & Handcock (2008), simulated from Add Health (Resnick et al. 1997) | statnet `ergm`, `data/faux.mesa.high.RData` | GPL-3 |
| `:wtc_police_calls` | Butts, Petrescu-Prahova & Cross (2007) | statnet `relevent` workshop, `relevent_workshop.Rdata` | GPL-2+ |
| `:s50` | West & Sweeting (1996); Michell & Amos (1997) | RSiena (`s501`–`s503`, `s50a`, `s50s`) | GPL-2+ |

## API

```@docs
load_dataset
```
