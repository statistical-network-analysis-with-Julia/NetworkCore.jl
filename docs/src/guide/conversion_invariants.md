# Conversion Invariants

The ecosystem has adapters running in every direction — `Network` ↔
`DynamicNetwork`, `DynamicNetwork` → snapshots, contact sequences and
relational-event sequences, `Network` → Siena panel matrices. Every one of
them maps a rich source object onto a target that can hold *less*.

The dangerous failure mode is the implicit one: returning a plausible-looking
object from which a field has quietly vanished. A snapshot of a partially
observed network that reports **zero** unobserved dyads is not a snapshot — it
is a fabrication, and every number computed from it is wrong in a way nothing
downstream can detect.

So the ecosystem states one rule for conversions, mirroring the
[missing-data contract](missing_data.md):

1. **Preserve what the target can represent.** A field that *can* survive the
   trip must survive it. Silently dropping representable metadata is never
   acceptable.
2. **Reject or policy-gate what it cannot.** A conversion whose target has no
   slot for a field either refuses the input, or drops the field under an
   explicit keyword the caller had to type. Missing-dyad masks are **rejected
   by default** (`missing=:error`), because "unobserved" silently becoming
   "absent" is the exact failure the missing-data contract exists to prevent.
   `missing=:face` is the auditable opt-in.
3. **Report what was dropped.** Every lossy conversion accepts `report=true`
   and then returns `(result, ::ConversionReport)` naming each dropped field
   and why.

```julia
using NetworkCore, DynamicNetworks

dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate_vertices!(dnet, [1, 2, 3, 4], 0.0, 10.0)
activate!(dnet, 0.0, 8.0; edge=(1, 2))

snap, rep = network_extract(dnet, 5.0; report=true)
is_lossless(rep)      # false — a snapshot has no time axis
dropped_fields(rep)   # [:spells, :observation_period]
rep                   # pretty-prints each dropped field with its reason
```

## The invariant table

Legend: **P** preserved · **T** transformed (representable, but re-encoded) ·
**D** dropped, and named in the `ConversionReport` · **R** rejected unless an
explicit policy is passed · **—** not applicable to this path.

| Field | `Network`→`DynamicNetwork` [^1] | `DynamicNetwork`→`Network` [^2] | `DynamicNetwork`→`ContactSequence` [^3] | `DynamicNetwork`→`EventSequence` [^4] | `Network`→Siena matrices [^5] | `Network`↔`Matrix` [^6] | `Network`→`EgoData` [^7] | `Network`s→`MultilayerNetwork` [^8] | `Network`→`RankNetwork` [^9] |
|---|---|---|---|---|---|---|---|---|---|
| Directedness | P | P | P | T (undirected → `(min,max)` sender/receiver) | P | T (keyword on the way back) | R (undirected only; a directed input is an `ArgumentError` naming the symmetrisation) | P (every layer must match the first; checked) | R (directed only: rankings are ego-specific; an undirected input is refused) |
| `loops` flag + self-loops | P | P | P | P | P (`allow_self_loops`) | T (keyword on the way back) | D, reported (`:loops`; the ego's own loop is never read, alter self-loops are not recorded) | P (`combine_networks`/`split_by_layer` carry the flag) | D, reported (`:self_loops`; the diagonal is ignored) |
| Vertex set / isolates | P | P (isolates kept; dropped only when `retain_all_vertices=false` and inactive) | P (`n_vertices`) | P (declared actor universe) | P (node-set size) | P | T (only the sampled egos and their alters; `population_size = nv`) | P (`n` checked per layer) | P |
| Vertex IDs | P | P, or T when `retain_all_vertices=false` (renumbered `1:k`; originals in `:vertex_pid`) | P | P | P | P | P (alter IDs are vertex IDs) | P (actor `a` is vertex `(l-1)·n + a` in `combine_networks`) | P |
| Edge set | T (each edge gets one spell) | P | T (one `Contact` per spell) | T (one `Event` per spell **onset**) | P | P | T (ego–alter and alter–alter ties per ego; ties outside the sample D, reported as `:edges` when `n_egos < nv`) | P (`combine_networks`: within-block only; cross-block dyads structurally empty) | T (integer ranks read from the edge attribute `attr=`; every ordered pair must carry one, else R) |
| Vertex attributes | P | P | D | D | D | D (only `:vertex_names`, via `names_eval`) | T (only `ego_attrs=`; every other attribute D, reported as `:vertex_attrs`; an unknown or partially set `ego_attrs` entry is R) | P in the layers; D in `combine_networks` (which adds `:layer`/`:actor`), reported per attribute and layer | D, reported |
| Edge attributes | P | P | D | D | D | T (one attribute at a time, via `attr=`) | D, reported (`:edge_attrs`) | P in the layers; D in `combine_networks`, reported | D, reported (all but `attr`) |
| Network attributes | P | P | D | D | D | D | D, reported (`:network_attrs`) | D in `combine_networks`, reported | D, reported |
| **Missing-dyad mask** | **P** | **P** (remapped; entries on vertices dropped by `retain_all_vertices=false` are D) | **R** (`missing=:error`) | **R** (`missing=:error`) | **R** (`missing=:error`) | **P** with `missing_value=` / `return_mask=` + `missing_mask=`; **D** by default | **R** (`missing=:error`, via `require_observed`); **P at face value** under the written opt-in `missing=:face`, reported as `:missing_dyads` | **P** (`as_multilayer` stores the layers as given; `combine_networks` masks layer `l`'s `(i, j)` at `((l-1)n+i, (l-1)n+j)`; `split_by_layer` remaps it back and **R**-ejects a masked cross-block dyad). No `missing=` keyword: the estimator `ergm_multi` refuses a masked layer | **R** (`missing=:error` only: an unobserved rank has no face value, so there is no `:face`) |
| Two-mode metadata | P | P, but D when `retain_all_vertices=false` (renumbering destroys the mode partition) | D | D | P (`:twomode`) | T (incidence matrix by default; `expand_bipartite=true` gives the square adjacency; flag not recoverable) | R (a `Network` with `bipartite` metadata or a `BipartiteNetwork` is refused) | R (`add_layer!`/`as_multilayer` refuse a two-mode layer: the within-layer dyad universe would enumerate the impossible within-mode pairs as observed non-ties) | R |
| Observation window | P (set to `[onset, terminus)`) | D | D | D | — | — | — | — | — |
| Spells (incl. overlapping) | T (one spell per element) | D (collapsed to presence/absence) | P (one contact per spell; overlaps stay separate) | D (only onsets survive; termini are not events) | — | — | — | — | — |
| Point spells `[t,t)` | — | T (active exactly at `t`) | P (duration 0) | P (an event at `t`) | — | — | — | — | — |
| Spell censoring flags | — | D | D | D — onset-censored spells are *skipped* unless `include_onset_censored=true` | — | — | — | — | — |
| Time-varying (TEA) attributes | — | D | D | D | — | — | — | — | — |
| Actor presence / composition | P (all vertices active for the window) | D | D | D (the risk set is flat over time) | — | — | — | — | — |

Every adapter in the table now honours the contract: a mask is preserved,
rejected, or read at face value only under a policy the caller typed, and
every dropped field is named in a `ConversionReport` under `report=true`.
(Before 0.2.0 the `EgoData` and `MultilayerNetwork` columns read the mask
at face value with no keyword; see [^10] for where the guards live now.)

[^1]: `DynamicNetworks.as_dynamic_network` — **lossless**. A `DynamicNetwork`
    wraps a `Network`, so the whole source object is carried into it by `copy`.
    `network_collapse(as_dynamic_network(net))` reproduces `net`.
[^2]: `DynamicNetworks.network_extract` / `network_collapse` (and
    `TSNA.t_aggregate`, which is `network_collapse` plus an aggregation rule).
    `network_collapse` keeps vertex IDs stable, so it never hits the
    `retain_all_vertices` caveats.
[^3]: `TSNA.as_contact_sequence`.
[^4]: `REM.EventSequence(::DynamicNetwork)` (the `REMDynamicNetworksExt`
    extension).
[^5]: `Siena.DependentNetwork` / `ConstantDyadCovariate` /
    `VaryingDyadCovariate` from `Network` waves (the `SienaNetworkExt`
    extension).
[^6]: `NetworkCore.as_matrix` / `network_from_matrix`. See
    [Missing Data](missing_data.md) for the lossless mask round-trip.
[^7]: `ERGMEgo.simulate_ego_sample(net::Network, n_egos; ego_attrs,
    missing=:error, report=false, rng)` — draws an egocentric sample from a
    complete network. It calls `require_observed(net, missing;
    context="simulate_ego_sample")` first (`supports_missing == true`,
    `missing_policies == (:error, :face)`), refuses directed and two-mode
    input and an unknown or partially set `ego_attrs` entry, and with
    `report=true` returns `(ed, ConversionReport(:Network, :EgoData))`
    recording `:edges`, `:vertex_attrs` (one entry per unrequested
    attribute), `:edge_attrs`, `:network_attrs`, `:loops` and, under
    `:face`, `:missing_dyads`; a census with every vertex attribute requested
    is `is_lossless`. Every cell is pinned by ERGMEgo's testset "Missing
    dyads and the conversion contract (simulate_ego_sample)".
[^8]: `ERGMMulti.as_multilayer(nets, names)` / `add_layer!` (the layers are
    stored by identity, masks and attributes in place; a two-mode layer is
    refused), `ERGMMulti.combine_networks(m; report)` (the block-diagonal
    `Layer()` network: masks preserved at their block position, `loops`
    carried, every layer vertex/edge/network attribute dropped and named in
    the report) and `ERGMMulti.split_by_layer(combined, n, L; report)`
    (masks remapped back; a masked cross-block dyad is an `ArgumentError`;
    the combined network's own attributes dropped and named).
    `supports_missing` is `true` for all three; none takes a `missing=`
    keyword because there is nothing to opt into — the guard stays on
    `ergm_multi`, which calls `require_observed(layer; face_ok=false)` on
    every layer. See ERGMMulti.jl's CHANGELOG, 0.2.0 Added, "The multilayer
    adapters honour the conversion contract".
[^9]: `ERGMRank.as_rank_network(net::Network; attr=:rank, missing=:error,
    report=false)` (and the matrix form `as_rank_network(mat)`, which has no
    mask to guard). `supports_missing == true`, `missing_policies ==
    (:error,)`. Pinned by ERGMRank's testset "Network → RankNetwork adapter
    honours the conversion contract".
[^10]: The two adapters that predated the contract gained `missing=` /
    `report=` with `require_observed` and a `ConversionReport` in 0.2.0
    (see the ERGMEgo and ERGMMulti CHANGELOGs). A masked
    present-face dyad no longer becomes an observed alter–alter tie in an
    `EgoData` (refused, or reported under `missing=:face`) nor an observed
    within-block tie in a combined network (masked there too).

### Paths that are not their own conversion

Some packages do not author a conversion — they *consume* one, and inherit its
row of the table wholesale:

- **NDTV.jl** renders snapshots via
  `network_extract(dnet, t; retain_all_vertices=true)` (its `_snapshot`), so it
  inherits column 2 with stable vertex IDs. It computes layouts and never
  returns a network, so it has no invariants of its own.
- **TSNA.jl's point measures** (`t_degree`, `t_betweenness`, …) do the same.
  The snapshot they measure *carries the missing-dyad mask*, so whether a
  partially observed slice is refused depends on what each measure delegates
  to. Today six of them — `t_degree`, `t_betweenness`, `t_closeness`,
  `t_transitivity`, `t_eigenvector`, `t_pagerank` — refuse through SNA.jl's
  `require_observed` (with SNA's context name in the message), and `t_density`
  and the `:density`/`:mean_degree`/`:n_edges` columns of `t_sna_stats` now
  refuse through `NetworkCore.network_density` (context `network_density`).
  `t_reciprocity`, the temporal-path family (`earliest_arrival!`,
  `forward_reachable_set`, `temporal_path`, `reachability_matrix`,
  `temporal_distance`), `t_edge_duration`, `t_edge_formation` /
  `t_edge_dissolution`, `t_edge_persistence`, `t_turnover` and `tie_decay`
  still compute at **face value**, and no `t_*` measure exposes a `missing=`
  keyword. Giving every `t_*` measure a `missing::Symbol=:error` keyword and
  one behaviour is planned; see `TSNA.jl/CHANGELOG.md` for the status.
- **TSNA.jl's `t_aggregate`** is `network_collapse` plus an aggregation rule
  (`:union` / `:intersection` / `:weighted`), and inherits column 2 exactly.

## Why masked dyads are rejected, not coerced

A masked dyad is **unobserved**, not **absent**. Three of the five adapters
target a representation with no way to say that:

- A `Contact` is a tie that happened at a time. There is no contact meaning
  "we do not know whether this pair ever met."
- An `Event` is an instant. Turning an unobserved dyad into a never-happened
  non-event biases the REM likelihood, which is *conditional on the risk set*.
- A Siena matrix has its own mask — but it records **structural** zeros and
  ones, i.e. ties that are *determined*. Coding an unobserved dyad as a
  structural zero tells the estimator the tie is known to be impossible. That
  is a different, and false, claim.

None of these has a faithful encoding, so all three raise on a masked input:

```julia-repl
julia> as_contact_sequence(dnet)
ERROR: ArgumentError: as_contact_sequence does not support missing (unobserved)
dyads, but the network has 2 masked dyads. A masked dyad is unobserved, not
absent, so reading its face value would silently invent data.
  • pass `missing=:face` to analyse the stored face values anyway …
```

The three ways forward are always the same: pass `missing=:face` to accept the
face values in writing, call `clear_missing_dyads!` to declare the dyads
observed, or use a routine that handles missingness
(`NetworkCore.supports_missing(f) == true`).

Where the target *can* hold a mask — `Network` ↔ `DynamicNetwork` in both
directions — it is carried, not dropped. An unobserved dyad of a base network
is unobserved in every snapshot of it.

## API

```@docs
ConversionReport
record_drop!
is_lossless
dropped_fields
```
