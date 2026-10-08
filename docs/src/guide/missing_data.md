# Missing Data (Unobserved Dyads)

Real network data are often incompletely observed: a survey respondent skips
the roster question, an archival record is illegible, a dyad simply was not
measured. Treating such dyads as "no tie" biases every statistic computed from
the network — R's statnet marks them as `NA` edges instead. NetworkCore.jl
provides the analogous facility: a **missing-dyad mask**.

## Semantics

A masked dyad means the tie status is **unobserved** — this is distinct from
"no tie":

- The backing graph always stores a *face value* for every dyad (edge present
  or absent). Structural queries — `has_edge`, `ne`, `degree`, `as_matrix`,
  and every other Graphs.jl generic — continue to report the face value and
  do **not** consult the mask.
- The mask is metadata for statistical procedures. ERGM.jl, for example,
  excludes masked dyads from the MPLE design matrix (they are not observed
  responses) and never toggles them during MCMC simulation.

For undirected networks the mask is symmetric (`(i, j)` and `(j, i)` are the
same dyad); for directed networks each ordered pair (arc) is masked
independently. Lookups are O(1), so samplers can consult the mask in hot
loops.

## API

```julia
using NetworkCore

net = network(5)
add_edge!(net, 1, 2)

set_missing_dyad!(net, 3, 4)      # mark 3→4 as unobserved
set_missing_dyad!(net, 1, 2)      # an existing edge can be unobserved too

is_missing_dyad(net, 3, 4)        # true
has_edge(net, 3, 4)               # false — face value, NOT an observation
has_edge(net, 1, 2)               # true  — face value

n_missing_dyads(net)              # 2
collect(missing_dyads(net))       # [(1, 2), (3, 4)] (order unspecified)

delete_missing_dyad!(net, 3, 4)   # declare it observed again -> true
clear_missing_dyads!(net)         # unmask everything
```

The mask survives `copy`, is remapped by `rem_vertex!` (swap-and-pop on a
one-mode network, stable renumbering on a two-mode one), `rem_vertices!`,
`get_induced_subgraph`, and `permute_vertices`, and is forwarded by
`BipartiteNetwork` (within-mode dyads of two-mode networks cannot be masked —
they are structurally impossible, not unobserved).

## Guarded routines in NetworkCore.jl

NetworkCore.jl itself computes one statistic, and it honours the contract:

- **`network_density`** refuses a masked network by default and takes
  `missing=:face` to count every edge at its stored face value over the full
  dyad count (R's `network.density(x, na.omit=FALSE)`). R's default
  `na.omit=TRUE` — drop the NA edges from the numerator, keep the denominator —
  is deliberately *not* offered, because it reads an unobserved tie as absent.
  `missing_policies(network_density) == (:error, :face)`.
- The **conversion adapters** are lossless on request rather than guarded:
  `as_matrix(net; missing_value=missing)` writes `missing` at masked dyads and
  `as_matrix(net; return_mask=true)` returns the mask beside the matrix, so the
  round trip through `network_from_matrix` preserves both the mask and the face
  values (see [Conversion Invariants](conversion_invariants.md)).

Every SNA.jl descriptive measure, every conversion adapter into another
representation (`DynamicNetwork` → snapshots, contact sequences, event
sequences, Siena matrices; and, since 0.2.0, `Network` →
`EgoData`, `RankNetwork` and the multilayer adapters — see the invariant
table) and `network_density` refuse a masked network unless
asked in writing; ERGM.jl's MPLE drops masked dyads from the pseudo-likelihood.

## Limitations

- The mask is **not persisted** by the I/O routines (Pajek, GraphML, CSV);
  those formats have no missing-tie representation.
- The *default* `as_matrix(net)` / `as_edgelist(net)` carry face values only
  (the mask is invisible in a plain `Matrix{Float64}`); use `missing_value=` or
  `return_mask=` when the mask must survive the trip. Face values reached this
  way are what the caller asked for, in writing.
- Downstream estimation support is deliberately scoped: ERGM.jl's MPLE
  excludes masked dyads exactly, while its MCMLE *conditions* on their face
  values only under the explicit `missing=:condition_on_face` policy (`:face`
  is not accepted there). Full statnet-style missing-data maximum likelihood
  is future work.
