# Changelog

All notable changes to NetworkCore.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

First public release. NetworkCore.jl provides the `Network` type (a Graphs.jl
`AbstractGraph` with vertex, edge and network attributes and a missing-dyad
mask), conversions and I/O, six bundled datasets, and the contracts, numerics
and result-presentation layer shared by the model packages built on it.

Entries marked **Changed** differ from the 0.1.0 development version, which
was never published; they matter only to code written against it.

### Renamed

- **The package and module are now `NetworkCore`** (directory and repository
  `NetworkCore.jl`); the package was previously developed as `Networks.jl`.
  The name `Networks` is registered in Julia's General registry by an unrelated
  package, so this package could never be registered under it. Write
  `using NetworkCore` / `import NetworkCore: …` where code said
  `using Networks` / `import Networks: …`, and `NetworkCore.foo` for qualified
  names. The UUID is unchanged, and so are the type and function names
  (`Network`, `BipartiteNetwork`, `network`, …).

### Network type and Graphs.jl interface

- `Network{T,D}`: directedness `D` is a type parameter, so
  `Graphs.is_directed` is truthful for instances and types and Graphs.jl's
  trait-dispatched algorithms use the right semantics. `Network{Int}` matches
  both directednesses and `Network{Int}(; n, directed)` picks `D`.
  `net.directed` is read-only. **Changed:** it was a mutable field, and
  `is_directed` returned `true` for every network.
- `BipartiteNetwork{T,D}` carries the same parameter, so its `network` field
  is concretely typed. `BipartiteNetwork(n₁, n₂; directed)` and
  `BipartiteNetwork{T}(…)` construct as before.
- Any integer vertex type works (`Network{Int32}`), and constructors accept
  any `Integer` vertex count. A negative count, `multiple=true` and
  `hyper=true` are `ArgumentError`s.
- `show(io, net)` prints one line (`Network{Int64,true}(50 vertices, 113
  edges)`), so containers of networks stay readable. The REPL shows an
  R-style block with directedness, loops, the two-mode sizes and sorted
  attribute names. **Changed:** code that parsed `string(net)` for
  `"Vertices:"` must use `repr(MIME"text/plain"(), net)`.
- `copy(net)` duplicates the graph, all attributes and the missing-dyad mask;
  `zero` is defined for network types and instances.
- Re-exported from Graphs.jl: `nv`, `ne`, `vertices`, `edges`, `src`, `dst`,
  `neighbors`, `inneighbors`, `outneighbors`, `degree`, `indegree`,
  `outdegree`, `has_edge`, `has_vertex`, `is_directed` and the mutators.
- `Graphs.induced_subgraph(net, vlist)`, `Graphs.complement(net)` (R's
  `!net`) and `Graphs.rem_vertices!(net, vs)` (R's `delete.vertices`, stable
  renumbering) return networks that keep attributes, the mask and the
  two-mode partition. `Graphs.SimpleGraph(net)` and `Graphs.SimpleDiGraph(net)`
  convert to plain Graphs.jl graphs.
- `Graphs.is_bipartite(net)` is the graph-theoretic two-colourability test;
  `is_two_mode(net)` reports the two-mode construction flag. **Changed:**
  `is_bipartite` used to return the flag.
- Undirected networks store both arcs in `net.graph`. Call Graphs.jl
  algorithms on the `Network`, never on `net.graph`; the documentation shows
  the wrong answers the field gives.

### Editing networks

- `add_edges!` and `network_from_edgelist` never drop an edge: an endpoint
  outside `1:nv(net)`, a self-loop under `loops=false` or a within-mode pair
  on a two-mode network is an `ArgumentError` naming the edge. A repeated edge
  counts zero. The single `add_edge!` keeps Graphs.jl's `Bool` return.
- `rem_edge!` removes the edge's attributes with it and is a no-op on an
  absent edge.
- `rem_vertex!` keeps attributes and mask attached to the right vertices. On a
  one-mode network it follows Graphs.jl's swap-and-pop; on a two-mode network
  it renumbers stably and lowers `bipartite` when a mode-1 vertex is removed,
  as R does.
- `get_induced_subgraph` refuses an out-of-range or repeated vertex, keeps the
  two-mode partition, and refuses a selection that interleaves the two modes.
- `permute_vertices` refuses a vector that is not a permutation, and on a
  two-mode network a permutation that moves a vertex across the modes.
- Two-mode networks reject within-mode edges.

### Attributes

- Vertex, edge and network attributes with `get_*`, `set_*!`, `delete_*!`,
  `list_*` and the `net[:v, :attr]` shorthand.
- Setters validate their targets: a vertex outside `1:nv(net)` is an
  `ArgumentError`, and an edge attribute needs an existing edge unless
  `require_edge=false`. The `Dict` forms check every key before writing any.
- Typed accessors return concrete containers:
  `get_vertex_attribute(net, attr, T)`, `get_edge_attribute(net, attr, T)` and
  `vertex_attribute_vector(net, attr, T; default)`. The last throws an
  `ArgumentError` naming the attribute, the vertex and the attributes that
  exist when a value is absent and no `default` is given.

### Missing data

- Missing-dyad mask: `set_missing_dyad!`, `is_missing_dyad`,
  `delete_missing_dyad!`, `clear_missing_dyads!`, `missing_dyads`,
  `n_missing_dyads`. A masked dyad is unobserved, not absent. The mask
  survives `copy` and every vertex edit.
- The missing-data contract for the ecosystem: `supports_missing(f)`,
  `require_observed(net, policy; context, hint)`, `MISSING_POLICIES` and
  `missing_policies(f)` (also `missing_policies(f, ResultType)` for shared
  verbs such as `gof`). A routine refuses a masked network unless the caller
  passes `missing=:face`.
- `network_density(net; missing=:error, diag=net.loops,
  discount_bipartite=false)` follows R's `network.density`, pinned by a golden
  fixture. It refuses a masked network by default and returns `NaN` below two
  vertices. `Network` and `BipartiteNetwork` share one implementation with R's
  full-dyad denominator; `discount_bipartite=true` uses the cross-mode dyads.

### Conversion and I/O

- `as_matrix`, `as_adjacency_matrix`, `as_edgelist`, `network_from_matrix`,
  `network_from_edgelist`; `as_dataframe` and `network_from_dataframe` through
  a DataFrames.jl package extension (calling them without DataFrames loaded
  says so).
- `as_matrix(net; missing_value=missing)` and `return_mask=true` carry the
  mask; with `network_from_matrix(A; missing_mask=)` the round-trip is
  lossless. `network_from_matrix` takes `store_values=` and warns on an
  asymmetric matrix with `directed=false`.
- Two-mode networks convert to the `n₁ × n₂` incidence matrix, or to the
  square adjacency with `expand_bipartite=true`. On a directed two-mode
  network the incidence matrix holds mode-1 → mode-2 arcs only; the others are
  dropped with a warning, as R drops them.
- `as_matrix(net; attr=)` and `as_edgelist(net; attrs=)` refuse an edge
  attribute the network does not have, listing the ones it has.
- `read_pajek` / `write_pajek`: `*Arcs`, `*Edges`, `*Arcslist`, `*Edgeslist`,
  weights, self-loops and the two-mode header `*Vertices n n₁`. A tie the
  network cannot hold is an `ArgumentError`, a repeated tie warns, a mixed
  `*Arcs`/`*Edges` file loads as directed with a warning, and `*Matrix` is
  refused. `write_graphml` and `write_edgelist_csv` (RFC 4180 quoting) write
  only.
- The conversion contract for adapters in other packages: `ConversionReport`,
  `record_drop!`, `is_lossless`, `dropped_fields`. A lossy conversion takes
  `report=true` and names what it dropped. The guide publishes the invariant
  table for every conversion path.

### Datasets

- `load_dataset(name)` for `:florentine_marriage`, `:florentine_business`,
  `:sampson` and `:faux_mesa_high` (networks), `:s50` (three friendship waves
  with alcohol and smoking) and `:wtc_police_calls` (481 ordinal events among
  37 actors). The last two return documented `NamedTuple`s.
- `data/README.md` records source, file hashes, upstream version, licence and
  citation for each dataset, and `data/r/extract_datasets.R` regenerates the
  files. The data come from GPL-licensed R packages; the README tabulates the
  upstream licences.
- The data directory is resolved at run time, so loading works from a
  relocated depot or a system image.

### Shared layer for model packages

- **Statistic protocol:** empty generics `compute`, `name`, `compute_all` and
  `gof`, which every model package extends, so loading several packages never
  leaves a verb undefined.
- **Result presentation:** `print_coeftable` (optional standard-error column,
  `p_floor=` for permutation tests), `CoefficientTable`, `format_pvalue`,
  `signif_code`, `SIGNIF_LEGEND`, `GOFStatistic`, `GOFResult`,
  `n_simulations`.
- **p-values:** `z_pvalues` computes two-sided normal p-values through `erfc`,
  floored at `floatmin` for finite `z` and `NaN`-aware. `mc_pvalue` is the
  `(1 + k)/(N + 1)` Monte-Carlo p-value; it returns `NaN` for an undefined
  observed statistic or when no simulated value is defined, and drops `NaN`
  draws.
- **StatsAPI:** `coeftable` and `coefnames` are re-exported from StatsAPI,
  and `check_statsapi(fit; required, strict)` reports which StatsAPI verbs a
  fitted-model type implements consistently. `coefnames` is checked whenever
  a fit defines it (one string per coefficient, equal to
  `coeftable(fit).names`) and required when listed in `required=`.
- **Result metadata:** `ResultMetadata`, `fit_metadata`, `estimand`,
  `objective`, `is_exact`, `se_method`, `missing_method`, `tie_method`,
  `approximations`.
- **Tied events:** `TIE_POLICIES` and `check_tie_policy` give relational-event
  packages one `ties=` vocabulary; an unsupported policy fails with the reason.
- **Numerics:**
  - `newton_fit` (public, not exported): Newton–Raphson with up to 30 step
    halvings. A step with a non-finite objective is never taken. Convergence
    is judged by the Newton decrement, so the verdict does not depend on the
    scale of the design. Standard errors are `NaN`, with a warning and
    `converged=false`, when the information matrix is not positive definite
    or not numerically identified (`information_rtol=`). Errors raised by the
    objective propagate. The result's `stop` field says why the iteration
    ended (`:converged`, `:maxiter`, `:singular_start`, `:singular`,
    `:no_ascent`, `:singular_information`), so a fit whose information is
    singular at the starting value is no longer indistinguishable from one
    that ran out of iterations. `converged` does not certify a finite maximiser; the
    separation verdict below does.
  - `logistic_derivatives` (public, not exported): the logistic
    log-likelihood, gradient and Hessian for binary rows (with `weights=` and
    `offset=`) or binomial rows; an evaluation allocates only what it returns.
  - `bootstrap_cov(refit, simulate, θ̂; n_boot, rng, threaded)`: the shared
    resampling loop, reproducible from `rng` and independent of the thread
    count. An exception thrown by `refit` propagates as itself, not wrapped
    in a `CompositeException`/`TaskFailedException`.
  - `spawn_all(f, n)` and `unwrap_task_exception(e)` (public, not
    exported): run `f(1), …, f(n)` on tasks, wait for all of them, and
    rethrow the first failure as the task's own exception. The
    `TaskFailedException` stays on the exception stack as its cause, so the
    backtrace from inside the failing task is printed and available through
    `current_exceptions()`. One helper for every threaded sampler.
- **Separation:** one verdict and one policy for every fitting function
  (public, not exported).
  - `logistic_separation` (binary or binomial rows), `poisson_separation`
    and `clogit_separation` (strata of alternatives, one or several chosen
    per stratum) decide whether the likelihood has a finite maximiser.
    `separation_from_margins` is the criterion they share. The verdict is a
    property of the data, not of where the optimizer stopped. It is decided
    by a linear programme and certified in exact arithmetic. A
    `SeparationVerdict` names an inclusion-minimal set of separated
    coefficients and the observations predicted perfectly.
  - The verdict's cost is that of a few floating-point linear programmes,
    separated or not: the direction is reduced to its minimal support before
    it is certified, and the certificate uses exact integer arithmetic only on
    margins that are numerically zero. A quasi-separated relational-event
    design of 640,000 margin rows takes well under a second and a few hundred
    MB. Proving separation never needs more than 200,000 big-integer row
    evaluations; past that the verdict reports `certified=false`.
  - The policy: warn (`separation_message`, `warn_separation`), report
    `converged=false`, flag the separated terms, and withhold z values,
    p-values and confidence intervals (`z_pvalues(est, se, verdict)` returns
    `NaN`; `separation_caveat` is the `approximations` entry). The warning
    names the objective after the estimate: "log-pseudo-likelihood" for an
    MPLE, CMPLE or swap-MPLE, "log-likelihood" otherwise.
  - `check_se(se, allowed; context)` validates an `se=` keyword.
- **Golden fixtures:** `load_golden`, `check_golden`, `golden_report`,
  `golden_tolerance`, `GoldenFixture`. A fixture without a `[provenance]`
  block, or whose generator script cannot be found, is refused.

### Performance

- `rem_edge!` allocates nothing; `add_edge!` allocates what Graphs.jl's own
  `add_edge!` does.
- `ne` is O(1), including on undirected networks with self-loops.
- `network_density` allocates nothing, and the whole-`Dict` attribute getters
  allocate nothing when the attribute exists.

### Requirements and quality

- Julia 1.12 or later. Dependencies beyond the standard library: Graphs,
  SpecialFunctions, StatsAPI; DataFrames is optional.
- Every exported and public name has a docstring with a runnable example, and
  the test suite checks this, freezes the export list and runs Aqua.jl.

### Known limitations

Kept in sync with the README "Not implemented" section.

- Multiplex (parallel) edges and hyperedges are not supported
  (`multiple=true`/`hyper=true` throw).
- The missing-dyad mask is not persisted by Pajek, GraphML or CSV output; use
  `as_matrix(net; return_mask=true)` / `network_from_matrix(A; missing_mask=)`.
- `read_pajek` refuses `*Matrix` sections; there is no GraphML reader.
- `network_density` does not offer R's `na.omit=TRUE` (unobserved read as
  absent); only `missing=:error` and `missing=:face`.
- `get_induced_subgraph` numbers vertices in the order of `vlist`, not in
  their original order as R does; pass `sort(vlist)` to match R. Interleaved
  two-mode selections are refused.
- `permute_vertices` refuses permutations that move a vertex across the modes
  of a two-mode network (R warns and proceeds).
- `newton_fit`'s `converged` does not certify a finite maximiser; the
  separation verdict (`logistic_separation`, `poisson_separation`,
  `clogit_separation`) does.

## [0.1.0] - 2026-02-09

Development version, never published.
