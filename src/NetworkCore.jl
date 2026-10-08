"""
    NetworkCore.jl - Network Data Structures for Julia

A Julia package providing network data structures compatible with the StatNet
ecosystem. Implements the Graphs.jl AbstractGraph interface while providing
StatNet-compatible features like vertex/edge attributes.

This is the foundational package upon which SNA.jl, ERGM.jl, and other
Statistical Network Analysis packages are built.

The **module** is `NetworkCore` and the **type** it exports is `Network`, so the
obvious thing works:

```julia
using NetworkCore
net = Network(5)          # the documented constructor, callable
```

The module was called `Network` until v0.2. Because a module and its type
cannot share a name in Julia, the exported type shadowed the module: `Network(5)`
raised `MethodError: objects of type Module are not callable`, qualified access
like `NetworkCore.require_observed(...)` resolved to *field access on the type* and
failed to precompile in dependent packages, and the docsystem needed a
`Base.Docs.doc!` workaround. Renaming the module removes all three problems and
leaves the type name — which appears in ~200 downstream signatures — untouched.
"""
module NetworkCore

using Graphs
using SparseArrays
using Printf
using Random
using LinearAlgebra: diag, dot, cholesky, Symmetric, mul!, norm, issuccess,
                     SingularException, LAPACKException, ZeroPivotException, PosDefException
using SpecialFunctions: erfc
using Statistics: mean, quantile, cov
# The StatsAPI generics. NetworkCore.jl defines NO method on any of them (it has
# no fitted-model type); it re-exports `coeftable` and `coefnames` so that every
# model package extends the ONE binding (`NetworkCore.coeftable ===
# StatsAPI.coeftable`, likewise `coefnames`) and `check_statsapi` can ask
# `hasmethod` about the rest.
import StatsAPI
import StatsAPI: coeftable, coefnames

# Core types
export AbstractNetwork, Network, BipartiteNetwork

# Graph interface (re-exported from Graphs.jl for convenience). `degree`,
# `indegree` and `outdegree` are the Graphs.jl bindings, re-exported so that a
# plain `using NetworkCore` (and the packages that mirror this list, ERGM/SNA)
# answers the first descriptive question a user asks.
export nv, ne, vertices, edges, has_vertex, has_edge
export neighbors, inneighbors, outneighbors
export degree, indegree, outdegree
export is_directed
# `edges(net)` yields Graphs.jl edges, whose endpoints are read with `src`/`dst`;
# re-exporting the iterator without its accessors left every term author who
# wrote `for e in edges(net); src(e) ...` with an `UndefVarError` inside
# `compute` (ERGMUserterms).
export src, dst

# Network construction
export network, network_initialize
export add_vertex!, add_vertices!, rem_vertex!
export add_edge!, add_edges!, rem_edge!

# Missing-dyad (unobserved tie) mask
export set_missing_dyad!, is_missing_dyad, delete_missing_dyad!
export clear_missing_dyads!, missing_dyads, n_missing_dyads

# Missing-data contract (see missing.jl): the trait a routine opts into, the
# guard it calls so unobserved ties are never read at face value by default,
# and the routine-level vocabulary (`missing_policies(f)`) that tooling such as
# the capability matrix prints instead of assuming `:face` everywhere
export supports_missing, require_observed, MISSING_POLICIES, missing_policies

# Attribute handling
export get_vertex_attribute, set_vertex_attribute!, delete_vertex_attribute!
export vertex_attribute_vector
export get_edge_attribute, set_edge_attribute!, delete_edge_attribute!
export get_network_attribute, set_network_attribute!, delete_network_attribute!
export list_vertex_attributes, list_edge_attributes, list_network_attributes

# Coercion and conversion
export as_matrix, as_adjacency_matrix, as_edgelist, as_dataframe
export network_from_matrix, network_from_edgelist

# Conversion contract (see conversion.jl): the metadata a lossy conversion
# returns under `report=true`, naming every field it could not carry across
export ConversionReport, record_drop!, is_lossless, dropped_fields

# I/O
export read_pajek, write_pajek, write_graphml, write_edgelist_csv

# Bundled teaching datasets
export load_dataset

# DataFrame conversion (methods provided by the NetworkDataFramesExt
# package extension; requires `using DataFrames`)
export network_from_dataframe

# Shared statistic protocol (see statistics.jl): the empty `compute`/`name`/
# `compute_all` generics that EVERY model package extends for its own statistic
# types, so that `using ERGM, REM` leaves the verbs usable unqualified instead
# of undefined by Julia's conflicting-export rule
export compute, name, compute_all

# Shared result presentation (used by the model packages' show methods).
# `z_pvalues` is the ONE z → p helper (floored, NaN-aware) behind every
# Pr(>|z|) column in the ecosystem.
export print_coeftable, format_pvalue, signif_code, SIGNIF_LEGEND, z_pvalues
export GOFStatistic, GOFResult, n_simulations, mc_pvalue, gof
# The generic coefficient table every `coeftable(fit)` returns, the StatsAPI
# verbs `coeftable` and `coefnames` themselves (same bindings as StatsAPI's),
# and the checker of the StatsAPI surface (see statsapi.jl).
export CoefficientTable, coeftable, coefnames, check_statsapi

# Shared result-metadata protocol (see results.jl): what a fit ACTUALLY did —
# which objective, exact or not, which SEs, how missing dyads and ties were
# handled. Model packages add methods for their own result types.
export ResultMetadata, fit_metadata
export estimand, objective, is_exact, se_method, missing_method, tie_method
export approximations
# The shared tie-handling vocabulary (`ties=`): one set of symbols, one meaning
# per symbol, for every model with an event clock; `check_tie_policy` refuses a
# policy a model cannot honour instead of letting it no-op.
export TIE_POLICIES, check_tie_policy
# The shared `se=` validator (see bootstrap.jl): one message shape for every
# fitter that offers a subset of (:hessian, :sandwich, :bootstrap).
export check_se

# Shared resampling infrastructure (see bootstrap.jl): the ONE bootstrap loop
# behind every `se=:bootstrap` option in the ecosystem — the parametric
# bootstrap of the pseudo-likelihood estimators (ERGM, ERGMCount, ERGMRank,
# ERGMMulti) and REM's repeated draw of the case-control risk set. It lives here
# rather than in ERGM.jl because REM.jl does not depend on ERGM.jl.
export bootstrap_cov

# Golden-fixture harness (see golden.jl): frozen reference-implementation output
# WITH provenance. Loading a fixture that cannot be regenerated is an error.
export GoldenFixture, load_golden, check_golden, golden_report, golden_tolerance

# Utilities
export network_size, network_density, network_edgecount
export is_two_mode
export permute_vertices, get_neighborhood, get_induced_subgraph

# ----------------------------------------------------------------------------
# Public API that is deliberately NOT exported (Julia ≥ 1.11 `public`).
#
# These names are part of the supported surface — documented, semver-covered,
# and safe to `import NetworkCore: <name>` from a dependant — but they are not
# pulled into the caller's namespace by `using NetworkCore`, either because they
# would collide with an existing binding in a dependant (ERGM.jl's test suite
# does `using ERGM; using NetworkCore` and calls `newton_fit` unqualified, so a
# second export of a different binding would make the name ambiguous there
# until ERGM aliases `import NetworkCore: newton_fit`), or because they are
# tooling that a user analysing a network never types.
#
# The shared optimizer (`newton_fit`, `logistic_derivatives`; see newton.jl)
# is hosted here and declared `public` — ERGM.jl re-exports it (`import
# NetworkCore: newton_fit, logistic_derivatives` + `export`), so `using ERGM`
# keeps working while there is exactly ONE definition. Note that `names(M)`
# lists public names too (Julia ≥ 1.11): the "Exported API inventory" testset
# freezes the exported and the public lists separately.
# ----------------------------------------------------------------------------
public newton_fit, logistic_derivatives
# The separation verdict and its policy (see separation.jl): one detector per
# likelihood family over ONE criterion, and the one warning, caveat and
# withheld inference every fitting function applies to a separated fit.
public SeparationVerdict, separation_from_margins
public logistic_separation, poisson_separation, clogit_separation
public separation_message, warn_separation, separation_caveat
# Shared task helpers (see tasks.jl): a threaded loop whose errors reach the
# caller as the task's own exception.
public spawn_all, unwrap_task_exception

# Include source files
include("types.jl")
include("attributes.jl")
include("missing.jl")
include("conversion.jl")
include("graphs_interface.jl")
include("coercion.jl")
include("io.jl")
include("datasets.jl")
include("statistics.jl")
include("results.jl")
include("statsapi.jl")
include("bootstrap.jl")
include("newton.jl")
include("separation.jl")
include("tasks.jl")
include("golden.jl")

# `as_dataframe` / `network_from_dataframe` are empty generics here; their
# methods live in the DataFrames package extension. Without DataFrames loaded a
# call is a bare `MethodError`, so say what is missing.
function _dataframes_error_hint(io::IO, exc::MethodError, argtypes, kwargs)
    if (exc.f === as_dataframe || exc.f === network_from_dataframe) &&
       Base.get_extension(@__MODULE__, :NetworkDataFramesExt) === nothing
        print(io, "\n", nameof(exc.f), " is provided by NetworkCore.jl's DataFrames ",
              "extension, which is not loaded: run `using DataFrames` ",
              "(after `Pkg.add(\"DataFrames\")` if needed) and call it again.")
    end
    return nothing
end

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(_dataframes_error_hint, MethodError)
    end
    return nothing
end

end # module
