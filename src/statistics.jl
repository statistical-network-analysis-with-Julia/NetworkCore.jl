# The shared statistic protocol.
#
# Every model package in the ecosystem has a notion of a *statistic* (ERGM.jl's
# terms, TERGM.jl's temporal terms, REM.jl's relational-event statistics,
# ERGMCount/ERGMMulti/ERGMRank/ERGMEgo/ERGMUserterms' term families, Revel's
# relational event effects). Each one needs the same two verbs:
#
#   compute(stat, ...) -> Float64   evaluate the statistic
#   name(stat)         -> String    the label it prints under
#
# plus a `compute_all(set, ...)` sweep over a collection of them.
#
# Those verbs live HERE, as empty generics, for the same reason `gof` does
# (see results.jl): if each package defined and exported its own `compute` and
# `name`, the bindings would be *different functions with the same name*, and
# Julia's rule for conflicting exports would leave `compute` and `name`
# UNDEFINED in a session that loaded two of them —
#
#   using ERGM, REM      # cross-sections + dynamics: the core statnet workflow
#   compute              # UndefVarError, before this file existed
#
# Because every package instead extends this ONE generic, `compute(term, net)`
# and `compute(stat, state, sender, receiver)` are simply methods of the same
# function, dispatched on their (different) signatures. Nothing is ambiguous,
# nothing has to be qualified, and users can add their own statistics with
# `import NetworkCore: compute, name` from anywhere.
#
# The generics are deliberately empty: NetworkCore.jl knows nothing about what a
# statistic *is*. Each package owns its own abstract statistic type, its own
# fallbacks, and its own methods.

"""
    compute(stat, args...) -> Float64

Evaluate a network statistic. The ecosystem-wide generic for the "value of this
statistic" verb; every model package adds methods for its own statistic types:

| package    | signature                                       |
|:-----------|:------------------------------------------------|
| ERGM.jl    | `compute(term::AbstractERGMTerm, net)`          |
| TERGM.jl   | `compute(term::TemporalTerm, net, prev_net)`    |
| REM.jl     | `compute(stat::AbstractStatistic, state, sender, receiver)` |

The methods have different signatures, so they coexist on one generic without
ambiguity — which is the point: `using ERGM, REM` leaves `compute` usable
unqualified, instead of undefined by Julia's conflicting-export rule.

Extend it for your own statistic with `import NetworkCore: compute`.

# Example
```julia
import NetworkCore: compute, name
struct EdgeCount end
compute(::EdgeCount, net::Network) = Float64(ne(net))
name(::EdgeCount) = "edges"

net = network(3); add_edges!(net, [(1, 2), (2, 3)])
compute(EdgeCount(), net)    # 2.0
```

See also [`name`](@ref), [`compute_all`](@ref).
"""
function compute end

"""
    name(stat) -> String

The display label of a network statistic — the string that identifies it in a
coefficient table. The ecosystem-wide generic; every model package adds methods
for its own statistic types (see [`compute`](@ref)).

Extend it for your own statistic with `import NetworkCore: name`.

# Example
```julia
import NetworkCore: name
struct EdgeCount end
name(::EdgeCount) = "edges"
name(EdgeCount())    # "edges"
```
"""
function name end

"""
    compute_all(set, args...) -> Vector{Float64}

Evaluate every statistic in a set, returning the vector of values in set order.
The ecosystem-wide generic; model packages add methods for their own set types
(ERGM.jl's `TermSet`, REM.jl's `StatisticSet`, ...).

# Example
```julia
import NetworkCore: compute, compute_all
struct EdgeCount end
struct VertexCount end
compute(::EdgeCount, net::Network) = Float64(ne(net))
compute(::VertexCount, net::Network) = Float64(nv(net))
compute_all(stats::Vector, net::Network) = [compute(s, net) for s in stats]

net = network(3); add_edge!(net, 1, 2)
compute_all(Any[EdgeCount(), VertexCount()], net)    # [1.0, 3.0]
```

See also [`compute`](@ref).
"""
function compute_all end
