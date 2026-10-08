# Graph Interface API Reference

This page documents the Graphs.jl `AbstractGraph` interface methods implemented by NetworkCore.jl. These functions let a `Network` be passed to Graphs.jl algorithms (`Graphs.core_number(net)`, `Graphs.betweenness_centrality(net)`, `Graphs.connected_components(net)`, …), which then see the network's true directedness.

!!! warning "Pass the `Network`, never `net.graph`"
    An undirected `Network` is stored as a `SimpleDiGraph` holding **both
    arcs** of every edge (the field `net.graph`). The interface methods below
    account for that; Graphs.jl itself does not. A Graphs.jl algorithm called
    on `net.graph` therefore sees a directed graph with every tie twice:

    | Florentine marriage network | on `flo` | on `flo.graph` |
    |---|---|---|
    | `core_number` | 1, 2, 2, … (max 2) | 2, 4, 4, … (doubled) |
    | maximum `degree` | 6 | 12 |
    | `global_clustering_coefficient` | 0.1915 | 0.0395 |
    | `maximal_cliques` | works | `MethodError` (directed) |

    Always call Graphs.jl on the `Network` itself. When an algorithm needs a
    concrete Graphs.jl type, convert explicitly: `Graphs.SimpleGraph(net)` is
    an undirected copy (each edge once) and `Graphs.SimpleDiGraph(net)` a
    directed copy of the storage. Never mutate `net.graph` either — the
    mutators keep edge attributes, the missing-dyad mask, the two-mode
    partition and the cached self-loop count in step with it.

```julia
using NetworkCore, Graphs
flo = load_dataset(:florentine_marriage)
core_number(flo) == core_number(SimpleGraph(flo))    # true — the right answer
core_number(flo.graph) == 2 .* core_number(flo)      # true — the trap
```

## Basic Properties

### nv

```@docs
nv
```

### ne

```@docs
ne
```

### vertices

```@docs
vertices
```

### edges

```@docs
edges
src
dst
```

### has\_vertex

```@docs
has_vertex
```

### has\_edge

```@docs
has_edge
```

### degree, indegree, outdegree

The Graphs.jl degree family, re-exported so that `using NetworkCore` alone
answers the first descriptive question a user asks. On an undirected network
all three agree; on a directed one `degree` is `indegree + outdegree`. A
self-loop counts once on an undirected network and twice (one in, one out) on
a directed one, as in Graphs.jl.

```@docs
degree
indegree
outdegree
```

### is\_directed

```@docs
is_directed
```

### is\_bipartite

`is_bipartite` is Graphs.jl's graph-theoretic two-colourability test, which
NetworkCore.jl extends for `Network`. It is *not* the statnet two-mode metadata
flag — for that, see [`is_two_mode`](@ref) below.

```@docs
Graphs.is_bipartite
```

### is\_two\_mode

```@docs
is_two_mode
```

### zero

```@docs
Base.zero
```

## Vertex Modification

### add\_vertex!

```@docs
add_vertex!
```

### add\_vertices!

```@docs
add_vertices!
```

### rem\_vertex!

```@docs
rem_vertex!
```

### rem\_vertices!

```@docs
Graphs.rem_vertices!(::Network, ::AbstractVector{<:Integer})
```

## Edge Modification

### add\_edge!

```@docs
add_edge!
```

### add\_edges!

```@docs
add_edges!
```

### rem\_edge!

```@docs
rem_edge!
```

## Neighbors

### neighbors

```@docs
neighbors
```

### inneighbors

```@docs
inneighbors
```

### outneighbors

```@docs
outneighbors
```

## Utility Functions

### network\_size

```@docs
network_size
```

### network\_edgecount

```@docs
network_edgecount
```

### network\_density

```@docs
network_density
```

### get\_neighborhood

```@docs
get_neighborhood
```

### get\_induced\_subgraph

```@docs
get_induced_subgraph
```

### permute\_vertices

```@docs
permute_vertices
```

## Graphs.jl constructions

Graphs.jl generics that build a new graph from a `Network`. NetworkCore.jl adds
methods so the result is a `Network` with its attributes, missing-dyad mask,
`loops` flag and two-mode partition — not a bare graph.

### induced\_subgraph

```@docs
Graphs.induced_subgraph(::Network, ::AbstractVector{<:Integer})
```

### complement

```@docs
Graphs.complement(::Network)
```
