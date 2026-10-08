# Core network types for the NetworkCore.jl package.
#
# Provides network data structures that implement Graphs.jl's AbstractGraph
# interface while supporting vertex/edge attributes like R's network package.

"""
    AbstractNetwork{T} <: Graphs.AbstractGraph{T}

Abstract base type for all network types. Extends Graphs.jl's AbstractGraph
to provide StatNet-compatible network functionality. [`Network`](@ref) and
[`BipartiteNetwork`](@ref) are its two concrete subtypes; write a method
against `AbstractNetwork` to accept either.

# Example
```julia
Network(3) isa AbstractNetwork              # true
BipartiteNetwork(2, 3) isa AbstractNetwork  # true
n_dyads(net::AbstractNetwork) = nv(net) * (nv(net) - 1)
n_dyads(Network(4))                         # 12
```
"""
abstract type AbstractNetwork{T<:Integer} <: Graphs.AbstractGraph{T} end

"""
    Network{T,D}

Core network data structure representing a network with optional vertex and edge attributes.

Mirrors R's network class from the statnet network package, while implementing
Graphs.jl's AbstractGraph interface for interoperability with the Julia graph ecosystem.

# Type Parameters
- `T`: Vertex ID type (typically Int)
- `D::Bool`: Compile-time directedness. `Network{T,true}` is directed,
  `Network{T,false}` undirected, so `Graphs.is_directed` is truthful at the
  type level and Graphs.jl trait-dispatched generics use the correct
  (un)directed semantics. In user code you rarely need to spell out `D`:
  `Network{Int}` matches both, and `network(n; directed=...)` /
  `Network{Int}(; directed=...)` pick `D` from the keyword.

# Storage design
The backing store is always a `SimpleDiGraph{T}` (field `graph`), regardless
of `D`. Undirected networks store each edge in *both* arc directions
(self-loops once); all Graphs.jl interface methods (`ne`, `edges`, `degree`,
`has_edge`, …) account for this so each undirected edge is reported exactly
once. This design was chosen over splitting the storage into
`SimpleGraph`/`SimpleDiGraph` because downstream packages (SNA.jl, NDTV.jl, …)
access `net.graph` directly and rely on it being a symmetric `SimpleDiGraph`
for undirected networks; keeping the storage uniform makes the directedness
change purely a type-level (trait) fix with no data-layout fallout.

Directedness is *not* a mutable runtime field: `net.directed` is a read-only
property returning `D`, and attempting `net.directed = x` throws.

# Fields
- `graph::SimpleDiGraph{T}`: Underlying graph storage (symmetric when `D == false`)
- `bipartite::Union{Nothing, Int}`: Number of vertices in first mode if two-mode, nothing otherwise
- `loops::Bool`: Whether self-loops are allowed
- `multiple::Bool`: Whether multiple edges between same vertices are allowed.
  NOT supported by the simple-graph backing store; passing `true` throws.
- `hyper::Bool`: Whether hyperedges are allowed.
  NOT supported; passing `true` throws.
- `vertex_attrs::Dict{Symbol, Dict{T, Any}}`: Vertex attributes by name
- `edge_attrs::Dict{Symbol, Dict{Tuple{T,T}, Any}}`: Edge attributes by vertex pair
- `network_attrs::Dict{Symbol, Any}`: Network-level attributes
- `missing_dyads::Set{Tuple{T,T}}`: Dyads whose tie status is UNOBSERVED
  (statnet-style NA ties). See [`set_missing_dyad!`](@ref) for the semantics;
  manipulate through the missing-dyad API, not directly.
- `n_self_loops::Int`: Number of self-loops stored (internal cache that makes
  `ne` O(1) on undirected networks with loops); maintained by the mutators.

!!! warning "Call Graphs.jl on the `Network`, never on `net.graph`"
    For an undirected network `net.graph` is a `SimpleDiGraph` holding **both**
    arcs of every edge, so a Graphs.jl algorithm run on it treats each tie as
    two directed arcs: `Graphs.core_number(net.graph)` doubles every core
    number, `degree(net.graph)` doubles every degree, betweenness and
    clustering are wrong the same way, and undirected-only algorithms
    (`maximal_cliques`, …) throw. Pass the `Network` itself
    (`Graphs.core_number(net)`), which reports the correct undirected
    semantics, or convert explicitly with `Graphs.SimpleGraph(net)` (an
    undirected copy) when an algorithm needs a concrete Graphs.jl type. Never
    mutate `net.graph` directly either: the mutators keep edge attributes, the
    missing-dyad mask, the two-mode partition and the self-loop count in step
    with it.

# Example
```julia
# Create an empty directed network with 5 vertices
net = Network(5)

# `network(5)` is the equivalent statnet-style alias
net = network(5)

# Create an undirected network
net = network(5; directed=false)

# Add edges
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

# Set vertex attributes
set_vertex_attribute!(net, :name, Dict(1 => "Alice", 2 => "Bob", 3 => "Carol"))

# Set edge attributes
set_edge_attribute!(net, :weight, Dict((1,2) => 1.5, (2,3) => 2.0))
```
"""
mutable struct Network{T<:Integer, D} <: AbstractNetwork{T}
    graph::SimpleDiGraph{T}
    bipartite::Union{Nothing, Int}
    loops::Bool
    multiple::Bool
    hyper::Bool
    vertex_attrs::Dict{Symbol, Dict{T, Any}}
    edge_attrs::Dict{Symbol, Dict{Tuple{T,T}, Any}}
    network_attrs::Dict{Symbol, Any}
    missing_dyads::Set{Tuple{T,T}}
    n_self_loops::Int

    function Network{T, D}(;
        n::Integer = 0,
        bipartite::Union{Nothing, Integer} = nothing,
        loops::Bool = false,
        multiple::Bool = false,
        hyper::Bool = false
    ) where {T<:Integer, D}
        D isa Bool || throw(ArgumentError(
            "the directedness type parameter D of Network{T,D} must be a Bool"))
        # R: "Network objects cannot be of negative order." A negative n used
        # to build an empty 0-vertex network silently.
        n >= 0 || throw(ArgumentError("network: n must be non-negative (got $n)"))
        # The SimpleDiGraph backing store cannot represent parallel or hyper
        # edges; reject the flags rather than silently collapsing the data
        # to a simple graph.
        multiple && throw(ArgumentError(
            "multiple (parallel) edges are not supported by NetworkCore.jl"))
        hyper && throw(ArgumentError(
            "hyperedges are not supported by NetworkCore.jl"))
        !isnothing(bipartite) && !(0 <= bipartite <= n) && throw(ArgumentError(
            "bipartite mode-1 size must be between 0 and n"))
        # Any integer vertex count is accepted (an `Int32` network's
        # `nv(net)::Int32` is what `copy` passes); stored as `Int`.
        n = Int(n)
        bipartite = bipartite === nothing ? nothing : Int(bipartite)
        g = SimpleDiGraph{T}(n)
        new{T, D}(
            g,
            bipartite,
            loops,
            multiple,
            hyper,
            Dict{Symbol, Dict{T, Any}}(),
            Dict{Symbol, Dict{Tuple{T,T}, Any}}(),
            Dict{Symbol, Any}(),
            Set{Tuple{T,T}}(),
            0
        )
    end
end

# Positional form `Network{T,D}(n)`: the constructor Graphs.jl's generic
# algorithms call on a graph type (`T(nv)` in `induced_subgraph`, …).
(::Type{Network{T, D}})(n::Integer) where {T<:Integer, D} = Network{T, D}(; n=n)

"""
    Network{T}(; n=0, directed=true, kwargs...) -> Network{T, directed}

Construct a network specifying only the vertex-ID type; the directedness
type parameter `D` is filled in from the `directed` keyword. This keeps the
familiar one-parameter form (`Network{Int}(; n=5, directed=false)`) working
even though directedness is now part of the type.
"""
function (::Type{Network{T}})(; n::Integer = 0, directed::Bool = true,
                              kwargs...) where {T<:Integer}
    return Network{T, directed}(; n=n, kwargs...)
end

# `net.directed` remains available as a read-only property computed from the
# type parameter D (the field was removed when directedness became a
# compile-time property). Branches on `net.directed` constant-fold.
function Base.getproperty(net::Network{T, D}, name::Symbol) where {T, D}
    name === :directed && return D
    return getfield(net, name)
end

function Base.setproperty!(net::Network, name::Symbol, x)
    name === :directed && throw(ArgumentError(
        "directedness is a type parameter of Network{T,D} and cannot be " *
        "mutated; construct a new network with the desired directedness"))
    val = setfield!(net, name, convert(fieldtype(typeof(net), name), x))
    # Replacing the backing graph wholesale (as `copy` does) re-derives the
    # self-loop count cached for `ne`.
    name === :graph && setfield!(net, :n_self_loops, _count_self_loops(val))
    return val
end

function _count_self_loops(g::SimpleDiGraph)
    c = 0
    for v in Graphs.vertices(g)
        Graphs.has_edge(g, v, v) && (c += 1)
    end
    return c
end

Base.propertynames(net::Network) = (fieldnames(typeof(net))..., :directed)

# Convenience constructor with default Int type
function Network(n::Integer=0; kwargs...)
    Network{Int}(; n=n, kwargs...)
end

"""
    network(n::Int; kwargs...) -> Network

Create a new network with `n` vertices. Alias for Network constructor
to match R's network() function.

# Keyword Arguments
- `directed::Bool=true`: Whether the network is directed
- `bipartite::Union{Nothing,Int}=nothing`: Number of vertices in first mode if bipartite
- `loops::Bool=false`: Whether self-loops are allowed
- `multiple::Bool=false`: Must be `false`; parallel edges are not supported
- `hyper::Bool=false`: Must be `false`; hyperedges are not supported

# Example
```julia
net = network(10; directed=false)
```
"""
network(n::Integer=0; kwargs...) = Network(n; kwargs...)

"""
    network_initialize(n::Int; kwargs...) -> Network

Initialize an empty network with `n` vertices. Alias matching R's
`network.initialize()`; takes the same keywords as [`network`](@ref).

# Example
```julia
net = network_initialize(10; directed=false)
nv(net), ne(net)    # (10, 0)
```
"""
network_initialize(n::Integer; kwargs...) = Network(n; kwargs...)

"""
    BipartiteNetwork{T,D}

A bipartite (two-mode) network with distinct vertex sets. `D` is the
directedness of the wrapped network, carried in the type exactly as in
`Network{T,D}`, so the `network` field is concretely typed; in user code
`BipartiteNetwork{T}` (or plain `BipartiteNetwork`) matches both and the
`directed` keyword picks `D`.

# Fields
- `network::Network{T,D}`: Underlying network
- `n_mode1::Int`: Number of vertices in the first mode
- `n_mode2::Int`: Number of vertices in the second mode

Vertices 1:n_mode1 are in the first mode, (n_mode1+1):(n_mode1+n_mode2) are in the second mode.

`BipartiteNetwork(n₁, n₂; directed=true, loops=false)` wraps
`network(n₁ + n₂; bipartite=n₁)`: the same two-mode semantics (within-mode
edges are rejected by `add_edge!`, `as_matrix` returns the `n₁ × n₂` incidence
matrix, [`network_density`](@ref) takes `discount_bipartite=`), with the mode
sizes carried as fields.

# Example
```julia
b = BipartiteNetwork(2, 3; directed=false)
add_edge!(b, 1, 3)          # mode 1 → mode 2: true
add_edge!(b, 1, 2)          # within mode 1: false
b.n_mode1, b.n_mode2        # (2, 3)
size(as_matrix(b))          # (2, 3)
```
"""
struct BipartiteNetwork{T<:Integer, D} <: AbstractNetwork{T}
    network::Network{T, D}
    n_mode1::Int
    n_mode2::Int

    function BipartiteNetwork{T, D}(n_mode1::Integer, n_mode2::Integer;
                                    kwargs...) where {T<:Integer, D}
        net = Network{T, D}(; n=n_mode1 + n_mode2, bipartite=n_mode1, kwargs...)
        new{T, D}(net, Int(n_mode1), Int(n_mode2))
    end
end

# One-parameter form: directedness from the `directed` keyword, as `Network{T}`.
function (::Type{BipartiteNetwork{T}})(n_mode1::Integer, n_mode2::Integer;
                                       directed::Bool=true, kwargs...) where {T<:Integer}
    return BipartiteNetwork{T, directed}(n_mode1, n_mode2; kwargs...)
end

BipartiteNetwork(n_mode1::Integer, n_mode2::Integer; kwargs...) = BipartiteNetwork{Int}(n_mode1, n_mode2; kwargs...)

# Copy the mutable state (graph and attribute dicts) of `src` into `dst`.
# Outer and inner attribute dicts are duplicated; attribute values themselves
# are shared by reference.
function _copy_state!(dst::Network{T, D}, src::Network{T, D}) where {T, D}
    dst.graph = copy(src.graph)
    dst.vertex_attrs = Dict{Symbol, Dict{T, Any}}(
        attr => copy(vals) for (attr, vals) in src.vertex_attrs)
    dst.edge_attrs = Dict{Symbol, Dict{Tuple{T,T}, Any}}(
        attr => copy(vals) for (attr, vals) in src.edge_attrs)
    dst.network_attrs = copy(src.network_attrs)
    dst.missing_dyads = copy(src.missing_dyads)
    return dst
end

"""
    copy(net::Network) -> Network

Create an independent copy of the network with the same vertices and edges.

The underlying graph, every attribute dictionary (vertex, edge, and
network-level; outer and inner dicts), and the missing-dyad mask (see
[`set_missing_dyad!`](@ref)) are duplicated, so mutating the copy's
structure, attributes, or mask does not affect the original. Attribute
values themselves are shared by reference, not deep-copied. The `directed`,
`bipartite`, and `loops` settings are preserved.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
net2 = copy(net)
add_edge!(net2, 2, 3)   # does not modify net
```
"""
function Base.copy(net::Network{T, D}) where {T, D}
    dup = Network{T, D}(; n=Graphs.nv(net.graph),
                        bipartite=net.bipartite, loops=net.loops)
    return _copy_state!(dup, net)
end

"""
    copy(bnet::BipartiteNetwork) -> BipartiteNetwork

Create an independent copy of the bipartite network. The wrapped `Network` is
copied as described for `copy(::Network)`; mode sizes are preserved.
"""
function Base.copy(bnet::BipartiteNetwork{T, D}) where {T, D}
    dup = BipartiteNetwork{T, D}(bnet.n_mode1, bnet.n_mode2;
                                 loops=bnet.network.loops)
    _copy_state!(dup.network, bnet.network)
    return dup
end

# Helper function to get canonical edge representation
function _canonical_edge(net::Network{T}, i::T, j::T) where T
    if net.directed
        return (i, j)
    else
        # For undirected networks, use canonical ordering
        return minmax(i, j)
    end
end

# Display methods.
#
# Two-arg `show` is the ONE-LINE form Base uses inside containers (a
# `Vector{Network}`, a `Dict{Symbol,Network}`, the `NamedTuple` that
# `load_dataset(:s50)` returns) and for `print`/`string`; the multi-line block
# R users expect from `print.network` is the `MIME"text/plain"` method the
# REPL calls for a top-level value. One multi-line 2-arg method garbled every
# container of networks.
function Base.show(io::IO, net::Network{T, D}) where {T, D}
    print(io, "Network{", T, ",", D, "}(", nv(net), " vertices, ", ne(net), " edges")
    net.loops && print(io, ", loops")
    isnothing(net.bipartite) || print(io, ", two-mode ", net.bipartite, "+", nv(net) - net.bipartite)
    print(io, ")")
end

# The multi-line block. Mirrors what R's `print.network` always states —
# directedness, `loops`, the two-mode size — because a `loops=true` network
# that prints identically to a `loops=false` one, or a bipartite network that
# prints only the word "bipartite", is the classic student confusion.
# Attribute names are sorted, as R lists them.
function _show_network_block(io::IO, net::Network{T}, header::AbstractString) where T
    println(io, header)
    println(io, "  Vertices: ", nv(net))
    println(io, "  Edges: ", ne(net))
    println(io, "  Directed: ", net.directed)
    println(io, "  Loops: ", net.loops)
    if !isnothing(net.bipartite)
        k = net.bipartite
        println(io, "  Two-mode: ", k, " mode-1 vertices (1:", k, "), ",
                nv(net) - k, " mode-2 vertices (", k + 1, ":", nv(net), ")")
    end
    n_vattrs = length(net.vertex_attrs)
    n_eattrs = length(net.edge_attrs)
    n_nattrs = length(net.network_attrs)
    if n_vattrs > 0
        println(io, "  Vertex attributes: ", join(sort!(collect(keys(net.vertex_attrs))), ", "))
    end
    if n_eattrs > 0
        println(io, "  Edge attributes: ", join(sort!(collect(keys(net.edge_attrs))), ", "))
    end
    if n_nattrs > 0
        println(io, "  Network attributes: ", join(sort!(collect(keys(net.network_attrs))), ", "))
    end
    if !isempty(net.missing_dyads)
        println(io, "  Missing (unobserved) dyads: ", length(net.missing_dyads))
    end
end

function Base.show(io::IO, ::MIME"text/plain", net::Network{T}) where T
    dir_str = net.directed ? "directed" : "undirected"
    mode_str = isnothing(net.bipartite) ? "" : " two-mode"
    _show_network_block(io, net, "Network{$T}:$mode_str $dir_str network")
end

function Base.show(io::IO, bnet::BipartiteNetwork{T}) where T
    net = bnet.network
    print(io, "BipartiteNetwork{", T, "}(", bnet.n_mode1, "+", bnet.n_mode2,
          " vertices, ", ne(net), " edges")
    net.loops && print(io, ", loops")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", bnet::BipartiteNetwork{T}) where T
    net = bnet.network
    dir_str = net.directed ? "directed" : "undirected"
    _show_network_block(io, net, "BipartiteNetwork{$T}: $dir_str two-mode network")
end
