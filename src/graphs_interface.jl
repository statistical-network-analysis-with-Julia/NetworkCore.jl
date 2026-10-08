# Graphs.jl AbstractGraph interface implementation for Network types.
#
# This allows Network objects to be used with all Graphs.jl algorithms
# and functions seamlessly.

import Graphs: nv, ne, vertices, edges, has_vertex, has_edge
import Graphs: is_directed, add_vertex!, add_vertices!, add_edge!, rem_vertex!, rem_edge!
import Graphs: neighbors, inneighbors, outneighbors
import Graphs: edgetype, SimpleEdge

# ============================================================================
# Basic Graph Properties
# ============================================================================

"""
    nv(net::Network) -> Int

Return the number of vertices in the network (R's `network.size`).

# Example
```julia
net = network(5)
nv(net)    # 5
```
"""
Graphs.nv(net::Network) = Graphs.nv(net.graph)
Graphs.nv(net::BipartiteNetwork) = Graphs.nv(net.network.graph)

"""
    ne(net::Network) -> Int

Return the number of edges in the network (R's `network.edgecount`).
For undirected networks, each edge is counted once (including self-loops).

# Example
```julia
net = network(3; directed=false)
add_edge!(net, 1, 2)
add_edge!(net, 2, 1)    # false — already present (undirected)
ne(net)                 # 1
```
"""
function Graphs.ne(net::Network)
    n = Graphs.ne(net.graph)
    net.directed && return n
    # Undirected edges are stored in both directions in the backing digraph,
    # except self-loops, which are stored once and must not be halved. The
    # self-loop count is cached on the network (maintained by `add_edge!`,
    # `rem_edge!`, `rem_vertex!` and by assigning `net.graph`), so this is O(1)
    # — it used to count the loops over every vertex on each call.
    nloops = getfield(net, :n_self_loops)
    return (n - nloops) ÷ 2 + nloops
end
Graphs.ne(net::BipartiteNetwork) = Graphs.ne(net.network)

"""
    vertices(net::Network) -> iterator

Return an iterator over all vertices, `1:nv(net)`.

# Example
```julia
net = network(3)
collect(vertices(net))    # [1, 2, 3]
```
"""
Graphs.vertices(net::Network) = Graphs.vertices(net.graph)
Graphs.vertices(net::BipartiteNetwork) = Graphs.vertices(net.network.graph)

"""
    UndirectedEdgeIter{T}

Edge iterator for undirected networks. The backing `SimpleDiGraph` stores
each undirected edge in both arc directions (self-loops once); this iterator
yields each undirected edge exactly once as a `SimpleEdge{T}` with
`src <= dst`, and carries a correct `length` (`ne(net)`) and `eltype` as
required by the Graphs.jl interface.
"""
struct UndirectedEdgeIter{T<:Integer}
    net::Network{T, false}
end

Base.length(it::UndirectedEdgeIter) = Graphs.ne(it.net)
Base.eltype(::Type{UndirectedEdgeIter{T}}) where T = SimpleEdge{T}
Base.IteratorSize(::Type{UndirectedEdgeIter{T}}) where T = Base.HasLength()
Base.IteratorEltype(::Type{UndirectedEdgeIter{T}}) where T = Base.HasEltype()

function Base.iterate(it::UndirectedEdgeIter, state=nothing)
    inner = Graphs.edges(getfield(it.net, :graph))
    next = state === nothing ? iterate(inner) : iterate(inner, state)
    while next !== nothing
        e, s = next
        src(e) <= dst(e) && return (e, s)
        next = iterate(inner, s)
    end
    return nothing
end

"""
    edges(net::Network) -> iterator

Return an iterator over all edges. For undirected networks each edge is
yielded exactly once (canonical `src <= dst` orientation) even though the
backing digraph stores both arc directions; the iterator has a correct
`length` and `eltype`.

# Example
```julia
net = network(3; directed=false)
add_edges!(net, [(2, 1), (2, 3)])
[(src(e), dst(e)) for e in edges(net)]    # [(1, 2), (2, 3)]
length(edges(net))                        # 2
```
"""
Graphs.edges(net::Network{T, true}) where T = Graphs.edges(net.graph)
Graphs.edges(net::Network{T, false}) where T = UndirectedEdgeIter{T}(net)
Graphs.edges(net::BipartiteNetwork) = Graphs.edges(net.network)

@doc """
    src(e) -> Int

Source (first endpoint) of an edge yielded by [`edges`](@ref) — the Graphs.jl
accessor, re-exported with the iterator so that `for e in edges(net);
src(e) …` works after a plain `using NetworkCore` (or `using ERGM`). On an
undirected network the edge is yielded in canonical orientation, so
`src(e) <= dst(e)`.

# Example
```julia
net = network(3; directed=true)
add_edge!(net, 3, 1)
e = first(edges(net))
src(e), dst(e)      # (3, 1)
```
""" Graphs.src

@doc """
    dst(e) -> Int

Destination (second endpoint) of an edge yielded by [`edges`](@ref) — the
Graphs.jl accessor, re-exported with the iterator; see [`src`](@ref).

# Example
```julia
net = network(3; directed=false)
add_edge!(net, 3, 1)
dst(first(edges(net)))    # 3 — canonical orientation (1, 3)
```
""" Graphs.dst

"""
    has_vertex(net::Network, v) -> Bool

Check if vertex `v` exists in the network, i.e. `1 <= v <= nv(net)`.

# Example
```julia
net = network(3)
has_vertex(net, 3)    # true
has_vertex(net, 4)    # false
```
"""
Graphs.has_vertex(net::Network, v) = Graphs.has_vertex(net.graph, v)
Graphs.has_vertex(net::BipartiteNetwork, v) = Graphs.has_vertex(net.network.graph, v)

"""
    has_edge(net::Network, i, j) -> Bool

Check if edge `(i, j)` exists in the network (R's `is.adjacent`). Symmetric
for undirected networks. Reads the face value: a masked (unobserved) dyad
reports its stored tie status — see [`is_missing_dyad`](@ref).

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
has_edge(net, 1, 2)    # true
has_edge(net, 2, 1)    # false — directed
```
"""
function Graphs.has_edge(net::Network, i, j)
    if net.directed
        return Graphs.has_edge(net.graph, i, j)
    else
        return Graphs.has_edge(net.graph, i, j) || Graphs.has_edge(net.graph, j, i)
    end
end
Graphs.has_edge(net::BipartiteNetwork, i, j) = Graphs.has_edge(net.network, i, j)
Graphs.add_edge!(net::BipartiteNetwork, i, j) = Graphs.add_edge!(net.network, i, j)
Graphs.rem_edge!(net::BipartiteNetwork, i, j) = Graphs.rem_edge!(net.network, i, j)

"""
    is_directed(net::Network) -> Bool

Check if the network is directed.

Directedness is a type parameter of `Network{T,D}`, so the type-level trait
`Graphs.is_directed(typeof(net))` is truthful and Graphs.jl generics that
dispatch on it use the correct (un)directed semantics.

# Example
```julia
is_directed(network(3))                    # true (the default)
is_directed(network(3; directed=false))    # false
```
"""
Graphs.is_directed(::Network{T, D}) where {T, D} = D
Graphs.is_directed(::Type{Network{T, D}}) where {T, D} = D
Graphs.is_directed(::BipartiteNetwork{T, D}) where {T, D} = D
Graphs.is_directed(::Type{BipartiteNetwork{T, D}}) where {T, D} = D

"""
    is_bipartite(net::Network) -> Bool

Check whether the network is bipartite in the graph-theoretic sense, i.e.
whether its vertices can be two-colored so that every edge joins the two
color classes (computed by Graphs.jl's algorithm on the underlying graph).

This is *not* the two-mode metadata flag set via `network(n; bipartite=k)`;
use [`is_two_mode`](@ref) for that.
"""
Graphs.is_bipartite(net::Network) = Graphs.is_bipartite(net.graph)
# BipartiteNetwork admits cross-mode edges only, so it is bipartite by
# construction.
Graphs.is_bipartite(net::BipartiteNetwork) = true

"""
    is_two_mode(net) -> Bool

Return whether the network is flagged as a two-mode (bipartite) network,
i.e. whether it carries the statnet-style `bipartite` metadata giving the
size of the first vertex mode (`network(n; bipartite=k)` or
`BipartiteNetwork`). This is pure metadata about how the vertex set is
partitioned; for the graph-theoretic property use `Graphs.is_bipartite`.

# Example
```julia
is_two_mode(network(5; bipartite=2))    # true
is_two_mode(BipartiteNetwork(2, 3))     # true
is_two_mode(network(5))                 # false
```
"""
is_two_mode(net::Network) = !isnothing(net.bipartite)
is_two_mode(net::BipartiteNetwork) = true

# ============================================================================
# Edge Type
# ============================================================================

Graphs.edgetype(net::Network{T}) where T = SimpleEdge{T}
Graphs.edgetype(net::BipartiteNetwork{T}) where T = SimpleEdge{T}

# ============================================================================
# Zero Graph (required by the Graphs.jl interface)
# ============================================================================

"""
    zero(::Type{Network{T,D}}) -> Network{T,D}

Return an empty network of the given type, as required by the Graphs.jl
interface. Directedness is encoded in the type parameter `D` and is
preserved. For the one-parameter form `zero(Network{T})` (directedness
unspecified) an empty *directed* network is returned, matching the
constructor default; use `zero(net)` on an instance to match its
directedness.
"""
Base.zero(::Type{Network{T, D}}) where {T, D} = Network{T, D}()
Base.zero(::Type{Network{T}}) where T = Network{T, true}()
Base.zero(::Type{Network}) = Network{Int, true}()

"""
    zero(net::Network) -> Network

Return an empty network with the same eltype, directedness, and `loops`
setting as `net`.
"""
Base.zero(net::Network{T, D}) where {T, D} = Network{T, D}(; loops=net.loops)

# ============================================================================
# Vertex Modification
# ============================================================================

"""
    add_vertex!(net::Network) -> Bool

Add a new vertex to the network (as vertex `nv(net) + 1`). Returns true if
successful.

# Example
```julia
net = network(3)
add_vertex!(net)    # true
nv(net)             # 4
```
"""
function Graphs.add_vertex!(net::Network)
    return Graphs.add_vertex!(net.graph)
end

"""
    add_vertices!(net::Network, n::Int) -> Int

Add n new vertices to the network. Returns the number of vertices added.

# Example
```julia
net = network(3)
add_vertices!(net, 2)    # 2
nv(net)                  # 5
```
"""
function Graphs.add_vertices!(net::Network, n::Int)
    added = 0
    for _ in 1:n
        if Graphs.add_vertex!(net.graph)
            added += 1
        end
    end
    return added
end

"""
    rem_vertex!(net::Network, v) -> Bool

Remove vertex v from the network. Returns true if successful.

All edges incident to `v` and their attributes are removed. On a **one-mode**
network, following the Graphs.jl swap-and-pop convention, the last vertex is
renumbered to `v`; its vertex attributes and the attributes of its incident
edges are remapped to the new ID, so attribute data stays attached to the
right vertices.

On a **two-mode** network swap-and-pop would move the last (mode-2) vertex into
a mode-1 slot and turn its edges into within-mode pairs, so the vertices are
instead renumbered **stably** — every vertex after `v` moves down by one — and
the mode-1 size `bipartite` drops by one when `v` is a mode-1 vertex. This is
R's `delete.vertices`; attributes and missing-dyad marks move with their
vertices. Use `Graphs.rem_vertices!(net, vs)` to delete several vertices at
once with the same stable renumbering.

# Example
```julia
net = network(3)
set_vertex_attribute!(net, :name, ["a", "b", "c"])
rem_vertex!(net, 1)                     # true; vertex 3 ("c") is renumbered to 1
nv(net)                                 # 2
get_vertex_attribute(net, :name, 1)     # "c"

b = network(5; bipartite=2, directed=false)
add_edges!(b, [(1, 3), (2, 5)])
rem_vertex!(b, 1)                       # stable: old 2..5 become 1..4
b.bipartite                             # 1
has_edge(b, 1, 4)                       # true — old edge (2, 5)
```
"""
function Graphs.rem_vertex!(net::Network{T}, v) where T
    if !has_vertex(net, v)
        return false
    end
    if is_two_mode(net)
        _rem_vertices_stable!(net, (v,))
        return true
    end

    vT = T(v)
    last = T(nv(net.graph))

    # Remove edges involving this vertex from edge attributes
    for attrs in values(net.edge_attrs)
        for edge in collect(keys(attrs))
            if edge[1] == vT || edge[2] == vT
                delete!(attrs, edge)
            end
        end
    end

    # Remove vertex from vertex attributes
    for attrs in values(net.vertex_attrs)
        delete!(attrs, vT)
    end

    # Drop missing-dyad marks involving this vertex
    for dyad in collect(net.missing_dyads)
        if dyad[1] == vT || dyad[2] == vT
            delete!(net.missing_dyads, dyad)
        end
    end

    # Graphs.rem_vertex! swaps the last vertex into slot v (swap-and-pop),
    # so remap the last vertex's attributes onto the new ID v
    if vT != last
        for attrs in values(net.vertex_attrs)
            if haskey(attrs, last)
                attrs[vT] = attrs[last]
                delete!(attrs, last)
            end
        end
        for attrs in values(net.edge_attrs)
            for edge in collect(keys(attrs))
                if edge[1] == last || edge[2] == last
                    val = attrs[edge]
                    delete!(attrs, edge)
                    a = edge[1] == last ? vT : edge[1]
                    b = edge[2] == last ? vT : edge[2]
                    attrs[_canonical_edge(net, a, b)] = val
                end
            end
        end
        # Remap missing-dyad marks of the renumbered last vertex
        for dyad in collect(net.missing_dyads)
            if dyad[1] == last || dyad[2] == last
                delete!(net.missing_dyads, dyad)
                a = dyad[1] == last ? vT : dyad[1]
                b = dyad[2] == last ? vT : dyad[2]
                push!(net.missing_dyads, _canonical_edge(net, a, b))
            end
        end
    end

    # Remove from graph (this also removes incident edges); a self-loop on
    # `v` leaves the cached count with it.
    had_loop = Graphs.has_edge(net.graph, v, v)
    result = Graphs.rem_vertex!(net.graph, v)
    (result && had_loop) && setfield!(net, :n_self_loops, getfield(net, :n_self_loops) - 1)
    return result
end

"""
    rem_vertices!(net::Network, vs) -> Vector

Remove the vertices `vs` (any iterable of vertex IDs; repeats are ignored)
from the network, with R's `delete.vertices` semantics: the remaining vertices
keep their relative order and are renumbered `1, 2, …` (stable renumbering,
never swap-and-pop), every edge incident to a removed vertex goes with its
attributes and missing-dyad marks, and on a two-mode network `bipartite` drops
by the number of removed mode-1 vertices. Returns `vmap`, the original ID of
each remaining vertex (`vmap[new] == old`), as Graphs.jl's `rem_vertices!`
does. A vertex outside `1:nv(net)` is an `ArgumentError` and leaves the
network untouched.

This is a method of Graphs.jl's `rem_vertices!` (not exported by NetworkCore.jl:
call it as `Graphs.rem_vertices!` or after `using Graphs`).

# Example
```julia
using Graphs: rem_vertices!
net = network(5; bipartite=2, directed=false)
add_edges!(net, [(1, 3), (2, 4), (2, 5)])
set_vertex_attribute!(net, :name, ["a", "b", "c", "d", "e"])
vmap = rem_vertices!(net, [1, 4])        # [2, 3, 5]
net.bipartite                            # 1
vertex_attribute_vector(net, :name, String)   # ["b", "c", "e"]
has_edge(net, 1, 3)                      # true — old edge (2, 5)
```
"""
function Graphs.rem_vertices!(net::Network{T}, vs::AbstractVector{<:Integer}) where T
    n = nv(net)
    for v in vs
        1 <= v <= n || throw(ArgumentError(
            "rem_vertices!: vertex $v is not a vertex of the network (it has $n vertices)"))
    end
    return _rem_vertices_stable!(net, vs)
end

# Stable (order-preserving) vertex deletion — R's `delete.vertices`. Rebuilds
# the backing graph and every vertex-keyed store through one old → new map, so
# the two-mode blocks stay contiguous and `bipartite` stays truthful.
function _rem_vertices_stable!(net::Network{T}, vs) where T
    n = nv(net)
    remove = falses(n)
    for v in vs
        remove[v] = true
    end
    newid = zeros(T, n)
    vmap = T[]
    for u in 1:n
        remove[u] && continue
        push!(vmap, T(u))
        newid[u] = T(length(vmap))
    end
    kept(i) = !remove[i]

    g = SimpleDiGraph{T}(length(vmap))
    for e in Graphs.edges(net.graph)
        a, b = Graphs.src(e), Graphs.dst(e)
        (kept(a) && kept(b)) && Graphs.add_edge!(g, Graphs.SimpleEdge{T}(newid[a], newid[b]))
    end

    for attr in collect(keys(net.vertex_attrs))
        vals = net.vertex_attrs[attr]
        net.vertex_attrs[attr] = Dict{T, Any}(newid[k] => x for (k, x) in vals if kept(k))
    end
    for attr in collect(keys(net.edge_attrs))
        vals = net.edge_attrs[attr]
        net.edge_attrs[attr] = Dict{Tuple{T,T}, Any}(
            _canonical_edge(net, newid[i], newid[j]) => x
            for ((i, j), x) in vals if kept(i) && kept(j))
    end
    net.missing_dyads = Set{Tuple{T,T}}(
        _canonical_edge(net, newid[i], newid[j])
        for (i, j) in net.missing_dyads if kept(i) && kept(j))
    if !isnothing(net.bipartite)
        k = net.bipartite
        net.bipartite = k - count(view(remove, 1:k))
    end
    net.graph = g          # re-derives the cached self-loop count
    return vmap
end

# ============================================================================
# Edge Modification
# ============================================================================

"""
    add_edge!(net::Network, i, j) -> Bool

Add edge (i, j) to the network. Returns true if successful.
For undirected networks, also adds (j, i).
For bipartite networks, edges must connect the two modes; within-mode
edges are rejected (returns false).

The other half of the sampler toggle (see [`rem_edge!`](@ref), which is 0 B).
Adds no allocation over Graphs.jl's own `add_edge!` on the wrapped
`SimpleDiGraph`: the edge is handed over as a concrete `SimpleEdge{T}` rather
than through the generic `add_edge!(g, x, y)`, whose `edgetype(g)(x, y)`
construction allocates, and the attribute store is skipped when empty. The
wrapped graph itself may allocate a few dozen bytes on alternate toggles after
a removal (a Base `insert!` regrowing the adjacency list after `deleteat!`),
and this method allocates exactly as much — the "allocation-free" testset pins
equality with the bare `SimpleDiGraph`, not zero.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)    # true
add_edge!(net, 1, 2)    # false — already present
add_edge!(net, 1, 1)    # false — loops=false
has_edge(net, 1, 2)     # true
```
"""
function Graphs.add_edge!(net::Network{T}, i, j) where T
    # Check self-loop constraint
    if i == j && !net.loops
        return false
    end

    # Bipartite networks only admit cross-mode edges
    if !isnothing(net.bipartite)
        k = net.bipartite
        ((i <= k) ⊻ (j <= k)) || return false
    end

    result = Graphs.add_edge!(net.graph, Graphs.SimpleEdge{T}(T(i), T(j)))

    # For undirected networks, add reverse edge too
    if !net.directed && result
        Graphs.add_edge!(net.graph, Graphs.SimpleEdge{T}(T(j), T(i)))
    end
    # Keep the self-loop count `ne` reads in step with the graph.
    (result && i == j) && setfield!(net, :n_self_loops, getfield(net, :n_self_loops) + 1)

    return result
end

"""
    add_edge!(net::Network, i, j, attrs::Dict) -> Bool

Add edge `(i, j)` and, if it was added, store each `attr => value` pair of
`attrs` on it. Returns whether the edge was added.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2, Dict(:weight => 2.5))    # true
get_edge_attribute(net, :weight, 1, 2)        # 2.5
```
"""
function add_edge!(net::Network{T}, i, j, attrs::Dict{Symbol}) where T
    result = Graphs.add_edge!(net, i, j)
    if result
        edge = _canonical_edge(net, T(i), T(j))
        for (attr, val) in attrs
            if !haskey(net.edge_attrs, attr)
                net.edge_attrs[attr] = Dict{Tuple{T,T}, Any}()
            end
            net.edge_attrs[attr][edge] = val
        end
    end
    return result
end

"""
    add_edges!(net::Network, edges) -> Int

Add multiple edges from an iterable of `(i, j)` pairs (R's `add.edges`).
Returns the number of edges newly added — an edge already present counts as
zero and is not an error.

**Every listed edge ends up in the network or the call throws.** An endpoint
outside `1:nv(net)` (an undersized `network(n)`, a mistyped ID), a self-loop
when `loops=false`, and a within-mode pair on a two-mode network are
`ArgumentError`s naming the offending edge, as R's `add.edges` errors on an
illegal vertex reference — never a silently smaller network. (The single
[`add_edge!`](@ref) keeps Graphs.jl's convention of returning `false`.)

# Example
```julia
net = network(4)
add_edges!(net, [(1, 2), (2, 3), (1, 2)])    # 2 — the repeat is not re-added
ne(net)                                      # 2
try add_edges!(net, [(3, 4), (1, 99)]) catch err; err end   # ArgumentError: edge (1, 99) references a vertex outside 1:4
```
"""
function add_edges!(net::Network, edge_list)
    n = nv(net)
    added = 0
    for (i, j) in edge_list
        (1 <= i <= n && 1 <= j <= n) ||
            throw(ArgumentError("add_edges!: edge ($i, $j) references a vertex outside " *
                                "1:$n; pass a larger n= to network() / add_vertices! " *
                                "first, or fix the edge list"))
        (net.loops || i != j) ||
            throw(ArgumentError("add_edges!: edge ($i, $j) is a self-loop but the " *
                                "network has loops=false; construct it with " *
                                "loops=true or drop the edge"))
        if !isnothing(net.bipartite)
            k = net.bipartite
            ((i <= k) ⊻ (j <= k)) ||
                throw(ArgumentError("add_edges!: edge ($i, $j) joins two vertices of " *
                                    "the same mode of a two-mode network (mode 1 is " *
                                    "1:$k, mode 2 is $(k + 1):$n); only cross-mode " *
                                    "edges are allowed"))
        end
        if Graphs.add_edge!(net, i, j)
            added += 1
        end
    end
    return added
end

"""
    rem_edge!(net::Network, i, j) -> Bool

Remove edge (i, j) from the network, together with every edge attribute stored
on it. Returns `true` if the edge existed and was removed, `false` (and touches
nothing) if it did not.

This is the toggle half of every Metropolis–Hastings sampler in the ecosystem,
so it is allocation-free: on a network with no edge attributes it never looks at
the attribute store, and on one with attributes it walks the attribute
dictionaries by value (one `delete!` per attribute, no key lookup). Pinned by
the "rem_edge!/add_edge! are allocation-free" testset.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 2.5)
rem_edge!(net, 1, 2)        # true
has_edge(net, 1, 2)         # false
get_edge_attribute(net, :weight, 1, 2) === nothing   # attribute went with it
rem_edge!(net, 1, 2)        # false — already absent
```
"""
function Graphs.rem_edge!(net::Network{T}, i, j) where T
    # An absent edge is a no-op: nothing to remove, no attributes to drop.
    Graphs.has_edge(net, i, j) || return false

    # Drop edge attributes. Iterate the attribute dictionaries by VALUE (no
    # `keys` + `getindex` round trip) and skip the store entirely when the
    # network carries no edge attributes — the common case inside a sampler.
    if !isempty(net.edge_attrs)
        edge = _canonical_edge(net, T(i), T(j))
        for attr in values(net.edge_attrs)
            delete!(attr, edge)
        end
    end

    result = Graphs.rem_edge!(net.graph, Graphs.SimpleEdge{T}(T(i), T(j)))

    # For undirected, also remove reverse edge
    if !net.directed
        Graphs.rem_edge!(net.graph, Graphs.SimpleEdge{T}(T(j), T(i)))
    end
    (result && i == j) && setfield!(net, :n_self_loops, getfield(net, :n_self_loops) - 1)

    return result
end

# ============================================================================
# Neighbors
# ============================================================================

"""
    neighbors(net::Network, v) -> iterator

Return neighbors of vertex v.
For directed networks, returns outgoing neighbors.

# Example
```julia
net = network(4)
add_edges!(net, [(1, 2), (1, 3), (4, 1)])
neighbors(net, 1)    # [2, 3]
```
"""
function Graphs.neighbors(net::Network, v::Integer)
    if net.directed
        return Graphs.outneighbors(net.graph, v)
    else
        # For undirected, outneighbors gives all neighbors since we store both directions
        return Graphs.outneighbors(net.graph, v)
    end
end

"""
    inneighbors(net::Network, v) -> iterator

Return incoming neighbors of vertex v (all neighbors on an undirected network).

# Example
```julia
net = network(4)
add_edges!(net, [(1, 2), (1, 3), (4, 1)])
inneighbors(net, 1)    # [4]
```
"""
Graphs.inneighbors(net::Network, v) = Graphs.inneighbors(net.graph, v)
Graphs.neighbors(net::BipartiteNetwork, v::Integer) = Graphs.neighbors(net.network, v)
Graphs.inneighbors(net::BipartiteNetwork, v) = Graphs.inneighbors(net.network, v)
Graphs.outneighbors(net::BipartiteNetwork, v) = Graphs.outneighbors(net.network, v)

"""
    outneighbors(net::Network, v) -> iterator

Return outgoing neighbors of vertex v (all neighbors on an undirected network).

# Example
```julia
net = network(4)
add_edges!(net, [(1, 2), (1, 3), (4, 1)])
outneighbors(net, 1)    # [2, 3]
```
"""
Graphs.outneighbors(net::Network, v) = Graphs.outneighbors(net.graph, v)

# ============================================================================
# Degree (Graphs.jl bindings, re-exported)
# ============================================================================
#
# `degree`, `indegree` and `outdegree` are Graphs.jl's generics, which work on
# a `Network` through the `AbstractGraph` interface above; NetworkCore.jl adds no
# methods. They are re-exported so that a plain
# `using NetworkCore` answers the first descriptive question a user asks, and the
# docstrings below — attached to the Graphs.jl bindings from this module — say
# what they mean on a `Network` specifically.

@doc """
    degree(net::Network, v) -> Int
    degree(net::Network) -> Vector{Int}

Degree of vertex `v` (or of every vertex), the Graphs.jl generic re-exported by
NetworkCore.jl. On a directed network `degree(net, v) == indegree(net, v) +
outdegree(net, v)`; on an undirected network all three agree. A self-loop
counts once on an undirected network and twice (one in, one out) on a
directed one, as in Graphs.jl.

Reads face values: a masked (unobserved) dyad contributes its stored tie
status. For a degree that honours the missing-data contract use
`SNA.degree_centrality(net; missing=...)`.

# Example
```julia
net = Network(3)
add_edge!(net, 1, 2)
add_edge!(net, 3, 1)
degree(net, 1)      # 2
degree(net)         # [2, 1, 1]
```
""" Graphs.degree

@doc """
    indegree(net::Network, v) -> Int
    indegree(net::Network) -> Vector{Int}

Number of edges pointing into `v` (Graphs.jl generic, re-exported). Equal to
[`degree`](@ref) on an undirected network. Reads face values at masked dyads.

# Example
```julia
net = Network(3)
add_edge!(net, 1, 2)
indegree(net, 2)    # 1
indegree(net, 1)    # 0
```
""" Graphs.indegree

@doc """
    outdegree(net::Network, v) -> Int
    outdegree(net::Network) -> Vector{Int}

Number of edges leaving `v` (Graphs.jl generic, re-exported). Equal to
[`degree`](@ref) on an undirected network. Reads face values at masked dyads.

# Example
```julia
net = Network(3)
add_edge!(net, 1, 2)
outdegree(net, 1)   # 1
outdegree(net, 2)   # 0
```
""" Graphs.outdegree

# ============================================================================
# Additional Utility Functions
# ============================================================================

"""
    network_size(net::Network) -> Int

Return the number of vertices. Alias for nv() to match R's network.size().

# Example
```julia
network_size(network(7))    # 7
```
"""
network_size(net::AbstractNetwork) = nv(net)

"""
    network_edgecount(net::Network) -> Int

Return the number of edges. Alias for ne() to match R's network.edgecount().

# Example
```julia
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
network_edgecount(net)    # 2
```
"""
network_edgecount(net::AbstractNetwork) = ne(net)

"""
    network_density(net; missing=:error, diag=net.loops, discount_bipartite=false) -> Float64

Network density: the number of edges divided by the number of dyads that could
carry one, following R's `network.density`. Returns `NaN` for networks with
fewer than two vertices (R does the same).

The denominator is chosen from the network's directedness and the `diag` and
`discount_bipartite` keywords:

| Network | `diag=false` | `diag=true` |
|---|---|---|
| Directed | `n (n - 1)` | `n²` |
| Undirected | `n (n - 1) / 2` | `n (n + 1) / 2` |
| Two-mode, `discount_bipartite=true`, undirected | `n₁ n₂` | `n₁ n₂` |
| Two-mode, `discount_bipartite=true`, directed | `2 n₁ n₂` | `2 n₁ n₂` |

# Keywords

- `missing::Symbol=:error` — the missing-dyad policy (see
  [`require_observed`](@ref)). By default a network with masked (unobserved)
  dyads is **refused**; `missing=:face` counts every edge at its stored face
  value, over a denominator that still contains every dyad. This is R's
  `network.density(x, na.omit=FALSE)`. R's default `na.omit=TRUE` — drop the
  NA edges from the numerator but keep the denominator — is deliberately **not**
  offered: it reads an unobserved tie as *absent*, which is the silent
  data-invention the missing-data contract exists to stop.
- `diag::Bool=net.loops` — whether self-loops count as possible dyads.
  Defaults to the network's `loops` flag, exactly as R keys the denominator on
  `has.loops(x)`, so an existing call gives the same value it always did.
  `diag=false` on a loops-allowed network mirrors `sna::gden(diag=FALSE)`:
  self-loops are dropped from **both** the numerator and the denominator.
- `discount_bipartite::Bool=false` — see below.

# Two-mode networks

A two-mode network (`network(n; bipartite=n₁)` or a [`BipartiteNetwork`](@ref))
is measured by default over **all** `n (n - 1)` (or `n (n - 1) / 2`) dyads,
including the within-mode dyads that can never carry a tie — this is R's
default, `discount.bipartite=FALSE`. Pass `discount_bipartite=true` to use only
the `n₁ n₂` cross-mode dyads (`2 n₁ n₂` for a directed two-mode network, one per
direction), which is the "density of the incidence matrix" most two-mode
analyses mean. The keyword is a no-op on a one-mode network. Both the flagged
`Network` and the `BipartiteNetwork` wrapper give the same answer for the same
edge set and keyword.

# Example
```julia
net = network(4; directed=false)
add_edges!(net, [(1, 2), (2, 3), (3, 4)])
network_density(net)                              # 0.5 == 3 / 6

b = network(5; bipartite=2, directed=false)
add_edges!(b, [(1, 3), (1, 4), (2, 5)])
network_density(b)                                # 0.3 — over all 10 dyads (R default)
network_density(b; discount_bipartite=true)       # 0.5 — over the 2·3 cross-mode dyads

l = network(3; loops=true)
add_edges!(l, [(1, 1), (1, 2)])
network_density(l)                                # 2/9   — loops count (diag=true)
network_density(l; diag=false)                    # 1/6   — self-loop dropped, sna::gden(diag=FALSE)

set_missing_dyad!(net, 1, 4)
try network_density(net) catch err; err end       # ArgumentError: masked dyad
network_density(net; missing=:face)               # 0.5 — face value, in writing
```

See also [`missing_policies`](@ref) (this routine declares `(:error, :face)`),
[`is_two_mode`](@ref).
"""
function network_density(net::Network; missing::Symbol=:error, diag::Bool=net.loops,
                         discount_bipartite::Bool=false)
    policy = missing   # the keyword, not `Base.missing`
    require_observed(net, policy; context="network_density")

    n = nv(net)
    n <= 1 && return NaN

    m = ne(net)
    # `diag=false` on a loops-allowed network: self-loops are not dyads, so
    # they leave the numerator as well as the denominator (sna::gden(diag=)).
    if !diag && net.loops
        m -= _n_self_loops(net)
    end

    if discount_bipartite && is_two_mode(net)
        n1 = net.bipartite::Int
        n2 = n - n1
        max_edges = net.directed ? 2 * n1 * n2 : n1 * n2
    elseif net.directed
        max_edges = diag ? n * n : n * (n - 1)
    else
        max_edges = diag ? n * (n + 1) ÷ 2 : n * (n - 1) ÷ 2
    end

    return m / max_edges
end

# Explicit keywords, not a `kwargs...` splat (the splat allocated 16–32 B per
# call). The wrapper's `network` field is a concrete `Network{T,D}`, so the
# forward is statically dispatched; the testset pins `@allocated == 0` for
# both methods. `diag` defaults to the wrapped network's `loops` flag exactly
# as the `Network` method does.
function network_density(bnet::BipartiteNetwork; missing::Symbol=:error,
                         diag::Union{Nothing,Bool}=nothing,
                         discount_bipartite::Bool=false)
    net = bnet.network
    return network_density(net; missing, diag=something(diag, net.loops), discount_bipartite)
end

# `network_density` is one of the routines that CAN read face values on
# request, so it declares the two policies its `missing=` keyword takes.
missing_policies(::typeof(network_density)) = (:error, :face)

# Number of self-loops stored in the network (each counted once): the cache
# the mutators maintain.
_n_self_loops(net::Network) = getfield(net, :n_self_loops)

"""
    get_neighborhood(net::Network, v, order::Int=1) -> Set

Get the neighborhood of vertex `v` up to the specified order (distance),
following out-neighbors on a directed network: the set of vertices reachable
from `v` in at most `order` steps, **including `v` itself** (R's
`get.neighborhood` returns the order-1 neighbors without `v`).

# Example
```julia
net = network(5)
add_edges!(net, [(1, 2), (2, 3), (3, 4)])
get_neighborhood(net, 1)       # Set([1, 2])
get_neighborhood(net, 1, 2)    # Set([1, 2, 3])
```
"""
function get_neighborhood(net::Network{T}, v, order::Int=1) where T
    neighborhood = Set{T}([T(v)])
    frontier = Set{T}([T(v)])

    for _ in 1:order
        new_frontier = Set{T}()
        for u in frontier
            for w in neighbors(net, u)
                if !(w in neighborhood)
                    push!(new_frontier, w)
                    push!(neighborhood, w)
                end
            end
        end
        frontier = new_frontier
        if isempty(frontier)
            break
        end
    end

    return neighborhood
end

"""
    get_induced_subgraph(net::Network, vlist) -> Network

Create a new network containing only the specified vertices and the edges
between them (R's `get.inducedSubgraph`). Vertex `vlist[i]` becomes vertex `i`
of the result; vertex, edge and network attributes and missing-dyad marks
between retained vertices are carried over. `vlist` must be distinct vertices
of `net` (`ArgumentError` naming the offending entry otherwise, as R's
"Illegal vertex selection").

On a **two-mode** network the result is two-mode too, with `bipartite` equal to
the number of selected mode-1 vertices (R's `get.inducedSubgraph` keeps the
partition the same way). Because vertex `vlist[i]` becomes vertex `i`, the
selected mode-1 vertices must come before the selected mode-2 ones in `vlist`
— otherwise the result could not be a contiguous two-mode network, and an
`ArgumentError` says so (`sort(vlist)` always qualifies) instead of silently
returning a one-mode network.

# Example
```julia
net = network(5)
add_edges!(net, [(1, 2), (2, 3), (3, 4)])
sub = get_induced_subgraph(net, [2, 3, 5])
nv(sub), ne(sub)         # (3, 1) — only 2→3 survives, as 1→2
has_edge(sub, 1, 2)      # true
```
"""
function get_induced_subgraph(net::Network{T}, vlist) where T
    for v in vlist
        1 <= v <= nv(net) ||
            throw(ArgumentError("get_induced_subgraph: vertex $v is not a vertex of " *
                                "the network (it has $(nv(net)) vertices)"))
    end
    allunique(vlist) ||
        throw(ArgumentError("get_induced_subgraph: vlist contains a repeated vertex " *
                            "($(first(v for v in vlist if count(==(v), vlist) > 1))); " *
                            "each vertex may be selected once"))
    vset = Set(vlist)
    n = length(vlist)

    # A two-mode selection stays two-mode when its mode-1 vertices come first.
    sub_bipartite = nothing
    if is_two_mode(net)
        k = net.bipartite
        n1 = count(<=(k), vlist)
        all(i -> (vlist[i] <= k) == (i <= n1), eachindex(vlist)) ||
            throw(ArgumentError("get_induced_subgraph: on a two-mode network the " *
                                "selected mode-1 vertices (1:$k) must precede the " *
                                "selected mode-2 vertices in vlist, so that the " *
                                "subgraph is a two-mode network with bipartite = $n1; " *
                                "pass sort(vlist)"))
        sub_bipartite = n1
    end

    # Create mapping from old to new vertex IDs
    old_to_new = Dict{T, T}(T(v) => T(i) for (i, v) in enumerate(vlist))

    # Create new network
    sub = Network{T}(; n=n, directed=net.directed, loops=net.loops,
                     bipartite=sub_bipartite)

    # Add edges (each is a legal edge of `sub` by construction; never drop one)
    for e in edges(net)
        i, j = src(e), dst(e)
        if i in vset && j in vset
            add_edge!(sub, old_to_new[i], old_to_new[j]) ||
                error("get_induced_subgraph: internal error, edge ($i, $j) could not be carried over")
        end
    end

    # Copy vertex attributes
    for (attr, vals) in net.vertex_attrs
        new_vals = Dict{T, Any}()
        for (v, val) in vals
            if v in vset
                new_vals[old_to_new[v]] = val
            end
        end
        if !isempty(new_vals)
            sub.vertex_attrs[attr] = new_vals
        end
    end

    # Copy edge attributes
    for (attr, vals) in net.edge_attrs
        new_vals = Dict{Tuple{T,T}, Any}()
        for ((i, j), val) in vals
            if i in vset && j in vset
                new_edge = _canonical_edge(sub, old_to_new[i], old_to_new[j])
                new_vals[new_edge] = val
            end
        end
        if !isempty(new_vals)
            sub.edge_attrs[attr] = new_vals
        end
    end

    # Copy network attributes
    for (attr, val) in net.network_attrs
        sub.network_attrs[attr] = val
    end

    # Carry over missing-dyad marks between retained vertices
    for (i, j) in net.missing_dyads
        if i in vset && j in vset
            push!(sub.missing_dyads, _canonical_edge(sub, old_to_new[i], old_to_new[j]))
        end
    end

    return sub
end

"""
    permute_vertices(net::Network, perm::Vector) -> Network

Create a new network with vertices permuted according to the given permutation
(R's `permute.vertexIDs`): old vertex `perm[i]` becomes new vertex `i`, and
every attribute and missing-dyad mark moves with its vertex. `perm` must be a
permutation of `1:nv(net)` (`ArgumentError` otherwise, as R's "Invalid
permutation vector").

On a **two-mode** network the permutation must map each mode onto itself
(`perm[1:n₁]` a permutation of `1:n₁`): a vertex moved across the mode boundary
would turn its cross-mode edges into impossible within-mode pairs. Such a
permutation is an `ArgumentError` (R's `permute.vertexIDs` only warns and
returns a network that violates its own partition).

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_vertex_attribute!(net, :name, ["a", "b", "c"])
p = permute_vertices(net, [3, 1, 2])    # old 3 → new 1, old 1 → new 2, old 2 → new 3
has_edge(p, 2, 3)                       # true
get_vertex_attribute(p, :name, 1)       # "c"
```
"""
function permute_vertices(net::Network{T}, perm::Vector) where T
    n = nv(net)
    (length(perm) == n && isperm(perm)) ||
        throw(ArgumentError("permute_vertices: perm must be a permutation of 1:$n " *
                            "(got a $(length(perm))-element vector" *
                            (length(perm) == n ? " that is not a permutation)" : ")")))

    if is_two_mode(net)
        k = net.bipartite
        bad = findfirst(i -> (i <= k) != (perm[i] <= k), 1:n)
        bad === nothing ||
            throw(ArgumentError("permute_vertices: on a two-mode network the " *
                                "permutation must keep each mode in place (mode 1 is " *
                                "1:$k, mode 2 is $(k + 1):$n), but it moves vertex " *
                                "$(perm[bad]) to position $bad; permute within the modes"))
    end

    # Create inverse permutation for mapping old -> new
    inv_perm = zeros(T, n)
    for (new_id, old_id) in enumerate(perm)
        inv_perm[old_id] = T(new_id)
    end

    # Create new network
    new_net = Network{T}(; n=n, directed=net.directed, loops=net.loops,
                         bipartite=net.bipartite, multiple=net.multiple, hyper=net.hyper)

    # Add permuted edges (a permutation of a valid network yields valid edges;
    # a `false` here would be a dropped edge, so it is never ignored)
    for e in edges(net)
        add_edge!(new_net, inv_perm[src(e)], inv_perm[dst(e)]) ||
            error("permute_vertices: internal error, edge ($(src(e)), $(dst(e))) " *
                  "could not be carried over")
    end

    # Permute vertex attributes
    for (attr, vals) in net.vertex_attrs
        new_vals = Dict{T, Any}(inv_perm[k] => v for (k, v) in vals)
        new_net.vertex_attrs[attr] = new_vals
    end

    # Permute edge attributes
    for (attr, vals) in net.edge_attrs
        new_vals = Dict{Tuple{T,T}, Any}()
        for ((i, j), v) in vals
            new_edge = _canonical_edge(new_net, inv_perm[i], inv_perm[j])
            new_vals[new_edge] = v
        end
        new_net.edge_attrs[attr] = new_vals
    end

    # Copy network attributes
    for (attr, val) in net.network_attrs
        new_net.network_attrs[attr] = val
    end

    # Permute missing-dyad marks
    for (i, j) in net.missing_dyads
        push!(new_net.missing_dyads, _canonical_edge(new_net, inv_perm[i], inv_perm[j]))
    end

    return new_net
end

# ============================================================================
# Graphs.jl generics that build a new graph
# ============================================================================

"""
    induced_subgraph(net::Network, vlist) -> (Network, vmap)

Graphs.jl's `induced_subgraph` on a `Network`: returns the subgraph built by
[`get_induced_subgraph`](@ref) — attributes, missing-dyad marks, `loops` and
the two-mode partition carried over — together with `vmap`, the original ID of
each vertex of the subgraph (`vmap[i] == vlist[i]`), as Graphs.jl returns it.
`vlist` may also be a `Bool` mask of length `nv(net)`. The Graphs.jl fallback
would build a bare network without any of the attributes.

# Example
```julia
using Graphs: induced_subgraph
net = network(5)
add_edges!(net, [(1, 2), (2, 3), (3, 4)])
sub, vmap = induced_subgraph(net, [2, 3, 5])
ne(sub), vmap        # (1, [2, 3, 5])
```
"""
function Graphs.induced_subgraph(net::Network, vlist::AbstractVector{U}) where U<:Integer
    return get_induced_subgraph(net, vlist), Vector{U}(vlist)
end

function Graphs.induced_subgraph(net::Network, vlist::AbstractVector{Bool})
    length(vlist) == nv(net) ||
        throw(ArgumentError("induced_subgraph: a Bool mask must have one entry per " *
                            "vertex ($(nv(net))), got $(length(vlist))"))
    return Graphs.induced_subgraph(net, findall(vlist))
end

"""
    complement(net::Network) -> Network

Graphs.jl's `complement` on a `Network`: the network with the same vertices in
which every *possible* dyad is a tie exactly when it is not one in `net` (R's
`!net`). "Possible" follows the network's own rules: self-loops are part of the
complement only when `loops=true`, and on a two-mode network only cross-mode
dyads are (within-mode pairs can never carry a tie, so the complement stays
two-mode with the same `bipartite`). Directedness, `loops`, `bipartite`, the
vertex and network attributes and the missing-dyad mask are kept — an
unobserved dyad stays unobserved, with its face value complemented — while
edge attributes are dropped, since the complement's edges are the original's
non-edges.

# Example
```julia
using Graphs: complement
net = network(3; directed=false)
add_edge!(net, 1, 2)
c = complement(net)
ne(c), has_edge(c, 1, 2), has_edge(c, 1, 3)    # (2, false, true)
```
"""
function Graphs.complement(net::Network{T, D}) where {T, D}
    n = nv(net)
    out = Network{T, D}(; n=n, loops=net.loops, bipartite=net.bipartite)
    k = net.bipartite
    for i in 1:n, j in (D ? (1:n) : (i:n))
        i == j && !net.loops && continue
        (k === nothing || ((i <= k) ⊻ (j <= k))) || continue
        Graphs.has_edge(net, i, j) && continue
        add_edge!(out, i, j) ||
            error("complement: internal error, dyad ($i, $j) could not be added")
    end
    for (attr, vals) in net.vertex_attrs
        out.vertex_attrs[attr] = copy(vals)
    end
    for (attr, val) in net.network_attrs
        out.network_attrs[attr] = val
    end
    union!(out.missing_dyads, net.missing_dyads)
    return out
end
