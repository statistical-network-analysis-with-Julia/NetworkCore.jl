# Creating Networks

This guide covers all aspects of creating, modifying, and inspecting networks in NetworkCore.jl.

## Network Constructors

### Basic Constructor

The primary way to create a network is with the `network` constructor:

```julia
using NetworkCore

# Empty network (no vertices)
net = network()

# Network with n vertices and no edges
net = network(10)

```

With an explicit type parameter, qualify the type with the module name
(after `using NetworkCore` in the REPL, the bare name `Network` refers to the
module, not the type):

<!-- skip-check -->
```julia
net = Network{Int64}(; n=10)
```

### R-Style Aliases

NetworkCore.jl provides aliases that match R's `network` package API:

```julia
# R's network(10)
net = network(10)

# R's network.initialize(10); equivalent to network(10)
net = network_initialize(10)
```

### Constructor Options

All constructors accept keyword arguments controlling network properties:

```julia
net = network(10;
    directed = true,         # Edge directionality
    bipartite = nothing,     # Set to Int for two-mode networks
    loops = false,           # Whether self-loops are allowed
    multiple = false,        # Whether multi-edges are allowed
    hyper = false,           # Whether hyperedges are allowed
)
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `directed` | `Bool` | `true` | If true, edges have direction (source to target) |
| `bipartite` | `Union{Nothing,Int}` | `nothing` | Number of mode-1 vertices for bipartite networks |
| `loops` | `Bool` | `false` | Whether vertices can have edges to themselves |
| `multiple` | `Bool` | `false` | Whether parallel edges between the same pair are allowed |
| `hyper` | `Bool` | `false` | Whether edges can connect more than two vertices |

## Directed vs. Undirected Networks

The `directed` flag fundamentally affects how the network behaves:

### Directed Networks

In directed networks, each edge has a source (sender) and a target (receiver). Edge `(1, 2)` is distinct from edge `(2, 1)`:

```julia
net = network(3; directed=true)
add_edge!(net, 1, 2)

has_edge(net, 1, 2)  # true
has_edge(net, 2, 1)  # false  -- direction matters
ne(net)              # 1
```

Directed networks are used for asymmetric relationships like email communication, citations, advice seeking, or Twitter follows.

### Undirected Networks

In undirected networks, edges are symmetric. Adding edge `(1, 2)` automatically creates `(2, 1)`:

```julia
net = network(3; directed=false)
add_edge!(net, 1, 2)

has_edge(net, 1, 2)  # true
has_edge(net, 2, 1)  # true  -- both directions present
ne(net)              # 1     -- counted once, not twice
```

Undirected networks are used for symmetric relationships like friendship (in some contexts), co-authorship, physical proximity, or shared membership.

### Implementation Detail

Internally, NetworkCore.jl always stores edges in a `SimpleDiGraph`. For undirected networks, both `(i, j)` and `(j, i)` are stored, but `ne()` divides by 2 and `edges()` only returns edges where `src <= dst` to avoid double-counting. Edge attributes use a canonical ordering via `minmax(i, j)` to ensure consistent lookups regardless of which direction you specify.

!!! warning "Call Graphs.jl on `net`, not on `net.graph`"
    Because of this storage, a Graphs.jl algorithm run on the internal field
    `net.graph` of an undirected network sees every tie twice, as two arcs:
    `Graphs.core_number(net.graph)` doubles every core number and
    `Graphs.degree(net.graph)` every degree. Pass the `Network` itself
    (`Graphs.core_number(net)`), or convert with `Graphs.SimpleGraph(net)`
    when you need a concrete undirected Graphs.jl graph. See the
    [Graph Interface](../api/graph_interface.md) page.

## Adding and Removing Vertices

### Adding Vertices

```julia
net = network(3)
nv(net)  # 3

# Add a single vertex
add_vertex!(net)
nv(net)  # 4

# Add multiple vertices
add_vertices!(net, 5)
nv(net)  # 9
```

New vertices are numbered sequentially starting after the current maximum vertex ID.

### Removing Vertices

```julia
net = network(5)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

# Remove vertex 2 (also removes edges incident to vertex 2)
rem_vertex!(net, 2)
nv(net)  # 4
```

When a vertex is removed:

1. All edges incident to that vertex are removed
2. All attributes for that vertex are deleted
3. Edge attributes for removed edges are deleted
4. On a one-mode network the last vertex is renumbered to fill the gap
   (Graphs.jl's swap-and-pop convention). On a two-mode network the vertices
   after the removed one each move down by one instead, and `bipartite` drops
   by one when a mode-1 vertex goes — R's `delete.vertices` — so the two modes
   stay contiguous blocks.

To remove several vertices at once with R's stable renumbering on any
network, use Graphs.jl's `rem_vertices!`, which returns the original ID of
each remaining vertex:

```julia
using Graphs: rem_vertices!
net = network(5)
add_edges!(net, [(1, 2), (2, 5), (4, 5)])
vmap = rem_vertices!(net, [1, 3])   # [2, 4, 5]
has_edge(net, 1, 3)                 # true — old edge (2, 5)
```

!!! note "Vertex Renumbering"
    When removing vertex `v`, vertex `nv(net)` is renumbered to `v` (the Graphs.jl swap-and-pop convention). NetworkCore.jl remaps the renumbered vertex's attributes and its incident edge attributes automatically, so attribute data stays attached to the right vertices — but any vertex IDs you hold in your own variables still refer to the old numbering.

## Adding and Removing Edges

### Adding Single Edges

```julia
net = network(5)

# Add an edge (returns true if successful, false if already exists)
add_edge!(net, 1, 2)   # true
add_edge!(net, 1, 2)   # false -- edge already exists

# Add an edge with attributes
add_edge!(net, 2, 3, Dict(:weight => 1.5, :type => :friendship))
```

### Self-Loops

By default, self-loops are not allowed:

```julia
net = network(5; loops=false)
add_edge!(net, 1, 1)   # false -- self-loop rejected

net = network(5; loops=true)
add_edge!(net, 1, 1)   # true -- self-loop allowed
```

### Adding Multiple Edges

```julia
net = network(5)

# From a list of tuples
edges_to_add = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
n_added = add_edges!(net, edges_to_add)  # returns 5
```

### Removing Edges

```julia
# Remove a single edge
rem_edge!(net, 1, 2)   # true if removed, false if not found

# For undirected networks, rem_edge! removes both directions
net = network(3; directed=false)
add_edge!(net, 1, 2)
rem_edge!(net, 1, 2)
has_edge(net, 2, 1)     # false -- both directions removed
```

When an edge is removed, all its attributes are also deleted.

## Bipartite Networks

Bipartite (two-mode) networks have two distinct vertex sets where edges only connect vertices from different sets.

### Creating Bipartite Networks

```julia
# 4 actors x 3 events
bnet = BipartiteNetwork(4, 3)

# Vertices 1-4 are mode 1 (actors)
# Vertices 5-7 are mode 2 (events)
bnet.n_mode1  # 4
bnet.n_mode2  # 3
nv(bnet)      # 7
```

### Adding Edges

Edges should connect vertices from different modes:

```julia
add_edge!(bnet, 1, 5)   # Actor 1 -- Event 5
add_edge!(bnet, 1, 6)   # Actor 1 -- Event 6
add_edge!(bnet, 2, 5)   # Actor 2 -- Event 5
add_edge!(bnet, 3, 7)   # Actor 3 -- Event 7
add_edge!(bnet, 4, 6)   # Actor 4 -- Event 6
add_edge!(bnet, 4, 7)   # Actor 4 -- Event 7
```

### Checking Bipartiteness

Two related but distinct notions are exposed:

- `Graphs.is_bipartite(net)` computes the **graph-theoretic** property:
  can the vertices be two-colored so that every edge joins the two color
  classes? (An edgeless graph is trivially bipartite.)
- `is_two_mode(net)` reports the **two-mode metadata flag**: was the
  network declared bipartite via `network(n; bipartite=k)` or built as a
  `BipartiteNetwork`?

```julia
using Graphs: is_bipartite

is_bipartite(bnet)  # true  -- only cross-mode edges are admitted
is_two_mode(bnet)   # true  -- carries mode metadata

# A one-mode network can still be bipartite as a graph...
net = network(4; directed=false)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
is_bipartite(net)   # true  -- a path is two-colorable
is_two_mode(net)    # false -- no mode metadata

# ...and stops being bipartite once it contains an odd cycle
add_edge!(net, 1, 3)
is_bipartite(net)   # false

# Marking a Network as two-mode via the constructor sets the metadata flag
net = network(7; bipartite=4)  # 4 vertices in mode 1
is_two_mode(net)    # true
```

### BipartiteNetwork Properties

```julia
bnet = BipartiteNetwork(4, 3)
println(bnet)      # BipartiteNetwork{Int64}(4+3 vertices, 0 edges)
bnet.n_mode1       # 4
bnet.n_mode2       # 3
display(bnet)
# BipartiteNetwork{Int64}: directed two-mode network
#   Vertices: 7
#   Edges: 0
#   Directed: true
#   Loops: false
#   Two-mode: 4 mode-1 vertices (1:4), 3 mode-2 vertices (5:7)
```

## Network Properties

### Size and Density

```julia
net = network(10)
add_edges!(net, [(1,2), (2,3), (3,4), (4,5), (1,5)])

nv(net)                 # 10 -- number of vertices
ne(net)                 # 5  -- number of edges
network_size(net)       # 10 -- alias for nv()
network_edgecount(net)  # 5  -- alias for ne()
network_density(net)    # 0.0556 -- 5 / (10*9) for directed
```

### Density Calculation

`network_density` follows R's `network.density`: the edge count `m` over the
number of dyads that could carry a tie. The denominator depends on the
network type and on two keywords:

| Network Type | Formula | Maximum Edges |
|-------------|---------|---------------|
| Directed, no loops | `m / (n * (n-1))` | `n * (n-1)` |
| Directed, with loops | `m / (n * n)` | `n^2` |
| Undirected, no loops | `m / (n * (n-1) / 2)` | `n * (n-1) / 2` |
| Undirected, with loops | `m / (n * (n+1) / 2)` | `n * (n+1) / 2` |
| Two-mode (`bipartite=n₁`), default | as the one-mode row above | all `n (n-1)` (or `/ 2`) dyads, within-mode ones included — R's `discount.bipartite=FALSE` |
| Two-mode, undirected, `discount_bipartite=true` | `m / (n₁ * n₂)` | `n₁ * n₂` cross-mode dyads |
| Two-mode, directed, `discount_bipartite=true` | `m / (2 * n₁ * n₂)` | `2 * n₁ * n₂` — one per direction |

- **`diag`** decides whether the diagonal (self-loops) is part of the dyad
  count. It defaults to the network's `loops` flag, exactly as R keys on
  `has.loops(x)`, so the "with loops" rows apply only to networks built with
  `loops=true`. `diag=false` on such a network mirrors `sna::gden(diag=FALSE)`:
  self-loops are dropped from both `m` and the denominator.
- **`discount_bipartite`** is a no-op on a one-mode network. The
  `BipartiteNetwork` wrapper and a `network(n; bipartite=n₁)` give the same
  answer for the same edge set.
- A network with **masked (unobserved) dyads is refused** by default; pass
  `missing=:face` to count edges at their stored face value (R's
  `na.omit=FALSE`). See [Missing Data](missing_data.md).

```julia
b = network(5; bipartite=2, directed=false)
add_edges!(b, [(1,3), (1,4), (2,5)])
network_density(b)                            # 0.3  -- 3 / 10, R's default
network_density(b; discount_bipartite=true)   # 0.5  -- 3 / (2*3)

l = network(3; loops=true)
add_edges!(l, [(1,1), (1,2)])
network_density(l)                            # 0.2222 -- 2 / 9
network_density(l; diag=false)                # 0.1667 -- 1 / 6, self-loop dropped
```

### Directedness

```julia
is_directed(net)    # true or false (a compile-time property of Network{T,D})
is_bipartite(net)   # graph-theoretic two-colorability
is_two_mode(net)    # two-mode metadata flag
```

### Vertex and Edge Iteration

```julia
# Iterate vertices
for v in vertices(net)
    println("Vertex: ", v)
end

# Iterate edges (src/dst are Graphs.jl accessors)
using Graphs: src, dst
for e in edges(net)
    println(src(e), " → ", dst(e))
end

# Existence checks
has_vertex(net, 3)    # true if vertex 3 exists
has_edge(net, 1, 2)   # true if edge (1,2) exists
```

## Neighborhoods and Subgraphs

### Neighborhoods

Find all vertices within a given distance of a focal vertex:

```julia
net = network(6; directed=false)
add_edges!(net, [(1,2), (2,3), (3,4), (4,5), (5,6)])

# Direct neighbors of vertex 1 (order=1, default)
nbhood1 = get_neighborhood(net, 1)
# Set{Int64} with 2 elements: {1, 2}

# 2-neighborhood (distance <= 2)
nbhood2 = get_neighborhood(net, 1, 2)
# Set{Int64} with 3 elements: {1, 2, 3}

# 3-neighborhood
nbhood3 = get_neighborhood(net, 1, 3)
# Set{Int64} with 4 elements: {1, 2, 3, 4}
```

The returned set always includes the focal vertex itself.

### Induced Subgraphs

Extract a subnetwork containing only specified vertices and edges between them:

```julia
net = network(6)
add_edges!(net, [(1,2), (1,3), (2,3), (3,4), (4,5), (5,6)])
set_vertex_attribute!(net, :name,
    Dict(1=>"A", 2=>"B", 3=>"C", 4=>"D", 5=>"E", 6=>"F"))
set_edge_attribute!(net, :weight,
    Dict((1,2)=>1.0, (1,3)=>2.0, (2,3)=>1.5, (3,4)=>3.0, (4,5)=>2.0, (5,6)=>1.0))

# Extract subgraph with vertices {1, 2, 3}
sub = get_induced_subgraph(net, [1, 2, 3])
nv(sub)  # 3
ne(sub)  # 3 (edges 1→2, 1→3, 2→3)

# Vertex and edge attributes are preserved (with remapped IDs)
get_vertex_attribute(sub, :name)  # Dict(1=>"A", 2=>"B", 3=>"C")
```

Note that vertex IDs are remapped to `1:length(vlist)` in the subgraph, in the order of `vlist` (R's `get.inducedSubgraph` keeps the original order instead; pass `sort(vlist)` to match it). Attributes are carried over with the new IDs. On a two-mode network the subgraph is two-mode as well, provided the selected mode-1 vertices come first in `vlist`; an interleaved selection is an `ArgumentError`.

Graphs.jl's own `induced_subgraph(net, vlist)` returns the same `Network`
together with the vertex map, and `complement(net)` the network of non-ties
(self-loops only when `loops=true`, cross-mode dyads only on a two-mode
network):

```julia
using Graphs: induced_subgraph, complement
sub, vmap = induced_subgraph(net, [1, 2, 3])
ne(complement(sub))   # 3 — the 6 ordered pairs minus the 3 arcs
```

## Vertex Permutation

Reorder vertex IDs according to a permutation:

```julia
net = network(4)
add_edges!(net, [(1,2), (2,3), (3,4)])
set_vertex_attribute!(net, :name, Dict(1=>"A", 2=>"B", 3=>"C", 4=>"D"))

# Reverse the vertex order
perm = [4, 3, 2, 1]
net_perm = permute_vertices(net, perm)

# Vertex 1 in the new network was vertex 4 in the old
get_vertex_attribute(net_perm, :name, 1)  # "D"
get_vertex_attribute(net_perm, :name, 4)  # "A"
```

The permutation vector `perm` specifies which old vertex occupies each new position: `perm[new_id] = old_id`.

## Display

A network has two displays, as R's `network` has `print` and its one-line
summary. At the REPL (or under `display`) it prints the block R users expect
from `print.network` — including `Directed:`, `Loops:` and, for a two-mode
network, the mode sizes, the three things most often mistaken — with
attribute names sorted:

```julia
net = network(5)
add_edges!(net, [(1,2), (2,3)])
set_vertex_attribute!(net, :name, Dict(1=>"Alice", 2=>"Bob"))
set_edge_attribute!(net, :weight, Dict((1,2) => 1.5))
set_network_attribute!(net, :title, "Test")

display(net)
# Network{Int64}: directed network
#   Vertices: 5
#   Edges: 2
#   Directed: true
#   Loops: false
#   Vertex attributes: name
#   Edge attributes: weight
#   Network attributes: title
```

Inside a container — a `Vector` of panel waves, a `Dict`, the `NamedTuple`
that `load_dataset(:s50)` returns — and from `println`/`string`, it is one
line, so a panel of three networks reads as three lines rather than a wall:

```julia
println(net)                       # Network{Int64,true}(5 vertices, 2 edges)
println(network(3; loops=true))    # Network{Int64,true}(3 vertices, 0 edges, loops)
println(network(6; bipartite=2))   # Network{Int64,true}(6 vertices, 0 edges, two-mode 2+4)
s50 = load_dataset(:s50)
s50.friendship                     # 3-element Vector{Network{Int64, true}}: one line per wave
```

## Creating Networks from External Data

Networks can also be created from matrices, edge lists, DataFrames, and files. See:

- [Conversion](conversion.md) for `network_from_matrix`, `network_from_edgelist`, and `network_from_dataframe`
- [I/O](io.md) for `read_pajek` and other file formats

## Random Graph Generation

For random graph generation, use the companion [SNA.jl](https://github.com/statistical-network-analysis-with-Julia/SNA.jl) package which provides `rgraph`, `rgnm`, and `rgnp` functions:

```julia
using SNA

# Erdos-Renyi G(n,p) random graph
net = rgnp(100, 0.05)

# G(n,m) random graph with exactly m edges
net = rgnm(100, 200)
```

Or use Graphs.jl's random graph generators and convert:

```julia
using Graphs

g = erdos_renyi(100, 0.05)
# Use g directly -- NetworkCore.jl is compatible via AbstractGraph
```
