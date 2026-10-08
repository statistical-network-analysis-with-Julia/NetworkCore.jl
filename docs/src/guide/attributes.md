# Attributes

NetworkCore.jl provides a rich attribute system for attaching metadata to vertices, edges, and the network itself. This mirrors the attribute functionality in R's `network` package, where vertex and edge attributes are central to statistical network analysis.

## Attribute Levels

Attributes exist at three levels:

| Level | Stored In | Key Type | Description |
|-------|-----------|----------|-------------|
| **Vertex** | `net.vertex_attrs` | `Symbol → Dict{T, Any}` | Per-vertex metadata (name, age, group) |
| **Edge** | `net.edge_attrs` | `Symbol → Dict{Tuple{T,T}, Any}` | Per-edge metadata (weight, type, date) |
| **Network** | `net.network_attrs` | `Symbol → Any` | Whole-network metadata (title, description) |

All attribute names are Julia `Symbol`s (e.g., `:name`, `:weight`, `:department`).

## Vertex Attributes

Vertex attributes store metadata about individual vertices (nodes) in the network.

### Setting Vertex Attributes

There are three ways to set vertex attributes:

#### From a Dict

Map vertex IDs to values:

```julia
using NetworkCore

net = network(5)

# Set names for specific vertices
set_vertex_attribute!(net, :name,
    Dict(1 => "Alice", 2 => "Bob", 3 => "Carol"))

# Set ages
set_vertex_attribute!(net, :age,
    Dict(1 => 25, 2 => 30, 3 => 28, 4 => 35, 5 => 22))
```

You do not need to provide values for all vertices. Missing vertices will return `nothing` when queried.

#### From a Vector

Provide one value per vertex in order:

```julia
# Values are assigned to vertices 1, 2, ..., n
set_vertex_attribute!(net, :score, [0.8, 0.6, 0.9, 0.3, 0.7])
```

The vector length must exactly match `nv(net)`.

#### For a Single Vertex

```julia
set_vertex_attribute!(net, :role, 1, "manager")
set_vertex_attribute!(net, :role, 2, "analyst")
set_vertex_attribute!(net, :role, 3, "engineer")
```

### Getting Vertex Attributes

#### Get All Values

Returns a `Dict` mapping vertex IDs to values:

```julia
names = get_vertex_attribute(net, :name)
# Dict(1 => "Alice", 2 => "Bob", 3 => "Carol")

# Iterate
for (v, name) in names
    println("Vertex $v: $name")
end
```

#### Get a Single Value

```julia
name = get_vertex_attribute(net, :name, 1)
# "Alice"

# Returns nothing if the vertex has no value for this attribute
get_vertex_attribute(net, :name, 4)
# nothing
```

### Listing Vertex Attributes

```julia
attrs = list_vertex_attributes(net)
# [:name, :age, :score, :role]
```

### Deleting Vertex Attributes

```julia
# Remove the :score attribute entirely
delete_vertex_attribute!(net, :score)

list_vertex_attributes(net)
# [:name, :age, :role]
```

### Vertex Attribute Patterns

#### Categorical Attributes

Store group memberships:

```julia
set_vertex_attribute!(net, :department,
    Dict(1 => "Engineering", 2 => "Engineering",
         3 => "Marketing", 4 => "Sales", 5 => "Sales"))
```

#### Numeric Attributes

Store continuous measurements:

```julia
set_vertex_attribute!(net, :centrality, [0.45, 0.32, 0.61, 0.28, 0.15])
```

#### Boolean Attributes

Store binary flags:

```julia
set_vertex_attribute!(net, :is_manager,
    Dict(1 => true, 2 => false, 3 => false, 4 => true, 5 => false))
```

## Edge Attributes

Edge attributes store metadata about individual edges (ties) in the network.

### Setting Edge Attributes

#### From a Dict

Map `(source, target)` tuples to values:

```julia
net = network(4)
add_edges!(net, [(1,2), (2,3), (3,4), (4,1)])

# Set edge weights
set_edge_attribute!(net, :weight,
    Dict((1,2) => 1.0, (2,3) => 2.5, (3,4) => 0.8, (4,1) => 1.2))

# Set edge types
set_edge_attribute!(net, :type,
    Dict((1,2) => :friendship, (2,3) => :collaboration,
         (3,4) => :mentorship, (4,1) => :friendship))
```

#### For a Single Edge

```julia
set_edge_attribute!(net, :weight, 1, 2, 3.5)
set_edge_attribute!(net, :created, 1, 2, "2024-01-15")
```

### Getting Edge Attributes

#### Get All Values

Returns a `Dict` mapping `(source, target)` tuples to values:

```julia
weights = get_edge_attribute(net, :weight)
# Dict((1,2) => 1.0, (2,3) => 2.5, (3,4) => 0.8, (4,1) => 1.2)

# Iterate
for ((i, j), w) in weights
    println("Edge ($i, $j): weight = $w")
end
```

#### Get a Single Value

```julia
w = get_edge_attribute(net, :weight, 1, 2)
# 1.0

# Returns nothing if the edge has no value for this attribute
get_edge_attribute(net, :weight, 2, 1)
# nothing (edge 2→1 does not exist in a directed network)
```

### Edge Attributes in Undirected Networks

For undirected networks, edge attributes use canonical ordering (`minmax(i, j)`), so the order of vertices does not matter:

```julia
net = network(3; directed=false)
add_edge!(net, 1, 2)

# Both orders refer to the same edge
set_edge_attribute!(net, :weight, 1, 2, 5.0)
get_edge_attribute(net, :weight, 2, 1)  # 5.0  -- same edge
get_edge_attribute(net, :weight, 1, 2)  # 5.0  -- same edge
```

### Listing Edge Attributes

```julia
attrs = list_edge_attributes(net)
# [:weight, :type, :created]
```

### Deleting Edge Attributes

```julia
delete_edge_attribute!(net, :created)

list_edge_attributes(net)
# [:weight, :type]
```

### Adding Edges with Attributes

You can add an edge and set its attributes in one call:

```julia
add_edge!(net, 1, 3, Dict(:weight => 2.0, :type => :friendship))
```

## Network Attributes

Network-level attributes store metadata about the entire network.

### Setting Network Attributes

```julia
set_network_attribute!(net, :title, "Florentine Marriage Network")
set_network_attribute!(net, :collected, 2024)
set_network_attribute!(net, :source, "Padgett & Ansell (1993)")
set_network_attribute!(net, :notes, "Historical marriage ties among Florentine families")
```

### Getting Network Attributes

```julia
title = get_network_attribute(net, :title)
# "Florentine Marriage Network"

# Returns nothing if the attribute does not exist
get_network_attribute(net, :nonexistent)
# nothing
```

### Listing Network Attributes

```julia
attrs = list_network_attributes(net)
# [:title, :collected, :source, :notes]
```

### Deleting Network Attributes

```julia
delete_network_attribute!(net, :notes)

list_network_attributes(net)
# [:title, :collected, :source]
```

## Convenience Indexing

NetworkCore.jl provides bracket-style indexing that mirrors R's `%v%` and `%e%` operators:

### Reading Attributes

```julia
# Vertex attributes (level = :v)
net[:v, :name]   # equivalent to get_vertex_attribute(net, :name)

# Edge attributes (level = :e)
net[:e, :weight]  # equivalent to get_edge_attribute(net, :weight)

# Network attributes (level = :n)
net[:n, :title]   # equivalent to get_network_attribute(net, :title)
```

### Writing Attributes

```julia
# Set vertex attribute
net[:v, :name] = Dict(1 => "Alice", 2 => "Bob")

# Set edge attribute
net[:e, :weight] = Dict((1,2) => 1.0, (1,3) => 2.0)   # both are edges of net

# Set network attribute
net[:n, :title] = "My Network"
```

### R Comparison

| R Syntax | Julia Equivalent |
|----------|-----------------|
| `net %v% "name"` | `net[:v, :name]` |
| `net %e% "weight"` | `net[:e, :weight]` |
| `set.vertex.attribute(net, "name", value)` | `net[:v, :name] = value` |
| `set.edge.attribute(net, "weight", value)` | `net[:e, :weight] = value` |

## Attributes and Graph Operations

### Subgraph Extraction

When extracting an induced subgraph, vertex and edge attributes are preserved with remapped vertex IDs:

```julia
net = network(5)
add_edges!(net, [(1,2), (2,3), (3,4), (4,5)])
set_vertex_attribute!(net, :name,
    Dict(1=>"A", 2=>"B", 3=>"C", 4=>"D", 5=>"E"))
set_edge_attribute!(net, :weight,
    Dict((1,2)=>1.0, (2,3)=>2.0, (3,4)=>3.0, (4,5)=>4.0))

sub = get_induced_subgraph(net, [2, 3, 4])
# Vertices are remapped: old 2→new 1, old 3→new 2, old 4→new 3

get_vertex_attribute(sub, :name)
# Dict(1=>"B", 2=>"C", 3=>"D")

get_edge_attribute(sub, :weight)
# Dict((1,2)=>2.0, (2,3)=>3.0)
```

Network-level attributes are also copied to the subgraph.

### Vertex Permutation

When permuting vertices, attributes follow the vertices:

```julia
net = network(3)
add_edges!(net, [(1,2), (2,3)])
set_vertex_attribute!(net, :name, Dict(1=>"A", 2=>"B", 3=>"C"))

net_perm = permute_vertices(net, [3, 1, 2])
# Old vertex 3 is now vertex 1, old 1 is now 2, old 2 is now 3

get_vertex_attribute(net_perm, :name)
# Dict(1=>"C", 2=>"A", 3=>"B")
```

### Edge Removal

When an edge is removed, its attributes are automatically deleted:

```julia
net = network(3)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 5.0)

rem_edge!(net, 1, 2)
get_edge_attribute(net, :weight, 1, 2)  # nothing
```

### Vertex Removal

When a vertex is removed, all its attributes and the attributes of incident edges are deleted:

```julia
net = network(3)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
set_vertex_attribute!(net, :name, Dict(1=>"A", 2=>"B", 3=>"C"))

rem_vertex!(net, 2)
# Vertex 2 and edges (1,2), (2,3) are all removed
# Vertex 3 is renumbered to 2
```

## Attribute Data Types

Attribute values can be any Julia type (the setters address existing
vertices and edges only — an unknown vertex ID or a dyad without a tie is an
`ArgumentError`, see [Common mistakes](@ref attr-common-mistakes)):

```julia
net = network(3)
add_edge!(net, 1, 2)

# Strings
set_vertex_attribute!(net, :name, 1, "Alice")

# Numbers
set_vertex_attribute!(net, :age, 1, 30)
set_edge_attribute!(net, :weight, 1, 2, 3.14)

# Symbols
set_edge_attribute!(net, :type, 1, 2, :friendship)

# Vectors
set_vertex_attribute!(net, :coordinates, 1, [40.7128, -74.0060])

# Dicts
set_network_attribute!(net, :metadata, Dict("year" => 2024, "complete" => true))

# Custom types
struct Location
    lat::Float64
    lon::Float64
end
set_vertex_attribute!(net, :location, 1, Location(40.7128, -74.0060))
```

Since attribute storage uses `Dict{..., Any}`, any Julia object can be stored. However, for serialization (writing to files), stick to basic types (strings, numbers, symbols).

## Typed Attribute Access (Fast Path)

Because attributes are stored as `Dict{..., Any}`, reading them inside a hot
loop is type-unstable: the compiler cannot know the value type, so every
access boxes and dynamically dispatches. Model code (ERGM change statistics,
SNA centralities weighted by an attribute, ...) should instead convert once
up front using the typed accessors, then loop over the concretely typed
result:

```julia
net = network(3)
add_edge!(net, 1, 2)
set_vertex_attribute!(net, :wealth, [10, 36, 55])
set_vertex_attribute!(net, :group, 1, "core")
set_edge_attribute!(net, :weight, 1, 2, 1.5)

# Dict form: Dict{Int, Float64}, converted in a single pass
wealth = get_vertex_attribute(net, :wealth, Float64)

# Vector form: dense Vector{Float64} indexed by vertex ID (1:nv)
wealth = vertex_attribute_vector(net, :wealth, Float64)

# Vertices without a value take `default`; without a default, a missing
# value throws an ArgumentError naming the attribute and the vertex (and
# listing the attributes that exist when the name is unknown) instead of
# silently inventing data
group = vertex_attribute_vector(net, :group, String; default="none")
# ["core", "none", "none"]

# Edge attributes too: Dict{Tuple{Int,Int}, Float64} with canonical keys
weights = get_edge_attribute(net, :weight, Float64)
```

All typed accessors perform exactly one conversion pass and throw if a
stored value cannot be converted to the requested type. The untyped
accessors (`get_vertex_attribute(net, :attr)` etc.) are unchanged and remain
the right choice for interactive exploration of heterogeneous data.

## [Common mistakes](@id attr-common-mistakes)

The setters and the typed fast path refuse the mistakes R's `network` refuses
— an unknown vertex ID, an edge attribute on a dyad that carries no tie, a
mistyped attribute name — with a message that names the fix, instead of
storing a value nobody can reach by a legal ID or throwing a bare `KeyError`:

```julia
net = network(3)
add_edge!(net, 1, 2)
set_vertex_attribute!(net, :group, [1, 1, 2])

# An unknown vertex ID (R: "vertex IDs out of range"); nothing is stored
try
    set_vertex_attribute!(net, :group, 9, 1)
catch e
    println(e.msg)   # set_vertex_attribute!(net, :group, ...): vertex 9 does not exist (network has 3 vertices)
end

# A dyad without a tie (R's set.edge.attribute can only address an edge that exists)
try
    set_edge_attribute!(net, :weight, 2, 3, 1.0)
catch e
    println(e.msg)   # ... (2, 3) is not an edge of the network; add_edge!(net, 2, 3) first, or pass require_edge=false ...
end
set_edge_attribute!(net, :weight, 2, 3, 1.0; require_edge=false)   # deliberately ahead of the edge, in writing

# A mistyped attribute name in the typed fast path lists what the network does have
try
    vertex_attribute_vector(net, :grup, Int)
catch e
    println(e.msg)   # vertex attribute :grup is not set for vertex 1 (the network has no attribute :grup; existing: [:group]); pass default= to fill gaps
end
```

## Best Practices

1. **Use meaningful attribute names**: Choose descriptive Symbol names like `:department`, `:friendship_strength`, `:year_established`.

2. **Be consistent with types**: All values for a given attribute should ideally have the same type (all strings, all numbers, etc.) to avoid surprises when processing.

3. **Document your attributes**: Use network-level attributes to store a description of each vertex/edge attribute:
   ```julia
   set_network_attribute!(net, :attribute_descriptions,
       Dict(:weight => "Interaction frequency per month",
            :type => "Relationship category"))
   ```

4. **Check attribute existence**: Use `list_vertex_attributes` before assuming an attribute exists.

5. **Use bracket syntax for quick access**: The `net[:v, :attr]` syntax is concise for interactive use.

6. **Prefer bulk operations**: When setting attributes for many vertices/edges, use the Dict or Vector form of `set_vertex_attribute!` / `set_edge_attribute!` rather than looping over individual calls. The Dict forms validate exactly as the single-element forms do — every key must be a vertex, and every edge-attribute key must be an existing edge unless `require_edge=false` — and they are all-or-nothing: one bad key means nothing is stored, so the mistake surfaces where it is made.
