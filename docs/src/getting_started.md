# Getting Started

Build a network, attach actor and edge attributes, and export it without losing track of what each value means. This tutorial uses small constructed networks before showing two-mode data and Graphs.jl interoperability.

!!! note "Before you begin"

    Static directed, undirected, and two-mode networks are supported. Edge values are attributes, with interpretation left to the analysis package. Multiple edges and hyperedges are not implemented. A missing dyad is unobserved, rather than an observed absence; consult the [missing-data guide](guide/missing_data.md) before analysis.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

Run the blocks below in order in that environment. They build on variables from earlier steps; stochastic examples use seeded random number generators where shown.

## Basic Workflow

The typical NetworkCore.jl workflow consists of four steps:

1. **Create a network** - Initialize with a number of vertices
2. **Add edges and attributes** - Build the structure and attach metadata
3. **Query the network** - Inspect properties, find neighbors, compute statistics
4. **Convert and export** - Output to matrices, DataFrames, or files

## Step 1: Create a Network

Networks are created by specifying the number of vertices:

```julia
using NetworkCore

# Create a directed network with 5 vertices — `Network(n)` is the primary constructor
net = Network(5)

# `network(5)` is the R-style alias (statnet's `network()`); both build the same object
net = network(5)

# Create an undirected network
net_undir = Network(5; directed=false)
```

### Directed vs. Undirected

By default, networks are directed (edges go from source to target). Set `directed=false` for undirected networks where edges are symmetric:

```julia
# Directed: edge 1→2 does NOT imply 2→1
net_dir = network(5; directed=true)
add_edge!(net_dir, 1, 2)
has_edge(net_dir, 1, 2)  # true
has_edge(net_dir, 2, 1)  # false

# Undirected: edge 1-2 implies both directions
net_undir = network(5; directed=false)
add_edge!(net_undir, 1, 2)
has_edge(net_undir, 1, 2)  # true
has_edge(net_undir, 2, 1)  # true
```

### Network Options

The constructor supports several options matching R's `network()`:

```julia
net = network(10;
    directed = true,       # Directed edges (default: true)
    loops = false,         # Allow self-loops (default: false)
    bipartite = nothing,   # Number of mode-1 vertices if two-mode (default: nothing)
    multiple = false,      # Must be false — parallel edges are not supported (true throws)
    hyper = false,         # Must be false — hyperedges are not supported (true throws)
)
```

`multiple` and `hyper` exist only so that R code translates keyword for
keyword: the `SimpleDiGraph` backing store cannot represent parallel or
hyper-edges, and the constructor throws rather than silently collapsing them.

### Alternative Constructors

```julia
# R-style alias
net = network(10; directed=false)

# R's network.initialize() alias
net = network_initialize(10)
```

## Step 2: Add Edges and Attributes

### Adding Edges

```julia
net = network(5)

# Add individual edges
add_edge!(net, 1, 2)  # returns true if added
add_edge!(net, 2, 3)
add_edge!(net, 3, 4)
add_edge!(net, 4, 5)
add_edge!(net, 5, 1)  # creates a cycle

# Add multiple edges at once
add_edges!(net, [(1, 3), (2, 4), (3, 5)])

# Add an edge with attributes
add_edge!(net, 1, 4, Dict(:weight => 2.5, :type => :friendship))
```

### Removing Edges

```julia
rem_edge!(net, 1, 2)  # returns true if removed
add_edge!(net, 1, 2)  # put it back: the attribute examples below address it
```

### Setting Vertex Attributes

Vertex attributes store metadata about each vertex (node):

```julia
# Set from a Dict
set_vertex_attribute!(net, :name,
    Dict(1 => "Alice", 2 => "Bob", 3 => "Carol", 4 => "David", 5 => "Eve"))

# Set from a Vector (one value per vertex, in order)
set_vertex_attribute!(net, :age, [25, 30, 28, 35, 22])

# Set for a single vertex
set_vertex_attribute!(net, :role, 1, "manager")
set_vertex_attribute!(net, :role, 2, "analyst")
```

### Setting Edge Attributes

Edge attributes store metadata about relationships:

```julia
# Set from a Dict of (source, target) => value
set_edge_attribute!(net, :weight,
    Dict((1,2) => 1.0, (2,3) => 2.0, (3,4) => 1.5))

# Set for a single edge
set_edge_attribute!(net, :type, 1, 2, :friendship)
set_edge_attribute!(net, :type, 2, 3, :collaboration)
```

### Setting Network Attributes

Network-level attributes describe the network as a whole:

```julia
set_network_attribute!(net, :title, "Office Communication Network")
set_network_attribute!(net, :collected, "2024-01-15")
set_network_attribute!(net, :organization, "Acme Corp")
```

## Step 3: Query the Network

### Basic Properties

```julia
nv(net)                # Number of vertices
ne(net)                # Number of edges
is_directed(net)       # Whether the network is directed
network_density(net)   # Edge density (proportion of possible edges)
```

### Vertex and Edge Queries

```julia
using Graphs: src, dst   # edge accessors from Graphs.jl

# Check membership
has_vertex(net, 3)     # true
has_edge(net, 1, 2)    # true

# Iterate over vertices and edges
for v in vertices(net)
    println("Vertex: ", v)
end

for e in edges(net)
    println("Edge: ", src(e), " → ", dst(e))
end
```

### Neighbors

```julia
# All neighbors (outgoing for directed networks)
collect(neighbors(net, 1))

# Incoming neighbors
collect(inneighbors(net, 1))

# Outgoing neighbors
collect(outneighbors(net, 1))
```

### Retrieving Attributes

```julia
# Get all values for an attribute
names = get_vertex_attribute(net, :name)     # Dict{Int, Any}
weights = get_edge_attribute(net, :weight)   # Dict{Tuple{Int,Int}, Any}

# Get a single value
get_vertex_attribute(net, :name, 1)          # "Alice"
get_edge_attribute(net, :weight, 1, 2)       # 1.0
get_network_attribute(net, :title)           # "Office Communication Network"

# List all attribute names
list_vertex_attributes(net)   # [:name, :age, :role]
list_edge_attributes(net)     # [:weight, :type]
list_network_attributes(net)  # [:title, :collected, :organization]
```

### Neighborhoods and Subgraphs

```julia
# Get all vertices within distance 2 of vertex 1
neighborhood = get_neighborhood(net, 1, 2)

# Extract a subgraph induced by a set of vertices
sub = get_induced_subgraph(net, [1, 2, 3])
```

## Step 4: Convert and Export

### To Adjacency Matrix

```julia
# Binary adjacency matrix
A = as_matrix(net)

# Weighted adjacency matrix (using an edge attribute)
W = as_matrix(net; attr=:weight)

# Sparse matrix (for large networks)
S = as_matrix(net; sparse=true)
```

### To Edge List

```julia
# Edge list as a matrix (n_edges x 2)
el = as_edgelist(net)

# With edge attributes
el, attrs = as_edgelist(net; attrs=[:weight, :type])
```

### To DataFrame

Requires DataFrames.jl to be loaded (the conversion methods live in a
package extension):

```julia
using DataFrames

# Edge DataFrame
edge_df = as_dataframe(net)

# Vertex DataFrame
vertex_df = as_dataframe(net; vertices=true)
```

### To File

```julia
# Write to Pajek format
write_pajek(net, "network.net")

# Write to GraphML
write_graphml(net, "network.graphml")

# Write to CSV edge list
write_edgelist_csv(net, "edges.csv")
```

## Complete Example

Here is a full example building and analyzing a small social network:

```julia
using NetworkCore

# Create a directed friendship network
net = network(6; directed=true)

# Add friendship ties
friendships = [
    (1, 2), (1, 3),       # Alice -> Bob, Carol
    (2, 1), (2, 3),       # Bob -> Alice, Carol
    (3, 4),                # Carol -> David
    (4, 5), (4, 6),       # David -> Eve, Frank
    (5, 6),                # Eve -> Frank
    (6, 4),                # Frank -> David
]
add_edges!(net, friendships)

# Add vertex names
set_vertex_attribute!(net, :name,
    Dict(1=>"Alice", 2=>"Bob", 3=>"Carol",
         4=>"David", 5=>"Eve", 6=>"Frank"))

# Add vertex department
set_vertex_attribute!(net, :department,
    Dict(1=>"Engineering", 2=>"Engineering", 3=>"Marketing",
         4=>"Marketing", 5=>"Sales", 6=>"Sales"))

# Add edge weights (interaction frequency)
set_edge_attribute!(net, :frequency,
    Dict((1,2)=>5, (1,3)=>3, (2,1)=>4, (2,3)=>2,
         (3,4)=>6, (4,5)=>3, (4,6)=>1, (5,6)=>7, (6,4)=>2))

# Query properties
println("Network Summary:")
println("  Vertices: ", nv(net))
println("  Edges: ", ne(net))
println("  Density: ", round(network_density(net), digits=3))
println("  Directed: ", is_directed(net))
println()

# Inspect a vertex
println("Alice's outgoing ties:")
for v in outneighbors(net, 1)
    name = get_vertex_attribute(net, :name, v)
    freq = get_edge_attribute(net, :frequency, 1, v)
    println("  → $name (frequency: $freq)")
end
println()

# Find David's 2-neighborhood
nbhood = get_neighborhood(net, 4, 2)
println("Vertices within distance 2 of David: ",
    [get_vertex_attribute(net, :name, v) for v in nbhood])

# Extract the engineering subnetwork
eng_vertices = [v for v in vertices(net)
                if get_vertex_attribute(net, :department, v) == "Engineering"]
eng_sub = get_induced_subgraph(net, eng_vertices)
println("\nEngineering subnetwork: $(nv(eng_sub)) vertices, $(ne(eng_sub)) edges")

# Convert to adjacency matrix
A = as_matrix(net)
println("\nAdjacency matrix:")
display(A)

# Export
write_pajek(net, "friendship.net")
println("\nNetwork written to friendship.net")
```

## Working with Bipartite Networks

Bipartite (two-mode) networks have two distinct sets of vertices where edges only connect vertices from different sets:

```julia
using Graphs   # for Graphs.is_bipartite below

# Create a bipartite network: 3 people x 4 events
bnet = BipartiteNetwork(3, 4)

# Vertices 1-3 are people, vertices 4-7 are events
# Add attendance edges
add_edge!(bnet, 1, 4)   # Person 1 attended Event 4
add_edge!(bnet, 1, 5)   # Person 1 attended Event 5
add_edge!(bnet, 2, 4)   # Person 2 attended Event 4
add_edge!(bnet, 2, 6)   # Person 2 attended Event 6
add_edge!(bnet, 3, 5)   # Person 3 attended Event 5
add_edge!(bnet, 3, 7)   # Person 3 attended Event 7

println("Bipartite network:")
println("  Mode 1 (people): ", bnet.n_mode1)
println("  Mode 2 (events): ", bnet.n_mode2)
println("  Edges: ", ne(bnet))
println("  Is bipartite: ", Graphs.is_bipartite(bnet))
```

## Using Graphs.jl Algorithms

Because `Network` implements `AbstractGraph`, you can use any Graphs.jl algorithm directly:

```julia
using Graphs

net = network(5; directed=false)
add_edges!(net, [(1,2), (2,3), (3,4), (4,5), (5,1), (1,3)])

# Shortest paths
dists = gdistances(net, 1)
println("Distances from vertex 1: ", dists)

# Connected components
comps = connected_components(net)
println("Components: ", comps)

# Clustering coefficient
cc = global_clustering_coefficient(net)
println("Clustering coefficient: ", cc)

# Degree centrality
for v in vertices(net)
    println("Vertex $v degree: ", degree(net, v))
end
```

## Creating Networks from Data

### From an Adjacency Matrix

```julia
A = [0 1 1 0;
     0 0 1 1;
     0 0 0 1;
     1 0 0 0]
net = network_from_matrix(A)
```

### From an Edge List

```julia
edges = [(1, 2), (2, 3), (3, 4), (4, 1)]
net = network_from_edgelist(edges)
```

### From a DataFrame

```julia
using DataFrames

df = DataFrame(
    source = [1, 1, 2, 3],
    target = [2, 3, 3, 4],
    weight = [1.0, 2.0, 1.5, 3.0]
)
net = network_from_dataframe(df)
# Edge weights are automatically stored as the :weight attribute
```

### From a Pajek File

<!-- skip-check -->
```julia
net = read_pajek("karate.net")
```

## Best Practices

1. **Choose directed vs. undirected early**: This affects edge counting, neighbor queries, and attribute storage. Changing later requires recreating the network.

2. **Use Symbol keys for attributes**: Attribute names are Julia `Symbol`s (e.g., `:name`, `:weight`), not strings.

3. **Leverage Graphs.jl**: Do not reimplement standard graph algorithms. Since `Network` is an `AbstractGraph`, use Graphs.jl's extensive library.

4. **Use sparse matrices for large networks**: When converting networks with thousands of vertices, pass `sparse=true` to `as_matrix()` to avoid memory issues.

5. **Prefer bulk operations**: Use `add_edges!` and the Dict/Vector forms of `set_vertex_attribute!` / `set_edge_attribute!` when adding many elements, rather than looping over individual calls. The bulk forms validate exactly as the single-element forms do (every endpoint must be a vertex, every edge-attribute key must be an edge unless `require_edge=false`), and all-or-nothing: a bad element means nothing is added or stored.

6. **Keep attribute names consistent**: Use snake_case Symbol names (e.g., `:vertex_name`, `:edge_weight`) for clarity.

## Next Steps

- Learn about [Creating Networks](guide/construction.md) in detail
- Explore [Attributes](guide/attributes.md) for rich metadata management
- Understand [Conversion](guide/conversion.md) between formats
- Work with [I/O](guide/io.md) for reading and writing files
