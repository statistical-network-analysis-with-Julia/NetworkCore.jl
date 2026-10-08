# I/O

NetworkCore.jl provides functions for reading and writing networks in several file formats. The primary supported format is Pajek (`.net`), which is widely used in social network analysis. Basic support is also provided for GraphML and CSV edge list output.

## Pajek Format

The [Pajek](http://vlado.fmf.uni-lj.si/pub/networks/pajek/) format is a widely used plain-text format for network data. It supports vertex labels, edge weights, and both directed and undirected networks.

### Pajek File Structure

A typical Pajek `.net` file looks like:

```text
*Vertices 5
1 "Alice"
2 "Bob"
3 "Carol"
4 "David"
5 "Eve"
*Arcs
1 2 1.0
1 3 2.0
2 3 1.5
3 4 3.0
4 5 0.5
```

Key sections:

| Section | Description |
|---------|-------------|
| `*Vertices n` | Declares n vertices, followed by optional vertex labels |
| `*Arcs` | Directed edges (source target weight) |
| `*Edges` | Undirected edges (vertex1 vertex2 weight) |

Lines starting with `%` are comments and are ignored.

### Reading Pajek Files

```julia
using NetworkCore

# A small Pajek file to work with
write("social_network.net", """
*Vertices 3
1 "alice"
2 "bob"
3 "carol"
*Arcs
1 2 1.0
2 3 1.0
""")

# Read a directed network
net = read_pajek("social_network.net")

# Check what was loaded
println(nv(net))           # Number of vertices
println(ne(net))           # Number of edges
println(is_directed(net))  # true if *Arcs, false if *Edges
```

#### What Gets Loaded

| Feature | Stored As |
|---------|-----------|
| Vertex labels | Vertex attribute `:vertex_names` |
| Edge weights | Edge attribute `:weight` (only if any weight differs from 1.0) |
| Directedness | `net.directed` (true for `*Arcs`, false for `*Edges`) |

#### Example: Loading and Inspecting

```julia
# A sample file with labeled vertices and weighted edges
write("karate.net", """
*Vertices 3
1 "Mr Hi"
2 "Actor 2"
3 "John A"
*Edges
1 2 4.0
2 3 5.0
""")

net = read_pajek("karate.net")

# Get vertex labels
names = get_vertex_attribute(net, :vertex_names)
for (v, name) in sort(collect(names))
    println("Vertex $v: $name")
end

# Check for edge weights
weights = get_edge_attribute(net, :weight)
if !isempty(weights)
    println("Network has edge weights")
    for ((i, j), w) in weights
        name_i = get(names, i, string(i))
        name_j = get(names, j, string(j))
        println("  $name_i → $name_j: $w")
    end
end
```

#### Handling Different Encodings

The reader expects UTF-8 encoded files. If your Pajek file uses a different encoding, convert it first:

<!-- skip-check -->
```julia
# Read with explicit encoding handling
lines = open("network_latin1.net", "r") do io
    readlines(io)
end
# Save as UTF-8 and re-read
```

### Writing Pajek Files

```julia
net = network(4)
add_edges!(net, [(1,2), (2,3), (3,4), (4,1)])

# Set vertex labels
set_vertex_attribute!(net, :vertex_names,
    Dict(1=>"Alice", 2=>"Bob", 3=>"Carol", 4=>"David"))

# Set edge weights
set_edge_attribute!(net, :weight,
    Dict((1,2)=>1.0, (2,3)=>2.5, (3,4)=>0.8, (4,1)=>1.2))

# Write to file
write_pajek(net, "output.net")
```

The resulting file:

```text
*Vertices 4
1 "Alice"
2 "Bob"
3 "Carol"
4 "David"
*Arcs
1 2 1.0
2 3 2.5
3 4 0.8
4 1 1.2
```

#### Custom Attribute Names

By default, `write_pajek` uses `:vertex_names` for labels and `:weight` for edge weights. You can specify different attributes:

```julia
# Store names under a different attribute
set_vertex_attribute!(net, :label,
    Dict(1=>"Node A", 2=>"Node B", 3=>"Node C", 4=>"Node D"))

# Store weights under a different attribute
set_edge_attribute!(net, :strength,
    Dict((1,2)=>5.0, (2,3)=>3.0, (3,4)=>7.0, (4,1)=>2.0))

# Write using custom attributes
write_pajek(net, "output.net"; vertex_names=:label, edge_weight=:strength)
```

#### Directed vs. Undirected

The writer uses `*Arcs` for directed networks and `*Edges` for undirected:

```julia
# Directed network → *Arcs section
net_dir = network(3; directed=true)
add_edges!(net_dir, [(1,2), (2,3)])
write_pajek(net_dir, "directed.net")
# Produces: *Arcs

# Undirected network → *Edges section
net_undir = network(3; directed=false)
add_edges!(net_undir, [(1,2), (2,3)])
write_pajek(net_undir, "undirected.net")
# Produces: *Edges
```

#### Vertices Without Labels

Vertices without a name in the `:vertex_names` attribute are written with just their numeric ID:

```text
*Vertices 4
1 "Alice"
2
3 "Carol"
4
```

### Function Signatures

<!-- skip-check -->
```julia
read_pajek(filepath::String) -> Network

write_pajek(net::Network, filepath::String;
    vertex_names::Symbol = :vertex_names,  # Vertex attribute for labels
    edge_weight::Symbol = :weight          # Edge attribute for weights
)
```

## GraphML Format

[GraphML](http://graphml.graphdrawing.org/) is an XML-based format for graph data that supports arbitrary attributes. NetworkCore.jl provides a basic writer for GraphML.

### Writing GraphML

```julia
net = network(3)
add_edges!(net, [(1,2), (2,3), (3,1)])

# Set attributes
set_vertex_attribute!(net, :name, Dict(1=>"Alice", 2=>"Bob", 3=>"Carol"))
set_vertex_attribute!(net, :age, Dict(1=>25, 2=>30, 3=>28))
set_edge_attribute!(net, :weight, Dict((1,2)=>1.0, (2,3)=>2.5, (3,1)=>0.8))

# Write
write_graphml(net, "network.graphml")
```

The resulting XML:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
  <key id="v_name" for="node" attr.name="name" attr.type="string"/>
  <key id="v_age" for="node" attr.name="age" attr.type="string"/>
  <key id="e_weight" for="edge" attr.name="weight" attr.type="string"/>
  <graph id="G" edgedefault="directed">
    <node id="n1">
      <data key="v_name">Alice</data>
      <data key="v_age">25</data>
    </node>
    <node id="n2">
      <data key="v_name">Bob</data>
      <data key="v_age">30</data>
    </node>
    <node id="n3">
      <data key="v_name">Carol</data>
      <data key="v_age">28</data>
    </node>
    <edge id="e0" source="n1" target="n2">
      <data key="e_weight">1.0</data>
    </edge>
    <edge id="e1" source="n2" target="n3">
      <data key="e_weight">2.5</data>
    </edge>
    <edge id="e2" source="n3" target="n1">
      <data key="e_weight">0.8</data>
    </edge>
  </graph>
</graphml>
```

#### Features and Limitations

- All vertex and edge attributes are included automatically
- Attribute values are converted to strings via `string()`
- Attribute types are all declared as `"string"` (basic implementation)
- Graph `edgedefault` is set to `"directed"` or `"undirected"` based on `is_directed(net)`
- Vertex IDs are written as `n1`, `n2`, etc.
- Edge IDs are written as `e0`, `e1`, etc.

!!! note "Read Support"
    NetworkCore.jl does not currently include a GraphML reader. For reading GraphML files, consider using [GraphIO.jl](https://github.com/JuliaGraphs/GraphIO.jl) and converting the result.

### Function Signature

<!-- skip-check -->
```julia
write_graphml(net::Network, filepath::String)
```

## CSV Edge List

Write the network as a comma-separated edge list, suitable for import into other tools.

### Writing CSV

```julia
net = network(4)
add_edges!(net, [(1,2), (2,3), (3,4)])
set_edge_attribute!(net, :weight, Dict((1,2)=>1.0, (2,3)=>2.5, (3,4)=>0.8))
set_edge_attribute!(net, :type, Dict((1,2)=>"friend", (2,3)=>"colleague", (3,4)=>"friend"))

# Write all edge attributes
write_edgelist_csv(net, "edges.csv")
```

The resulting CSV:

```text
source,target,weight,type
1,2,1.0,friend
2,3,2.5,colleague
3,4,0.8,friend
```

### Selecting Attributes

Write only specific edge attributes:

```julia
write_edgelist_csv(net, "edges_weights.csv"; attrs=[:weight])
```

Result:

```text
source,target,weight
1,2,1.0
2,3,2.5
3,4,0.8
```

### Without Attributes

If no edge attributes are set, only `source` and `target` columns are written:

```julia
net = network(3)
add_edges!(net, [(1,2), (2,3)])
write_edgelist_csv(net, "simple_edges.csv")
```

Result:

```text
source,target
1,2
2,3
```

### Missing Values

Edges that lack a value for an included attribute will have an empty field:

```text
source,target,weight
1,2,1.0
2,3,
3,4,0.8
```

### Function Signature

<!-- skip-check -->
```julia
write_edgelist_csv(net::Network, filepath::String;
    attrs::Vector{Symbol} = Symbol[]  # Edge attributes to include (empty = all)
)
```

## File Format Comparison

| Feature | Pajek (.net) | GraphML (.graphml) | CSV Edge List |
|---------|-------------|-------------------|---------------|
| Read support | Yes | No | No (use DataFrame) |
| Write support | Yes | Yes | Yes |
| Vertex labels | Yes | Yes | No |
| Edge weights | Yes | Yes | Yes |
| Arbitrary attributes | No (names + weights only) | Yes | Edge attributes only |
| Directed/undirected | Yes | Yes | Via convention |
| Human-readable | Yes | Semi (XML) | Yes |
| Interoperability | Pajek, UCINET, R | Gephi, yEd, Cytoscape | Universal |

## Reading from Other Formats

### From CSV

To read a CSV edge list, use DataFrames and the conversion functions:

```julia
using NetworkCore, DataFrames, CSV

# A CSV edge list on disk (here written by write_edgelist_csv)
tmp = mktempdir()
edges_file = joinpath(tmp, "edges.csv")
src_net = network(4)
add_edges!(src_net, [(1, 2), (2, 3), (3, 4)])
write_edgelist_csv(src_net, edges_file)

# Read CSV
df = CSV.read(edges_file, DataFrame)

# Convert to network
net = network_from_dataframe(df; source=:source, target=:target)
ne(net)   # 3
```

### From Other Graph Libraries

Convert from Graphs.jl types:

```julia
using Graphs

# Create a Graphs.jl graph
g = erdos_renyi(100, 0.05)

# Convert to Network via adjacency matrix
A = adjacency_matrix(g)
net = network_from_matrix(Matrix(A); directed=is_directed(g))
```

## Complete I/O Example

```julia
using NetworkCore
using DataFrames

# Create a network
net = network(5)
add_edges!(net, [(1,2), (1,3), (2,3), (3,4), (4,5)])

# Add rich metadata
set_vertex_attribute!(net, :vertex_names,
    Dict(1=>"Alice", 2=>"Bob", 3=>"Carol", 4=>"David", 5=>"Eve"))
set_vertex_attribute!(net, :department,
    Dict(1=>"Eng", 2=>"Eng", 3=>"Mkt", 4=>"Mkt", 5=>"Sales"))

set_edge_attribute!(net, :weight,
    Dict((1,2)=>3.0, (1,3)=>1.0, (2,3)=>2.0, (3,4)=>4.0, (4,5)=>1.5))
set_edge_attribute!(net, :type,
    Dict((1,2)=>:strong, (1,3)=>:weak, (2,3)=>:strong, (3,4)=>:strong, (4,5)=>:weak))

set_network_attribute!(net, :title, "Office Network")

# Export in multiple formats
write_pajek(net, "office.net")
write_graphml(net, "office.graphml")
write_edgelist_csv(net, "office_edges.csv")

# Also export vertex attributes separately (requires `using DataFrames`)
vertex_df = as_dataframe(net; vertices=true)
println("Vertex DataFrame:")
println(vertex_df)

# Read back from Pajek
net2 = read_pajek("office.net")
println("\nRe-loaded from Pajek:")
println("  Vertices: ", nv(net2))
println("  Edges: ", ne(net2))
println("  Names: ", get_vertex_attribute(net2, :vertex_names))
```

## Best Practices

1. **Use Pajek for network exchange**: It is the most widely supported format in social network analysis software (UCINET, Pajek, R's `network` package).

2. **Use GraphML for visualization tools**: GraphML is the preferred import format for Gephi, yEd, and Cytoscape.

3. **Use CSV for data pipelines**: CSV edge lists are easy to process with pandas, R, Excel, and other data tools.

4. **Save attributes separately when needed**: Pajek only preserves vertex names and edge weights. For full attribute preservation, also export a vertex DataFrame.

5. **Use DataFrames for reading CSV**: Rather than implementing a CSV reader, NetworkCore.jl relies on the mature CSV.jl and DataFrames.jl ecosystem.

6. **Verify round-trip fidelity**: When saving and reloading, check that the number of vertices, edges, and key attributes match the original.
