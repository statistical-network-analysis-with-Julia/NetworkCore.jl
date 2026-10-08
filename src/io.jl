# I/O functions for Network objects.
#
# Provides functions to read and write networks in various formats,
# including Pajek format commonly used in network analysis.

# ============================================================================
# Pajek Format
# ============================================================================

"""
    read_pajek(filepath::String) -> Network

Read a network from a Pajek (.net) file.

# Pajek Format
```
*Vertices n
1 "vertex1"
2 "vertex2"
...
*Edges (or *Arcs for directed)
1 2 weight
2 3 weight
...
```

An `*Arcs` section gives a directed network, `*Edges` an undirected one;
vertex labels are stored as the `:vertex_names` vertex attribute and edge
weights (when present) as the `:weight` edge attribute — the same names
[`write_pajek`](@ref) reads, so the two round-trip.

Also read: the list sections `*Arcslist` / `*Edgeslist` (`1 2 3 4` is the three
ties 1→2, 1→3, 1→4; no weights), the two-mode header `*Vertices n n₁` (a
two-mode network with `bipartite = n₁`), and self-loops (the network is
created with `loops=true` when the file contains one). A tie the network
cannot hold — an endpoint outside `1:n`, a within-mode pair under a two-mode
header — is an `ArgumentError`, never a silently smaller network; a tie listed
twice is kept once (with the last weight), with a warning, because parallel
edges are not supported. `*Matrix` sections are refused with an
`ArgumentError`; other `*` sections (`*Network`, `*Partition`, …) are skipped.

# Example
```julia
path = tempname() * ".net"
write(path, join(["*Vertices 3", "1 a", "2 b", "3 c", "*Arcs", "1 2 1", "2 3 2"], "\\n"))
net = read_pajek(path)
nv(net), ne(net), is_directed(net)      # (3, 2, true)
get_edge_attribute(net, :weight, 2, 3)  # 2.0
```
"""
function read_pajek(filepath::String)
    lines = readlines(filepath)

    n = 0
    n1 = nothing          # two-mode header `*Vertices n n1`
    saw_arcs = false
    saw_edges = false
    # (src, dst, weight, from_edges_section)
    edge_records = Tuple{Int, Int, Float64, Bool}[]
    has_weights = false  # any line carried an explicit weight column
    vertex_names = Dict{Int, String}()

    mode = :none
    for line in lines
        line = strip(line)
        isempty(line) && continue
        startswith(line, "%") && continue  # Comment

        lower_line = lowercase(line)

        if startswith(lower_line, "*")
            # Section headers. The `*arcslist`/`*edgeslist` checks must come
            # before `*arcs`/`*edges`, whose prefixes they share (the line
            # `1 2 3 4` of an arcs list used to be read as the arc 1→2 of
            # weight 3).
            if startswith(lower_line, "*vertices")
                mode = :vertices
                parts = split(line)
                length(parts) >= 2 && (n = parse(Int, parts[2]))
                length(parts) >= 3 && (n1 = parse(Int, parts[3]))
            elseif startswith(lower_line, "*arcslist")
                mode = :arcslist; saw_arcs = true
            elseif startswith(lower_line, "*edgeslist")
                mode = :edgeslist; saw_edges = true
            elseif startswith(lower_line, "*arcs")
                mode = :arcs; saw_arcs = true
            elseif startswith(lower_line, "*edges")
                mode = :edges; saw_edges = true
            elseif startswith(lower_line, "*matrix")
                throw(ArgumentError("read_pajek: $filepath has a *Matrix section, which " *
                                    "is not supported; convert it to *Arcs/*Edges"))
            else
                mode = :skip   # *Network title, *Partition, *Vector, ...
            end
        elseif mode == :vertices
            # Parse vertex line: id "label" or just id
            parts = split(line)
            if !isempty(parts)
                id = parse(Int, parts[1])
                if length(parts) >= 2
                    # Extract label (may be quoted)
                    label = join(parts[2:end], " ")
                    label = strip(label, ['"', '\''])
                    vertex_names[id] = label
                end
            end
        elseif mode in (:arcs, :edges)
            # Parse edge line: source target [weight]
            parts = split(line)
            if length(parts) >= 2
                src = parse(Int, parts[1])
                dst = parse(Int, parts[2])
                weight = 1.0
                if length(parts) >= 3
                    weight = parse(Float64, parts[3])
                    has_weights = true
                end
                push!(edge_records, (src, dst, weight, mode == :edges))
            end
        elseif mode in (:arcslist, :edgeslist)
            # Parse list line: source target1 target2 ...
            parts = split(line)
            src = parse(Int, parts[1])
            for t in parts[2:end]
                push!(edge_records, (src, parse(Int, t), 1.0, mode == :edgeslist))
            end
        end
    end

    # A file with any *Arcs section is directed; *Edges lines in a mixed
    # file are added as reciprocal arc pairs
    directed = saw_arcs
    if saw_arcs && saw_edges
        @warn "Pajek file mixes *Arcs and *Edges; reading as directed with *Edges lines as reciprocal arcs"
    end
    loops = any(r -> r[1] == r[2], edge_records)

    # Create network
    net = Network(n; directed=directed, loops=loops, bipartite=n1)

    # Add vertex names
    if !isempty(vertex_names)
        set_vertex_attribute!(net, :vertex_names, vertex_names)
    end

    # Add edges. The :weight attribute is stored whenever the file carried
    # an explicit weight column (even if every weight equals 1.0). Every tie
    # must land: a `false` from `add_edge!` is either a repeat (warned) or a
    # tie the network cannot hold (refused).
    duplicates = 0
    function place!(i, j, weight)
        (1 <= i <= n && 1 <= j <= n) ||
            throw(ArgumentError("read_pajek: tie ($i, $j) in $filepath references a " *
                                "vertex outside 1:$n"))
        if !add_edge!(net, i, j)
            has_edge(net, i, j) ||
                throw(ArgumentError("read_pajek: tie ($i, $j) in $filepath joins two " *
                                    "vertices of the same mode of a two-mode network " *
                                    "(*Vertices $n $(something(n1, ""))); only cross-mode " *
                                    "ties are allowed"))
            duplicates += 1
        end
        has_weights && set_edge_attribute!(net, :weight, i, j, weight)
    end
    for (src, dst, weight, from_edges) in edge_records
        place!(src, dst, weight)
        if directed && from_edges && src != dst
            place!(dst, src, weight)
        end
    end
    duplicates > 0 &&
        @warn "read_pajek: $duplicates tie(s) in $filepath were listed more than once; " *
              "parallel edges are not supported, so each was kept once" *
              (has_weights ? " (with the last weight listed)" : "")

    return net
end

"""
    write_pajek(net::Network, filepath::String; vertex_names::Symbol=:vertex_names,
                edge_weight::Symbol=:weight)

Write a network to a Pajek (.net) file.

# Arguments
- `net::Network`: The network to write
- `filepath::String`: Output file path
- `vertex_names::Symbol`: Vertex attribute to use as labels
- `edge_weight::Symbol`: Edge attribute to use as weights

# Example
```julia
net = network(3; directed=false)
add_edges!(net, [(1, 2), (2, 3)])
set_vertex_attribute!(net, :vertex_names, ["a", "b", "c"])
path = tempname() * ".net"
write_pajek(net, path)
ne(read_pajek(path))    # 2
```
"""
function write_pajek(net::Network, filepath::String;
                     vertex_names::Symbol=:vertex_names,
                     edge_weight::Symbol=:weight)
    open(filepath, "w") do io
        n = nv(net)

        # Write vertices (`*Vertices n n1` for a two-mode network, which
        # read_pajek reads back as `bipartite = n1`)
        println(io, is_two_mode(net) ? "*Vertices $n $(net.bipartite)" : "*Vertices $n")
        names = get_vertex_attribute(net, vertex_names)
        for v in 1:n
            if haskey(names, v)
                println(io, "$v \"$(names[v])\"")
            else
                println(io, v)
            end
        end

        # Write edges
        if net.directed
            println(io, "*Arcs")
        else
            println(io, "*Edges")
        end

        weights = get_edge_attribute(net, edge_weight)
        for e in edges(net)
            i, j = src(e), dst(e)
            edge = _canonical_edge(net, i, j)
            if haskey(weights, edge)
                println(io, "$i $j $(weights[edge])")
            else
                println(io, "$i $j")
            end
        end
    end
end

# ============================================================================
# GraphML Format (basic support)
# ============================================================================

"""
    write_graphml(net::Network, filepath::String)

Write a network to GraphML format (basic implementation): every vertex and
edge attribute becomes a `<key>`/`<data>` pair, and `edgedefault` follows the
network's directedness. There is no reader; the file is for igraph, Gephi,
NetworkX and the like.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_vertex_attribute!(net, :name, ["a", "b", "c"])
path = tempname() * ".graphml"
write_graphml(net, path)
occursin("<graphml", read(path, String))    # true
```
"""
function write_graphml(net::Network, filepath::String)
    open(filepath, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">")

        # Define attribute keys
        for attr in list_vertex_attributes(net)
            println(io, "  <key id=\"v_$attr\" for=\"node\" attr.name=\"$attr\" attr.type=\"string\"/>")
        end
        for attr in list_edge_attributes(net)
            println(io, "  <key id=\"e_$attr\" for=\"edge\" attr.name=\"$attr\" attr.type=\"string\"/>")
        end

        # Graph element
        edgedefault = net.directed ? "directed" : "undirected"
        println(io, "  <graph id=\"G\" edgedefault=\"$edgedefault\">")

        # Vertices
        for v in vertices(net)
            println(io, "    <node id=\"n$v\">")
            for attr in list_vertex_attributes(net)
                val = get_vertex_attribute(net, attr, v)
                if !isnothing(val)
                    println(io, "      <data key=\"v_$attr\">$val</data>")
                end
            end
            println(io, "    </node>")
        end

        # Edges
        edge_id = 0
        for e in edges(net)
            i, j = src(e), dst(e)
            println(io, "    <edge id=\"e$edge_id\" source=\"n$i\" target=\"n$j\">")
            for attr in list_edge_attributes(net)
                val = get_edge_attribute(net, attr, i, j)
                if !isnothing(val)
                    println(io, "      <data key=\"e_$attr\">$val</data>")
                end
            end
            println(io, "    </edge>")
            edge_id += 1
        end

        println(io, "  </graph>")
        println(io, "</graphml>")
    end
end

# ============================================================================
# Edge List CSV
# ============================================================================

# RFC-4180-style CSV field: quote when the value contains a comma, quote,
# or newline, doubling any inner quotes
function _csv_field(v)
    s = ismissing(v) ? "" : string(v)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

"""
    write_edgelist_csv(net::Network, filepath::String; attrs::Vector{Symbol}=Symbol[])

Write network as a CSV edge list with columns `source,target` plus one
column per requested edge attribute. Fields containing commas, quotes, or
newlines are quoted RFC-4180 style.

# Example
```julia
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, 1, 2, 2.5)
path = tempname() * ".csv"
write_edgelist_csv(net, path; attrs=[:weight])
readlines(path)    # ["source,target,weight", "1,2,2.5", "2,3,"]
```
"""
function write_edgelist_csv(net::Network{T}, filepath::String; attrs::Vector{Symbol}=Symbol[]) where T
    edge_list = as_edgelist(net)
    n_edges = size(edge_list, 1)

    open(filepath, "w") do io
        println(io, join(vcat(["source", "target"], string.(attrs)), ","))
        for r in 1:n_edges
            i, j = edge_list[r, 1], edge_list[r, 2]
            fields = Any[i, j]
            for attr in attrs
                push!(fields, get_edge_attribute(net, attr, i, j))
            end
            println(io, join((_csv_field(f === nothing ? missing : f) for f in fields), ","))
        end
    end
end
