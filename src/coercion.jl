# Coercion and conversion functions for Network objects.
#
# Provides functions to convert networks to/from matrices, edge lists,
# and DataFrames, similar to R's as.matrix.network, as.edgelist, etc.
# (The DataFrame conversions live in the NetworkDataFramesExt package
# extension and require DataFrames.jl to be loaded.)

# ============================================================================
# Network to Matrix
# ============================================================================

"""
    as_matrix(net::Network; attr=nothing, sparse=false, missing_value=nothing,
              return_mask=false, expand_bipartite=false)

Convert network to an adjacency matrix.

# Two-mode networks

For a two-mode network (`network(n; bipartite=n₁)` or a `BipartiteNetwork`)
the default is the `n_mode1 × n_mode2` two-mode **incidence matrix** (rows =
mode-1 vertices, columns = mode-2 vertices), matching R's `as.matrix.network`
for bipartite networks. Pass `expand_bipartite=true` (R's `expand.bipartite=`)
to get the full `(n₁ + n₂)`-square **adjacency matrix** instead — the mode-1
vertices in rows/columns `1:n₁`, the mode-2 vertices after them, and the
incidence block sitting at `A[1:n₁, n₁+1:end]` (mirrored into the lower block
for an undirected network). That square matrix is what every one-mode matrix
routine (eigenvector centrality, structural equivalence, QAP, ...) expects,
and it round-trips through `network_from_matrix`. `attr=`, `sparse=`,
`missing_value=` and `return_mask=` all apply to either shape.

On a **directed** two-mode network the incidence matrix can only represent
arcs from mode 1 to mode 2 (`B[i, j]` is the arc `i → k + j`). An arc in the
other direction (mode 2 → mode 1), and a mask on such a dyad, has no cell:
exactly as R's `as.matrix.network` does, it is **dropped** from the incidence
matrix (never written into the mirror cell as if it pointed the other way),
and NetworkCore.jl additionally emits a `@warn` naming how many were dropped.
Pass `expand_bipartite=true` for the square matrix that keeps every arc. On
an undirected two-mode network the question does not arise.

# Validation

`attr=` must name an existing edge attribute: `as_matrix(net; attr=:wt)` on a
network whose only edge attribute is `:weight` throws an `ArgumentError`
listing the attributes that exist, instead of returning the all-zero matrix
that would feed silently into QAP or an eigenvector centrality downstream (R's
`as.sociomatrix` warns "There is no edge attribute named wt"). An edge that
exists but carries no value for `attr` contributes `0`.

# Missing dyads

By default the matrix carries only the *face value* of each dyad, so a masked
(unobserved) dyad is indistinguishable from an absent tie and the mask is lost
on conversion. Two keywords make the conversion lossless:

- `missing_value=missing` marks each masked dyad with `missing`, giving an
  element type of `Union{Missing,Float64}`. [`network_from_matrix`](@ref) reads
  such a matrix back into a mask, so `network_from_matrix(as_matrix(net;
  missing_value=missing))` round-trips. Any other value (e.g. `NaN`, `-1`) is
  written literally instead.
- `return_mask=true` returns `(A, mask)` where `mask[i,j]` is `true` for the
  masked dyads, leaving `A` numeric.

# Arguments
- `net::Network`: The network to convert
- `attr::Union{Symbol,Nothing}=nothing`: Edge attribute to use as values (nothing = binary)
- `sparse::Bool=false`: Return a sparse matrix instead of dense
- `missing_value=nothing`: value written at masked dyads (`nothing` = face value)
- `return_mask::Bool=false`: also return the `BitMatrix` of masked dyads
- `expand_bipartite::Bool=false`: on a two-mode network, return the
  `(n₁ + n₂)`-square adjacency matrix instead of the `n₁ × n₂` incidence matrix

# Example
```julia
net = network(4)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, 1, 2, 2.0)
set_missing_dyad!(net, 3, 4)                 # 3→4 is unobserved
A = as_matrix(net)                           # face values: A[3,4] == 0
W = as_matrix(net; attr=:weight)             # W[1,2] == 2.0, W[2,3] == 0.0
M = as_matrix(net; missing_value=missing)    # lossless: M[3,4] === missing
A, mask = as_matrix(net; return_mask=true)   # mask[3,4] == true

b = network(5; bipartite=2, directed=false)
add_edge!(b, 1, 3)
size(as_matrix(b))                           # (2, 3) — incidence matrix
size(as_matrix(b; expand_bipartite=true))    # (5, 5) — square adjacency matrix

d = network(5; bipartite=2, directed=true)
add_edges!(d, [(1, 3), (4, 2)])
as_matrix(d)                                 # 2×3 with B[1,1] = 1 only; warns that 4→2 was dropped
as_matrix(d; expand_bipartite=true)[4, 2]    # 1.0 — the square matrix keeps it
```
"""
function as_matrix(net::Network{T}; attr::Union{Symbol,Nothing}=nothing, sparse::Bool=false,
                   missing_value=nothing, return_mask::Bool=false,
                   expand_bipartite::Bool=false) where T
    _check_edge_attr(net, attr, "as_matrix")
    if missing_value !== nothing || return_mask
        return _as_matrix_missing_aware(net; attr=attr, sparse=sparse,
                                        missing_value=missing_value,
                                        return_mask=return_mask,
                                        expand_bipartite=expand_bipartite)
    end
    if !isnothing(net.bipartite) && !expand_bipartite
        return _as_incidence_matrix(net; attr=attr, sparse=sparse)
    end

    n = nv(net)

    if sparse
        I = Int[]
        J = Int[]
        V = Float64[]

        for e in Graphs.edges(net.graph)
            i, j = src(e), dst(e)
            push!(I, i)
            push!(J, j)
            if isnothing(attr)
                push!(V, 1.0)
            else
                edge = _canonical_edge(net, T(i), T(j))
                val = get(get(net.edge_attrs, attr, Dict()), edge, 0.0)
                push!(V, Float64(val))
            end
        end

        return SparseArrays.sparse(I, J, V, n, n)
    else
        A = zeros(Float64, n, n)

        for e in Graphs.edges(net.graph)
            i, j = src(e), dst(e)
            if isnothing(attr)
                A[i, j] = 1.0
            else
                edge = _canonical_edge(net, T(i), T(j))
                A[i, j] = Float64(get(get(net.edge_attrs, attr, Dict()), edge, 0.0))
            end
        end

        return A
    end
end

# A mistyped edge-attribute name used to yield an all-zero matrix (or a column
# of `missing` from `as_edgelist`) with no diagnostic — the "missing attribute"
# common mistake, refused here where it is made. R's `as.sociomatrix` warns
# "There is no edge attribute named ...".
function _check_edge_attr(net::Network, attr, fname::AbstractString)
    attr === nothing && return nothing
    haskey(net.edge_attrs, attr) ||
        throw(ArgumentError("$fname: the network has no edge attribute :$attr; " *
                            "existing: $(sort!(collect(keys(net.edge_attrs))))"))
    return nothing
end

# A directed two-mode network can carry arcs from mode 2 to mode 1 (and masks
# on such dyads) that the n₁ × n₂ incidence matrix has no cell for. R's
# `as.matrix.network` drops them silently; the conversion contract says a drop
# must at least be reported, so say how many went and where the lossless
# shape is. Never mirror them into the mode-1 → mode-2 cell: that reports an
# arc that points the other way.
function _warn_dropped_reverse(net::Network, n_arcs::Int, n_masks::Int)
    (n_arcs > 0 || n_masks > 0) || return nothing
    k = net.bipartite::Int
    parts = String[]
    n_arcs > 0 && push!(parts, "$n_arcs arc(s) from mode 2 to mode 1")
    n_masks > 0 && push!(parts, "$n_masks masked dyad(s) from mode 2 to mode 1")
    @warn "as_matrix: the incidence matrix of a directed two-mode network holds only " *
          "mode-1 → mode-2 arcs (vertices 1:$k → $(k + 1):$(nv(net))); dropped " *
          join(parts, " and ") * ", as R's as.matrix.network does. Pass " *
          "expand_bipartite=true for the square adjacency matrix that keeps them."
    return nothing
end

# Two-mode incidence matrix for a bipartite network: rows are mode-1
# vertices 1:k, columns are mode-2 vertices (k+1):n mapped to 1:(n-k).
# Undirected: an edge is stored canonically with `src <= dst`, so `src` is the
# mode-1 endpoint already; the swap below is defensive. Directed: an arc with
# `src > k` points from mode 2 to mode 1 and has no cell — dropped and warned.
function _as_incidence_matrix(net::Network{T}; attr::Union{Symbol,Nothing}=nothing,
                              sparse::Bool=false) where T
    k = net.bipartite::Int
    n2 = nv(net) - k
    directed = is_directed(net)
    dropped = 0

    value(i, j) = isnothing(attr) ? 1.0 :
        Float64(get(get(net.edge_attrs, attr, Dict()),
                    _canonical_edge(net, T(i), T(j)), 0.0))

    if sparse
        I = Int[]
        J = Int[]
        V = Float64[]
        for e in edges(net)
            i, j = src(e), dst(e)
            if i > k
                directed && (dropped += 1; continue)
                (i, j) = (j, i)          # undirected: normalize to the mode-1 endpoint
            end
            push!(I, i)
            push!(J, j - k)
            push!(V, value(i, j))
        end
        _warn_dropped_reverse(net, dropped, 0)
        return SparseArrays.sparse(I, J, V, k, n2)
    else
        B = zeros(Float64, k, n2)
        for e in edges(net)
            i, j = src(e), dst(e)
            if i > k
                directed && (dropped += 1; continue)
                (i, j) = (j, i)
            end
            B[i, j - k] = value(i, j)
        end
        _warn_dropped_reverse(net, dropped, 0)
        return B
    end
end

as_matrix(net::BipartiteNetwork; kwargs...) = as_matrix(net.network; kwargs...)

# Missing-aware matrix conversion. Builds the ordinary face-value matrix first
# (so the two paths cannot drift), then overwrites the masked dyads.
#
# The mask is stored canonically — `(i,j)` for directed, `minmax(i,j)` for
# undirected — while the matrix is symmetric for undirected networks, so an
# undirected masked dyad has to be written to both triangles.
function _as_matrix_missing_aware(net::Network{T}; attr::Union{Symbol,Nothing},
                                  sparse::Bool, missing_value, return_mask::Bool,
                                  expand_bipartite::Bool=false) where T
    if sparse && missing_value !== nothing
        throw(ArgumentError(
            "as_matrix(...; sparse=true, missing_value=...) is not supported: a " *
            "sparse matrix cannot represent masked dyads distinctly from " *
            "structural zeros. Use `return_mask=true` (which returns a numeric " *
            "sparse matrix plus a separate mask) or `sparse=false`."))
    end

    A = as_matrix(net; attr=attr, sparse=sparse, expand_bipartite=expand_bipartite)
    # Only the INCIDENCE shape needs the two-mode index remap; the expanded
    # square adjacency of a two-mode network is indexed like a one-mode one.
    bip = !isnothing(net.bipartite) && !expand_bipartite
    k = bip ? net.bipartite::Int : 0

    # Map a canonical mask key onto its matrix index/indices. Bipartite
    # networks return an incidence matrix, so the mode-1 endpoint is the row
    # and the mode-2 endpoint is shifted into column space. A mask on a
    # mode-2 → mode-1 dyad of a DIRECTED two-mode network has no cell (like
    # the arc itself): dropped and warned, never reported as a mask on the
    # reverse dyad.
    directed = is_directed(net)
    function indices(i::T, j::T)
        if bip
            if i > k
                directed && return ()
                (i, j) = (j, i)
            end
            return ((Int(i), Int(j) - k),)
        elseif directed
            return ((Int(i), Int(j)),)
        else
            return i == j ? ((Int(i), Int(i)),) :
                            ((Int(i), Int(j)), (Int(j), Int(i)))
        end
    end

    mask = falses(size(A)...)
    dropped_masks = 0
    for (i, j) in net.missing_dyads
        idx = indices(i, j)
        isempty(idx) && (dropped_masks += 1)
        for (r, c) in idx
            mask[r, c] = true
        end
    end
    _warn_dropped_reverse(net, 0, dropped_masks)

    if missing_value !== nothing
        M = convert(Matrix{Union{typeof(missing_value),eltype(A)}}, A)
        for idx in eachindex(mask)
            mask[idx] && (M[idx] = missing_value)
        end
        A = M
    end

    return return_mask ? (A, mask) : A
end

"""
    as_adjacency_matrix(net::Network; kwargs...) -> Matrix

Alias for [`as_matrix`](@ref), taking the same keywords (`attr=`, `sparse=`,
`missing_value=`, `return_mask=`, `expand_bipartite=`). Matches R's
`as.matrix.network` with `matrix.type="adjacency"`.

# Example
```julia
net = network(3; directed=false)
add_edge!(net, 1, 2)
as_adjacency_matrix(net) == as_matrix(net)    # true
as_adjacency_matrix(net)[1, 2]                # 1.0
```
"""
as_adjacency_matrix(net::AbstractNetwork; kwargs...) = as_matrix(net; kwargs...)

# ============================================================================
# Network to Graphs.jl types
# ============================================================================

"""
    SimpleDiGraph(net::Network) -> SimpleDiGraph

Convert a network to a Graphs.jl `SimpleDiGraph` (a copy of the internal
storage). For undirected networks the result contains both arc directions.
"""
Graphs.SimpleDiGraph(net::Network) = copy(net.graph)

"""
    SimpleGraph(net::Network) -> SimpleGraph

Convert a network to an undirected Graphs.jl `SimpleGraph`. Directed edges
are collapsed to undirected ones.
"""
Graphs.SimpleGraph(net::Network) = Graphs.SimpleGraph(net.graph)

# ============================================================================
# Network to Edge List
# ============================================================================

"""
    as_edgelist(net::Network; attrs::Vector{Symbol}=Symbol[]) -> Matrix

Convert network to an edge list (n_edges x 2 matrix), one row per edge in
`edges(net)` order (each undirected edge once, `src <= dst`), matching R's
`as.edgelist`.

# Arguments
- `net::Network`: The network to convert
- `attrs::Vector{Symbol}`: Edge attributes to include as additional columns.
  Every name must be an existing edge attribute (`ArgumentError` listing the
  attributes that exist otherwise — a mistyped name used to return a column
  of `missing` silently); an edge that carries no value for an attribute
  gets `missing`.

# Returns
If no attrs: Matrix{Int} of size (n_edges, 2)
If attrs: Returns a tuple (edgelist::Matrix{Int}, attr_values::Dict{Symbol, Vector})

# Example
```julia
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, 1, 2, 2.5)
as_edgelist(net)                         # [1 2; 2 3]
el, vals = as_edgelist(net; attrs=[:weight])
vals[:weight]                            # Any[2.5, missing]
try as_edgelist(net; attrs=[:wt]) catch err; err end   # ArgumentError: no edge attribute :wt
```
"""
function as_edgelist(net::Network{T}; attrs::Vector{Symbol}=Symbol[]) where T
    for attr in attrs
        _check_edge_attr(net, attr, "as_edgelist")
    end
    edge_list = Tuple{T, T}[]

    for e in edges(net)
        push!(edge_list, (src(e), dst(e)))
    end

    n_edges = length(edge_list)

    if isempty(attrs)
        result = Matrix{T}(undef, n_edges, 2)
        for (idx, (i, j)) in enumerate(edge_list)
            result[idx, 1] = i
            result[idx, 2] = j
        end
        return result
    else
        # Include attributes
        result = Matrix{T}(undef, n_edges, 2)
        attr_vals = Dict{Symbol, Vector}()

        for attr in attrs
            attr_vals[attr] = Vector{Any}(undef, n_edges)
        end

        for (idx, (i, j)) in enumerate(edge_list)
            result[idx, 1] = i
            result[idx, 2] = j
            edge = _canonical_edge(net, i, j)
            for attr in attrs
                attr_vals[attr][idx] = get(get(net.edge_attrs, attr, Dict()), edge, missing)
            end
        end

        return (result, attr_vals)
    end
end

# ============================================================================
# Network to DataFrame (implemented in the NetworkDataFramesExt extension)
# ============================================================================

"""
    as_dataframe(net::Network; vertices::Bool=false, attrs::Vector{Symbol}=Symbol[]) -> DataFrame

Convert network to a DataFrame.

Requires DataFrames.jl to be loaded (`using DataFrames`); the method is
provided by the `NetworkDataFramesExt` package extension.

# Arguments
- `net::Network`: The network to convert
- `vertices::Bool=false`: If true, return vertex DataFrame instead of edge DataFrame
- `attrs::Vector{Symbol}`: Attributes to include (defaults to all)

# Example
```julia
using NetworkCore, DataFrames
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, 1, 2, 2.5)
set_vertex_attribute!(net, :name, Dict(1 => "a", 2 => "b", 3 => "c"))

edge_df = as_dataframe(net)                    # columns source, target, weight
edge_df.weight                                 # [2.5, missing]
vertex_df = as_dataframe(net; vertices=true)   # columns vertex, name
vertex_df.name                                 # ["a", "b", "c"]
```
"""
function as_dataframe end

# ============================================================================
# Matrix to Network
# ============================================================================

"""
    network_from_matrix(A::AbstractMatrix; directed::Bool=true, loops::Bool=false,
                        ignore_eval=nothing, names_eval=nothing,
                        missing_mask=nothing) -> Network

Create a network from an adjacency matrix.

# Missing dyads

Masked (unobserved) dyads can be restored in two ways, which differ in whether
the dyad's *face value* survives:

- **`missing` entries in `A`** are read back as masked dyads. Because the
  `missing` overwrote the face value in the matrix, the restored dyad has an
  *absent* face value. Mask-preserving, face-value-lossy — enough for
  `n_missing_dyads(network_from_matrix(as_matrix(net; missing_value=missing)))`
  to round-trip.
- **`missing_mask`** — a `Bool` matrix the size of `A` — masks those dyads while
  reading the face value from `A` itself. This is the fully lossless path, and
  pairs with `as_matrix(net; return_mask=true)`:

```julia
net = network(3)
add_edge!(net, 1, 2)
set_missing_dyad!(net, 2, 3)
A, mask = as_matrix(net; return_mask=true)
net2 = network_from_matrix(A; missing_mask=mask, directed=is_directed(net))
n_missing_dyads(net2)    # 1
```

# Arguments
- `A::AbstractMatrix`: Adjacency matrix (may contain `missing`)
- `directed::Bool=true`: Whether to create a directed network
- `loops::Bool=false`: Whether to allow self-loops
- `ignore_eval`: Function to test if a value should be ignored (default: ==(0))
- `names_eval`: Optional vertex names
- `missing_mask`: Optional `Bool` matrix marking masked dyads

# Example
```julia
A = [0 1 1; 1 0 0; 1 0 0]
net = network_from_matrix(A; directed=false)
```
"""
function network_from_matrix(A::AbstractMatrix;
                             directed::Bool=true,
                             loops::Bool=false,
                             ignore_eval=nothing,
                             names_eval=nothing,
                             attr::Symbol=:weight,
                             store_values::Union{Nothing,Bool}=nothing,
                             missing_mask::Union{Nothing,AbstractMatrix}=nothing)
    n = size(A, 1)
    size(A, 1) == size(A, 2) || throw(ArgumentError("Matrix must be square"))

    if !isnothing(missing_mask)
        size(missing_mask) == size(A) || throw(ArgumentError(
            "missing_mask must have the same size as A " *
            "(got $(size(missing_mask)) vs $(size(A)))"))
    end

    # Default: ignore zeros. `missing` is never an edge — it is an unobserved
    # dyad, recorded in the mask below rather than as a tie.
    base_ignore = isnothing(ignore_eval) ? (x -> x == 0) : ignore_eval
    should_ignore = x -> ismissing(x) || base_ignore(x)

    # For undirected networks only the upper triangle is read; warn if the
    # input is asymmetric so silent symmetrization doesn't go unnoticed
    # (isequal, so `missing` entries compare equal to themselves)
    if !directed && !isequal(A, transpose(A))
        @warn "Asymmetric matrix passed with directed=false; only the upper triangle is used"
    end

    net = Network(n; directed=directed, loops=loops)

    # Add vertex names if provided
    if !isnothing(names_eval)
        set_vertex_attribute!(net, :vertex_names, names_eval)
    end

    # Store matrix values as an edge attribute? Default: only when the
    # matrix is not purely 0/1. Pass store_values=true to force storage of
    # a genuinely weighted matrix whose values happen to all be 0/1.
    # (`missing` entries are dyad metadata, not values, so they don't count.)
    keep_values = isnothing(store_values) ?
        !all(x -> ismissing(x) || x == 0 || x == 1, A) : store_values

    for i in 1:n
        j_start = directed ? 1 : i  # For undirected, only upper triangle
        for j in j_start:n
            val = A[i, j]
            if !should_ignore(val)
                if !loops && i == j
                    continue
                end
                add_edge!(net, i, j)

                if keep_values
                    set_edge_attribute!(net, attr, i, j, val)
                end
            end
        end
    end

    # Restore the missing-dyad mask. Both sources are honoured: `missing`
    # entries in A (face value already lost) and an explicit mask matrix
    # (face value read from A above, so nothing is lost).
    for i in 1:n
        j_start = directed ? 1 : i
        for j in j_start:n
            (!loops && i == j) && continue
            masked = ismissing(A[i, j]) ||
                     (!isnothing(missing_mask) && Bool(missing_mask[i, j]))
            masked && set_missing_dyad!(net, i, j)
        end
    end

    return net
end

# ============================================================================
# Edge List to Network
# ============================================================================

"""
    network_from_edgelist(edges; n::Union{Int,Nothing}=nothing,
                          directed::Bool=true, loops::Bool=false) -> Network

Create a network from an edge list.

Every listed edge ends up in the network or the call throws: an endpoint
outside `1:n` (an undersized `n=`, or a mistyped ID) and a self-loop when
`loops=false` are `ArgumentError`s naming the offending edge, as R's
`network(el, matrix.type="edgelist")` errors on an illegal edgelist reference
— never a silently smaller network. Repeated edges collapse to one (parallel
edges are not supported).

# Arguments
- `edges`: Edge list as Matrix (n x 2) or Vector of tuples
- `n::Union{Int,Nothing}`: Number of vertices (auto-detected as the largest
  endpoint if nothing; pass it explicitly to keep trailing isolates)
- `directed::Bool=true`: Whether to create a directed network
- `loops::Bool=false`: Whether to allow self-loops

# Example
```julia
edges = [(1, 2), (2, 3), (3, 1)]
net = network_from_edgelist(edges)
```
"""
function network_from_edgelist(edges;
                               n::Union{Int,Nothing}=nothing,
                               directed::Bool=true,
                               loops::Bool=false)
    # Convert to vector of tuples if matrix
    if edges isa AbstractMatrix
        edge_tuples = [(edges[i, 1], edges[i, 2]) for i in 1:size(edges, 1)]
    else
        edge_tuples = collect(edges)
    end

    # Auto-detect n if not provided. Note auto-detection cannot see isolated
    # vertices with IDs beyond the largest endpoint; pass n explicitly then.
    if isnothing(n)
        n = isempty(edge_tuples) ? 0 :
            maximum(max(Int(i), Int(j)) for (i, j) in edge_tuples)
    end

    net = Network(n; directed=directed, loops=loops)

    for (i, j) in edge_tuples
        (1 <= i <= n && 1 <= j <= n) ||
            throw(ArgumentError("network_from_edgelist: edge ($i, $j) references a " *
                                "vertex outside 1:$n; pass a larger n= or fix the edge list"))
        (loops || i != j) ||
            throw(ArgumentError("network_from_edgelist: edge ($i, $j) is a self-loop " *
                                "but loops=false; pass loops=true or drop it"))
        add_edge!(net, i, j)
    end

    return net
end

"""
    network_from_dataframe(df::DataFrame; source::Symbol=:source, target::Symbol=:target,
                           directed::Bool=true, loops::Bool=false) -> Network

Create a network from a DataFrame with source and target columns.
Columns other than `source`/`target` become edge attributes.

Requires DataFrames.jl to be loaded (`using DataFrames`); the method is
provided by the `NetworkDataFramesExt` package extension.

# Arguments
- `df::DataFrame`: DataFrame with edge information
- `source::Symbol`: Column name for source vertices
- `target::Symbol`: Column name for target vertices
- `directed::Bool`: Whether to create a directed network
- `loops::Bool`: Whether to allow self-loops

# Example
```julia
using NetworkCore, DataFrames
df = DataFrame(source = [1, 2, 3], target = [2, 3, 1], weight = [1.0, 2.0, 1.5])
net = network_from_dataframe(df)
ne(net), get_edge_attribute(net, :weight, 3, 1)   # (3, 1.5)
```
"""
function network_from_dataframe end
