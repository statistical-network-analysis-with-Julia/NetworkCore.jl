# Attribute handling for Network objects.
#
# Provides functions to get, set, list, and delete vertex, edge, and network-level
# attributes, similar to R's network package attribute functions.

# ============================================================================
# Vertex Attributes
# ============================================================================

"""
    get_vertex_attribute(net::Network, attr::Symbol) -> Dict

Get a vertex attribute by name. Returns a Dict mapping vertex IDs to values.

# Example
```julia
net = load_dataset(:florentine_marriage)
names = get_vertex_attribute(net, :name)    # Dict{Int64, Any}: 1 => "Acciaiuoli", ...
names[9]                                    # "Medici"
```
"""
function get_vertex_attribute(net::Network{T}, attr::Symbol) where T
    # Lazy default: the eager `get(d, k, Dict())` built (and discarded) the
    # empty Dict on EVERY call, 80 B even when the attribute exists — paid per
    # sampler step by any term reading the whole-Dict getter (ERGMUserterms).
    return get(() -> Dict{T, Any}(), net.vertex_attrs, attr)
end

"""
    get_vertex_attribute(net::Network, attr::Symbol, v::Integer) -> Any

Get the attribute value for a specific vertex.
"""
function get_vertex_attribute(net::Network{T}, attr::Symbol, v::Integer) where T
    attrs = get(net.vertex_attrs, attr, nothing)
    if isnothing(attrs)
        return nothing
    end
    return get(attrs, T(v), nothing)
end

"""
    get_vertex_attribute(net::Network{T}, attr::Symbol, ::Type{V}) -> Dict{T,V}

Get a vertex attribute as a concretely typed `Dict{T,V}`, converting every
value to `V` in a single pass.

Attribute storage is untyped (`Dict{T,Any}`), which makes every read in a
tight loop type-unstable. This accessor is the fast path for model code:
convert once up front, then work with the concretely typed dictionary.
Throws if any stored value cannot be converted to `V`. Returns an empty
`Dict{T,V}` if the attribute does not exist.

See also [`vertex_attribute_vector`](@ref) for a dense `Vector{V}` in
vertex order.

# Example
```julia
net = load_dataset(:florentine_marriage)
wealth = get_vertex_attribute(net, :wealth, Float64)  # Dict{Int,Float64}
wealth[1]                                             # 10.0 (Acciaiuoli)
```
"""
function get_vertex_attribute(net::Network{T}, attr::Symbol, ::Type{V}) where {T, V}
    attrs = get(net.vertex_attrs, attr, nothing)
    out = Dict{T, V}()
    isnothing(attrs) && return out
    sizehint!(out, length(attrs))
    for (k, v) in attrs
        out[k] = convert(V, v)::V
    end
    return out
end

"""
    vertex_attribute_vector(net::Network, attr::Symbol, ::Type{V}; default) -> Vector{V}

Get a vertex attribute as a dense, concretely typed `Vector{V}` with one
entry per vertex (in vertex order `1:nv(net)`), converting every value to
`V` in a single pass.

Like the typed `get_vertex_attribute(net, attr, V)`, this is the fast path
for model code that reads attributes in hot loops: the returned vector is
concretely typed and indexable by vertex ID. Vertices without a stored
value take `default`; if no `default` is given and any vertex lacks a
value (or the attribute does not exist), an `ArgumentError` naming the
attribute and the vertex is thrown — and, when the attribute does not exist
at all, listing the attributes the network does have (so a mistyped name is
caught where it is made rather than three calls later).

# Example
```julia
net = load_dataset(:florentine_marriage)
wealth = vertex_attribute_vector(net, :wealth, Float64)             # 16-element Vector{Float64}
set_vertex_attribute!(net, :group, Dict(9 => "Medici"))
group  = vertex_attribute_vector(net, :group, String; default="?")  # gaps filled with "?"
try vertex_attribute_vector(net, :group, String) catch err; err end  # ArgumentError: not set for vertex 1
```
"""
function vertex_attribute_vector(net::Network{T}, attr::Symbol, ::Type{V};
                                 default=nothing) where {T, V}
    attrs = get(net.vertex_attrs, attr, nothing)
    n = nv(net)
    out = Vector{V}(undef, n)
    for i in Base.OneTo(n)
        val = isnothing(attrs) ? nothing : get(attrs, T(i), nothing)
        if val === nothing
            default === nothing && throw(_missing_vertex_value_error(net, attr, i, attrs))
            out[i] = convert(V, default)::V
        else
            out[i] = convert(V, val)::V
        end
    end
    return out
end

# The one message for a vertex-attribute value that is not there: names the
# attribute AND the vertex, lists the attributes that do exist when the name
# itself is unknown, and points at the fix. A bare `KeyError(i)` (the previous
# behaviour) was not actionable for the canonical mistyped-attribute mistake.
function _missing_vertex_value_error(net::Network, attr::Symbol, i, attrs)
    msg = "vertex attribute :$attr is not set for vertex $i"
    if attrs === nothing
        have = sort!(list_vertex_attributes(net))
        msg *= " (the network has no attribute :$attr; existing: $(have))"
    end
    return ArgumentError(msg * "; pass default= to fill gaps")
end

# Vertex IDs are validated on every attribute write, mirroring R's `network`
# ("vertex IDs out of range"): a typo'd ID would otherwise store a value for a
# vertex that does not exist and read back as plausible, wrong data.
# The message is only built on failure, so the setters stay cheap in the
# samplers that call them per proposal (ERGMCount's valued sampler).
@inline function _check_vertex(net::Network, v::Integer, fname::Symbol, attr::Symbol)
    1 <= v <= nv(net) ||
        throw(ArgumentError("$fname(net, :$attr, ...): vertex $v does not exist " *
                            "(network has $(nv(net)) vertices)"))
    return v
end

"""
    set_vertex_attribute!(net::Network, attr::Symbol, values::Dict)

Set a vertex attribute from a Dict mapping vertex IDs to values. Every key
must be an existing vertex (`ArgumentError` otherwise, nothing stored).

# Example
```julia
net = network(3)
set_vertex_attribute!(net, :name, Dict(1 => "Alice", 2 => "Bob"))
get_vertex_attribute(net, :name, 2)    # "Bob"
```
"""
function set_vertex_attribute!(net::Network{T}, attr::Symbol, values::Dict) where T
    for k in keys(values)
        _check_vertex(net, k, :set_vertex_attribute!, attr)
    end
    if !haskey(net.vertex_attrs, attr)
        net.vertex_attrs[attr] = Dict{T, Any}()
    end
    for (k, v) in values
        net.vertex_attrs[attr][T(k)] = v
    end
    return net
end

"""
    set_vertex_attribute!(net::Network, attr::Symbol, v::Integer, value)

Set the attribute value for a specific vertex. Throws an `ArgumentError` when
`v` is not a vertex of the network (`1 <= v <= nv(net)`), as R's
`set.vertex.attribute` does for an out-of-range ID, instead of storing a
value nobody can reach by a legal vertex ID.

# Example
```julia
net = network(3)
set_vertex_attribute!(net, :name, 1, "Alice")
get_vertex_attribute(net, :name, 1)          # "Alice"
try set_vertex_attribute!(net, :name, 9, "Zed") catch err; err end  # ArgumentError: vertex 9 does not exist
```
"""
function set_vertex_attribute!(net::Network{T}, attr::Symbol, v::Integer, value) where T
    _check_vertex(net, v, :set_vertex_attribute!, attr)
    if !haskey(net.vertex_attrs, attr)
        net.vertex_attrs[attr] = Dict{T, Any}()
    end
    net.vertex_attrs[attr][T(v)] = value
    return net
end

"""
    set_vertex_attribute!(net::Network, attr::Symbol, values::Vector)

Set a vertex attribute from a Vector (assumes values are in vertex order 1:n).
"""
function set_vertex_attribute!(net::Network{T}, attr::Symbol, values::Vector) where T
    length(values) == nv(net) || throw(ArgumentError("Vector length must match number of vertices"))
    net.vertex_attrs[attr] = Dict{T, Any}(T(i) => v for (i, v) in enumerate(values))
    return net
end

"""
    delete_vertex_attribute!(net::Network, attr::Symbol)

Delete a vertex attribute (every vertex's value). A no-op when the attribute
does not exist.

# Example
```julia
net = network(3)
set_vertex_attribute!(net, :name, ["a", "b", "c"])
delete_vertex_attribute!(net, :name)
list_vertex_attributes(net)    # Symbol[]
```
"""
function delete_vertex_attribute!(net::Network, attr::Symbol)
    delete!(net.vertex_attrs, attr)
    return net
end

"""
    list_vertex_attributes(net::Network) -> Vector{Symbol}

List all vertex attribute names (unsorted).

# Example
```julia
net = network(3)
set_vertex_attribute!(net, :age, [30, 40, 50])
set_vertex_attribute!(net, :name, ["a", "b", "c"])
sort(list_vertex_attributes(net))    # [:age, :name]
```
"""
function list_vertex_attributes(net::Network)
    return collect(keys(net.vertex_attrs))
end

# ============================================================================
# Edge Attributes
# ============================================================================

"""
    get_edge_attribute(net::Network, attr::Symbol) -> Dict

Get an edge attribute by name. Returns a Dict mapping (source, target) tuples to values.
"""
function get_edge_attribute(net::Network{T}, attr::Symbol) where T
    # Lazy default, as in `get_vertex_attribute`: a hit allocates nothing.
    return get(() -> Dict{Tuple{T,T}, Any}(), net.edge_attrs, attr)
end

"""
    get_edge_attribute(net::Network, attr::Symbol, i::Integer, j::Integer) -> Any

Get the attribute value for a specific edge.
"""
function get_edge_attribute(net::Network{T}, attr::Symbol, i::Integer, j::Integer) where T
    attrs = get(net.edge_attrs, attr, nothing)
    if isnothing(attrs)
        return nothing
    end
    edge = _canonical_edge(net, T(i), T(j))
    return get(attrs, edge, nothing)
end

"""
    get_edge_attribute(net::Network{T}, attr::Symbol, ::Type{V}) -> Dict{Tuple{T,T},V}

Get an edge attribute as a concretely typed `Dict{Tuple{T,T},V}`, converting
every value to `V` in a single pass. Keys are canonical edges (`(i,j)` as
stored; `minmax(i,j)` for undirected networks).

This is the fast path for model code reading edge attributes in hot loops;
see the typed `get_vertex_attribute(net, attr, V)` for rationale. Throws if
any stored value cannot be converted to `V`. Returns an empty dict if the
attribute does not exist.

# Example
```julia
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, Dict((1, 2) => 1, (2, 3) => 2.5))
weights = get_edge_attribute(net, :weight, Float64)  # Dict{Tuple{Int,Int},Float64}
weights[(1, 2)]                                       # 1.0
```
"""
function get_edge_attribute(net::Network{T}, attr::Symbol, ::Type{V}) where {T, V}
    attrs = get(net.edge_attrs, attr, nothing)
    out = Dict{Tuple{T,T}, V}()
    isnothing(attrs) && return out
    sizehint!(out, length(attrs))
    for (k, v) in attrs
        out[k] = convert(V, v)::V
    end
    return out
end

"""
    set_edge_attribute!(net::Network, attr::Symbol, values::Dict; require_edge::Bool=true)

Set an edge attribute from a Dict mapping (source, target) tuples to values.

Validates exactly as the single-edge form does, and **all-or-nothing**: every
key is checked before anything is written. Both endpoints of every key must be
vertices of the network, and — by default — every key must be an existing
edge; an `ArgumentError` names the first offending key otherwise and the
attribute is left untouched. Pass `require_edge=false` to store values ahead
of the edges deliberately (the same opt-out as the single-edge form).

# Example
```julia
net = network(3)
add_edges!(net, [(1, 2), (2, 3)])
set_edge_attribute!(net, :weight, Dict((1,2) => 1.5, (2,3) => 2.0))
get_edge_attribute(net, :weight, 2, 3)                          # 2.0
try set_edge_attribute!(net, :w, Dict((1,2) => 1.0, (1,3) => 2.0)) catch err; err end  # ArgumentError: (1, 3) is not an edge
:w in list_edge_attributes(net)                                 # false — nothing was stored
set_edge_attribute!(net, :w, Dict((1,3) => 2.0); require_edge=false)   # stored ahead of the edge
```
"""
function set_edge_attribute!(net::Network{T}, attr::Symbol, values::Dict;
                             require_edge::Bool=true) where T
    for (i, j) in keys(values)
        _check_vertex(net, i, :set_edge_attribute!, attr)
        _check_vertex(net, j, :set_edge_attribute!, attr)
        require_edge && !has_edge(net, i, j) && throw(_not_an_edge_error(attr, i, j))
    end
    if !haskey(net.edge_attrs, attr)
        net.edge_attrs[attr] = Dict{Tuple{T,T}, Any}()
    end
    for ((i, j), v) in values
        edge = _canonical_edge(net, T(i), T(j))
        net.edge_attrs[attr][edge] = v
    end
    return net
end

# The one message for an edge attribute addressed at a dyad that carries no
# tie, shared by the single-edge and Dict forms of `set_edge_attribute!`.
_not_an_edge_error(attr::Symbol, i, j) =
    ArgumentError("set_edge_attribute!(net, :$attr, ...): ($i, $j) is not an " *
                  "edge of the network; add_edge!(net, $i, $j) first, or pass " *
                  "require_edge=false to store the value ahead of the edge")

"""
    set_edge_attribute!(net::Network, attr::Symbol, i::Integer, j::Integer, value;
                        require_edge::Bool=true)

Set the attribute value for a specific edge.

Both endpoints must be vertices of the network, and — by default — `(i, j)`
must be an existing edge: an `ArgumentError` is thrown otherwise, matching R's
`set.edge.attribute`, which can only address an edge that exists. Storing a
value on a dyad that carries no tie is almost always a typo'd endpoint or a
forgotten `add_edge!`, and the value would otherwise sit unreachable until an
edge happened to appear there. Pass `require_edge=false` to store the value
ahead of the edge deliberately (the attribute is then reported by
`get_edge_attribute(net, attr)` and applies as soon as the edge is added).

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 2.5)
get_edge_attribute(net, :weight, 1, 2)        # 2.5
try set_edge_attribute!(net, :weight, 2, 3, 1.0) catch err; err end  # ArgumentError: (2, 3) is not an edge
set_edge_attribute!(net, :weight, 2, 3, 1.0; require_edge=false)   # stored ahead
```
"""
function set_edge_attribute!(net::Network{T}, attr::Symbol, i::Integer, j::Integer, value;
                             require_edge::Bool=true) where T
    _check_vertex(net, i, :set_edge_attribute!, attr)
    _check_vertex(net, j, :set_edge_attribute!, attr)
    require_edge && !has_edge(net, i, j) && throw(_not_an_edge_error(attr, i, j))
    if !haskey(net.edge_attrs, attr)
        net.edge_attrs[attr] = Dict{Tuple{T,T}, Any}()
    end
    edge = _canonical_edge(net, T(i), T(j))
    net.edge_attrs[attr][edge] = value
    return net
end

"""
    delete_edge_attribute!(net::Network, attr::Symbol)

Delete an edge attribute (every edge's value). A no-op when the attribute
does not exist.

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 2.5)
delete_edge_attribute!(net, :weight)
get_edge_attribute(net, :weight, 1, 2) === nothing    # true
```
"""
function delete_edge_attribute!(net::Network, attr::Symbol)
    delete!(net.edge_attrs, attr)
    return net
end

"""
    list_edge_attributes(net::Network) -> Vector{Symbol}

List all edge attribute names (unsorted).

# Example
```julia
net = network(3)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 2.5)
list_edge_attributes(net)    # [:weight]
```
"""
function list_edge_attributes(net::Network)
    return collect(keys(net.edge_attrs))
end

# ============================================================================
# Network Attributes
# ============================================================================

"""
    get_network_attribute(net::Network, attr::Symbol) -> Any

Get a network-level attribute by name; `nothing` when it is not set. The
Julia counterpart of R's `net %n% "attr"`.

# Example
```julia
net = network(3)
set_network_attribute!(net, :title, "Advice network")
get_network_attribute(net, :title)     # "Advice network"
get_network_attribute(net, :year)      # nothing
```
"""
function get_network_attribute(net::Network, attr::Symbol)
    return get(net.network_attrs, attr, nothing)
end

"""
    set_network_attribute!(net::Network, attr::Symbol, value)

Set a network-level attribute.

# Example
```julia
net = network(3)
set_network_attribute!(net, :title, "Friendship Network")
get_network_attribute(net, :title)    # "Friendship Network"
```
"""
function set_network_attribute!(net::Network, attr::Symbol, value)
    net.network_attrs[attr] = value
    return net
end

"""
    delete_network_attribute!(net::Network, attr::Symbol)

Delete a network-level attribute. A no-op when the attribute does not exist.

# Example
```julia
net = network(3)
set_network_attribute!(net, :title, "Advice network")
delete_network_attribute!(net, :title)
get_network_attribute(net, :title)    # nothing
```
"""
function delete_network_attribute!(net::Network, attr::Symbol)
    delete!(net.network_attrs, attr)
    return net
end

"""
    list_network_attributes(net::Network) -> Vector{Symbol}

List all network-level attribute names (unsorted).

# Example
```julia
net = network(3)
set_network_attribute!(net, :title, "Advice network")
list_network_attributes(net)    # [:title]
```
"""
function list_network_attributes(net::Network)
    return collect(keys(net.network_attrs))
end

# ============================================================================
# Convenience indexing (similar to R's %v% and %e% operators)
# ============================================================================

# net[:v, :attr] for vertex attributes
# net[:e, :attr] for edge attributes
# net[:n, :attr] for network attributes

function Base.getindex(net::Network, ::Val{:v}, attr::Symbol)
    return get_vertex_attribute(net, attr)
end

function Base.getindex(net::Network, ::Val{:e}, attr::Symbol)
    return get_edge_attribute(net, attr)
end

function Base.getindex(net::Network, ::Val{:n}, attr::Symbol)
    return get_network_attribute(net, attr)
end

function Base.setindex!(net::Network, value, ::Val{:v}, attr::Symbol)
    set_vertex_attribute!(net, attr, value)
end

function Base.setindex!(net::Network, value, ::Val{:e}, attr::Symbol)
    set_edge_attribute!(net, attr, value)
end

function Base.setindex!(net::Network, value, ::Val{:n}, attr::Symbol)
    set_network_attribute!(net, attr, value)
end

# Shorthand using symbols: net[:v, :name], net[:e, :weight], net[:n, :title]
# (only :v, :e, and :n have Val methods)
Base.getindex(net::Network, level::Symbol, attr::Symbol) = Base.getindex(net, Val(level), attr)
Base.setindex!(net::Network, value, level::Symbol, attr::Symbol) = Base.setindex!(net, value, Val(level), attr)
