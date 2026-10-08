# DataFrames.jl integration for NetworkCore.jl, loaded automatically when both
# Network and DataFrames are in the environment (package extension).
module NetworkDataFramesExt

using NetworkCore
using DataFrames
import NetworkCore: as_dataframe, network_from_dataframe, _canonical_edge

function as_dataframe(net::Network{T}; vertices::Bool=false,
                      attrs::Vector{Symbol}=Symbol[]) where T
    if vertices
        # Vertex DataFrame
        n = nv(net)
        df = DataFrame(vertex = collect(1:n))

        attrs_to_include = isempty(attrs) ? list_vertex_attributes(net) : attrs
        for attr in attrs_to_include
            vals = get_vertex_attribute(net, attr)
            df[!, attr] = [get(vals, T(i), missing) for i in 1:n]
        end

        return df
    else
        # Edge DataFrame
        edge_list = as_edgelist(net)
        n_edges = size(edge_list, 1)

        df = DataFrame(
            source = edge_list[:, 1],
            target = edge_list[:, 2]
        )

        attrs_to_include = isempty(attrs) ? list_edge_attributes(net) : attrs
        for attr in attrs_to_include
            edge_vals = get_edge_attribute(net, attr)
            df[!, attr] = [get(edge_vals, _canonical_edge(net, edge_list[i, 1], edge_list[i, 2]), missing)
                          for i in 1:n_edges]
        end

        return df
    end
end

function network_from_dataframe(df::DataFrame;
                                source::Symbol=:source,
                                target::Symbol=:target,
                                directed::Bool=true,
                                loops::Bool=false)
    edge_pairs = [(df[i, source], df[i, target]) for i in 1:nrow(df)]
    n = isempty(edge_pairs) ? 0 : maximum(max(e[1], e[2]) for e in edge_pairs)

    net = Network(n; directed=directed, loops=loops)

    # Add edges
    for (i, j) in edge_pairs
        add_edge!(net, i, j)
    end

    # Add edge attributes from other columns
    other_cols = setdiff(names(df), [string(source), string(target)])
    for col in other_cols
        attr = Symbol(col)
        for (idx, (i, j)) in enumerate(edge_pairs)
            set_edge_attribute!(net, attr, i, j, df[idx, col])
        end
    end

    return net
end

end # module
