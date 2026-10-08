using NetworkCore
using Test
using TOML
using Graphs
using DataFrames
using Random
using LinearAlgebra: dot, norm, eigvals, Symmetric, diag
import LinearAlgebra
using Distributions: Normal, ccdf, cdf   # reference for z_pvalues only
using Aqua
using Graphs: complement, induced_subgraph, adjacency_matrix
import StatsAPI   # the StatsAPI verbs the dummy fits below implement


# Dummy fitted-model types for the "CoefficientTable / coeftable" and
# "check_statsapi" testsets (a struct cannot be declared inside a @testset).
# `DummyFullFit` implements all ten verbs of the StatsAPI surface on its OWN
# type, the way a model package does; `DummyPartialFit` implements three.
struct DummyFullFit
    θ::Vector{Float64}
    V::Matrix{Float64}
end
StatsAPI.coef(f::DummyFullFit) = f.θ
StatsAPI.vcov(f::DummyFullFit) = f.V
StatsAPI.stderror(f::DummyFullFit) = sqrt.(diag(f.V))
StatsAPI.confint(f::DummyFullFit; level::Real=0.95) =
    hcat(f.θ .- 1.96 .* StatsAPI.stderror(f), f.θ .+ 1.96 .* StatsAPI.stderror(f))
StatsAPI.loglikelihood(::DummyFullFit) = -12.5
StatsAPI.nobs(::DummyFullFit) = 120
StatsAPI.dof(f::DummyFullFit) = length(f.θ)
StatsAPI.aic(f::DummyFullFit) = -2 * StatsAPI.loglikelihood(f) + 2 * StatsAPI.dof(f)
StatsAPI.bic(f::DummyFullFit) =
    -2 * StatsAPI.loglikelihood(f) + log(StatsAPI.nobs(f)) * StatsAPI.dof(f)
StatsAPI.coeftable(f::DummyFullFit) =
    CoefficientTable(["a", "b"][1:min(2, length(f.θ))], f.θ[1:min(2, length(f.θ))],
                     StatsAPI.stderror(f)[1:min(2, length(f.θ))])
Base.show(io::IO, f::DummyFullFit) = show(io, StatsAPI.coeftable(f))

struct DummyPartialFit
    θ::Vector{Float64}
    se::Vector{Float64}
end
StatsAPI.coef(f::DummyPartialFit) = f.θ
StatsAPI.stderror(f::DummyPartialFit) = f.se
StatsAPI.coeftable(f::DummyPartialFit) =
    CoefficientTable(["a", "b"][1:length(f.θ)], f.θ, f.se[1:length(f.θ)])

struct DummyDataFrameCoeftable end
StatsAPI.coeftable(::DummyDataFrameCoeftable) = Dict(:a => 1.0)   # not a CoefficientTable

# A function with a recognisable name for the "Shared task helpers" testset: its
# frame must appear in the backtrace `spawn_all` keeps for a failing task.
@noinline _spawn_all_inner_fail(c) = throw(ArgumentError("inner failure $c"))

struct DummyThrowingFit end
StatsAPI.coef(::DummyThrowingFit) = error("coefficients not available")

# `coefnames`, the eleventh verb: one fit whose names agree with its table, one
# whose names disagree, one with the wrong count and one returning symbols.
StatsAPI.coefnames(f::DummyFullFit) = ["a", "b"][1:min(2, length(f.θ))]
struct DummyNamedFit
    θ::Vector{Float64}
    names::AbstractVector
end
StatsAPI.coef(f::DummyNamedFit) = f.θ
StatsAPI.coefnames(f::DummyNamedFit) = f.names
StatsAPI.coeftable(f::DummyNamedFit) =
    CoefficientTable(["a", "b"], f.θ, fill(0.1, length(f.θ)))

# A 1×1 "Hessian" whose `\` raises a caller-side error, for the "Shared Newton
# optimizer" testset: newton_fit must let it propagate, not report it as
# non-convergence (only a failed factorisation is "no Newton direction").
struct ExplodingHessian <: AbstractMatrix{Float64} end
Base.size(::ExplodingHessian) = (1, 1)
Base.getindex(::ExplodingHessian, ::Int, ::Int) = -1.0
Base.:-(h::ExplodingHessian) = h
Base.:\(::ExplodingHessian, ::AbstractVector) = error("boom from the caller's code")

# ---------------------------------------------------------------------------
# Brute-force oracle for the separation verdict, independent of its linear
# programme: a design with margin rows M is separated iff the cone
# {d : M d ≥ 0} contains a d with M d ≠ 0. Modulo its lineality space (the
# null space of M) that cone is pointed, so it is non-trivial iff it has an
# extreme ray, and every extreme ray is the one-dimensional null space of
# r − 1 independent rows of M plus the lineality basis (r = rank M).
# Enumerating every such subset, in exact rational arithmetic, decides the
# question exhaustively for small designs.
# ---------------------------------------------------------------------------
const OracleQ = Rational{BigInt}

function oracle_rref(A::Matrix{OracleQ})
    A = copy(A); m, n = size(A); piv = Int[]; r = 0
    for j in 1:n
        r == m && break
        i = findfirst(k -> !iszero(A[k, j]), r+1:m)
        i === nothing && continue
        i += r
        A[[r + 1, i], :] = A[[i, r + 1], :]
        r += 1
        A[r, :] ./= A[r, j]
        for k in 1:m
            k != r && !iszero(A[k, j]) && (A[k, :] .-= A[k, j] .* A[r, :])
        end
        push!(piv, j)
    end
    return A[1:r, :], piv
end

function oracle_nullspace(A::Matrix{OracleQ}, n::Int)
    size(A, 1) == 0 && return Matrix{OracleQ}(LinearAlgebra.I, n, n)
    R, piv = oracle_rref(A)
    free = setdiff(1:n, piv)
    N = zeros(OracleQ, n, length(free))
    for (c, f) in enumerate(free)
        N[f, c] = 1
        for (r, p) in enumerate(piv)
            N[p, c] = -R[r, f]
        end
    end
    return N
end

oracle_combinations(v, k) = k == 0 ? [Int[]] :
    [vcat(x, rest) for (i, x) in enumerate(v) for rest in oracle_combinations(v[i+1:end], k - 1)]

function oracle_separated(M::AbstractMatrix{<:Real})
    m, q = size(M)
    m == 0 && return false
    MQ = OracleQ.(M)
    L = oracle_nullspace(MQ, q)
    r = q - size(L, 2)
    r == 0 && return false
    for S in oracle_combinations(collect(1:m), r - 1)
        N = oracle_nullspace(vcat(MQ[S, :], permutedims(L)), q)
        size(N, 2) == 1 || continue
        for s in (1, -1)
            v = MQ * (s .* N[:, 1])
            all(>=(0), v) && any(>(0), v) && return true
        end
    end
    return false
end

# The margin rows of each likelihood family, written from its definition
# (every case against EVERY alternative for the conditional logit; no
# reduction), for the oracle.
_rows(blocks, q) = isempty(blocks) ? zeros(Int, 0, q) : reduce(vcat, blocks)
oracle_logit_rows(X, y) = _rows([(y[i] ? 1 : -1) .* X[i:i, :] for i in axes(X, 1)], size(X, 2))
oracle_binomial_rows(X, nt, n1) = _rows(vcat(
    [X[r:r, :] for r in axes(X, 1) if nt[r] > 0 && n1[r] > 0],
    [-X[r:r, :] for r in axes(X, 1) if nt[r] > 0 && n1[r] < nt[r]]), size(X, 2))
oracle_poisson_rows(X, y) = _rows(vcat([-X[i:i, :] for i in axes(X, 1)],
                                       [X[i:i, :] for i in axes(X, 1) if y[i] > 0]), size(X, 2))
oracle_clogit_rows(X, ch, st) = _rows(
    [X[c:c, :] .- X[j:j, :] for c in eachindex(st) if ch[c] for j in eachindex(st)
     if st[j] == st[c] && j != c], size(X, 2))

# The ```julia blocks of NetworkCore's docstrings, collected by the docstring
# testset and executed by the last testset of the file.
const DOC_BLOCKS = Tuple{Symbol,String}[]

@testset "NetworkCore.jl" begin
    @testset "Network Construction" begin
        # Basic construction
        net = network(5)
        @test nv(net) == 5
        @test ne(net) == 0
        @test is_directed(net) == true

        # Undirected network
        net_undir = network(5; directed=false)
        @test is_directed(net_undir) == false

        # Using network() function
        net2 = network(10; directed=false)
        @test nv(net2) == 10
        @test is_directed(net2) == false

        # network_initialize alias
        net3 = network_initialize(3)
        @test nv(net3) == 3
    end

    @testset "Edge Operations" begin
        net = network(5)

        # Add edges
        @test add_edge!(net, 1, 2) == true
        @test add_edge!(net, 2, 3) == true
        @test add_edge!(net, 3, 4) == true

        @test ne(net) == 3
        @test has_edge(net, 1, 2) == true
        @test has_edge(net, 2, 1) == false  # Directed

        # Remove edge
        @test rem_edge!(net, 1, 2) == true
        @test has_edge(net, 1, 2) == false
        @test ne(net) == 2

        # Self-loops (not allowed by default)
        @test add_edge!(net, 1, 1) == false

        # Self-loops (when allowed)
        net_loops = network(3; loops=true)
        @test add_edge!(net_loops, 1, 1) == true
        @test has_edge(net_loops, 1, 1) == true
    end

    @testset "Undirected Edge Operations" begin
        net = network(5; directed=false)

        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)

        # Both directions should work
        @test has_edge(net, 1, 2) == true
        @test has_edge(net, 2, 1) == true

        # Edge count should be 2 (not 4)
        @test ne(net) == 2
    end

    @testset "Vertex Attributes" begin
        net = network(3)

        # Set via Dict
        set_vertex_attribute!(net, :name, Dict(1 => "Alice", 2 => "Bob", 3 => "Carol"))

        @test get_vertex_attribute(net, :name, 1) == "Alice"
        @test get_vertex_attribute(net, :name, 2) == "Bob"

        # Set via Vector
        set_vertex_attribute!(net, :age, [25, 30, 35])
        @test get_vertex_attribute(net, :age, 1) == 25
        @test get_vertex_attribute(net, :age, 3) == 35

        # Set single vertex
        set_vertex_attribute!(net, :score, 1, 100.0)
        @test get_vertex_attribute(net, :score, 1) == 100.0

        # List attributes
        @test :name in list_vertex_attributes(net)
        @test :age in list_vertex_attributes(net)

        # Delete attribute
        delete_vertex_attribute!(net, :score)
        @test !(:score in list_vertex_attributes(net))
    end

    @testset "Edge Attributes" begin
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)

        # Set via Dict
        set_edge_attribute!(net, :weight, Dict((1,2) => 1.5, (2,3) => 2.0))

        @test get_edge_attribute(net, :weight, 1, 2) == 1.5
        @test get_edge_attribute(net, :weight, 2, 3) == 2.0

        # Set single edge
        set_edge_attribute!(net, :label, 1, 2, "friendship")
        @test get_edge_attribute(net, :label, 1, 2) == "friendship"

        # List attributes
        @test :weight in list_edge_attributes(net)
    end

    @testset "Network Attributes" begin
        net = network(3)

        set_network_attribute!(net, :title, "Test Network")
        set_network_attribute!(net, :year, 2024)

        @test get_network_attribute(net, :title) == "Test Network"
        @test get_network_attribute(net, :year) == 2024

        @test :title in list_network_attributes(net)

        delete_network_attribute!(net, :year)
        @test isnothing(get_network_attribute(net, :year))
    end

    @testset "Neighbors" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 4, 1)

        @test Set(outneighbors(net, 1)) == Set([2, 3])
        @test Set(inneighbors(net, 1)) == Set([4])
        @test Set(neighbors(net, 1)) == Set([2, 3])  # outneighbors for directed
    end

    @testset "Matrix Conversion" begin
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 1)

        A = as_matrix(net)
        @test size(A) == (3, 3)
        @test A[1, 2] == 1.0
        @test A[2, 3] == 1.0
        @test A[3, 1] == 1.0
        @test A[1, 1] == 0.0

        # Weighted matrix
        set_edge_attribute!(net, :weight, Dict((1,2) => 0.5, (2,3) => 1.0, (3,1) => 0.75))
        W = as_matrix(net; attr=:weight)
        @test W[1, 2] == 0.5
        @test W[2, 3] == 1.0
    end

    @testset "Edge List Conversion" begin
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)

        el = as_edgelist(net)
        @test size(el) == (2, 2)
        @test (1, 2) in [(el[i, 1], el[i, 2]) for i in 1:2]
        @test (2, 3) in [(el[i, 1], el[i, 2]) for i in 1:2]
    end

    @testset "Network from Matrix" begin
        A = [0 1 0; 0 0 1; 1 0 0]
        net = network_from_matrix(A)

        @test nv(net) == 3
        @test ne(net) == 3
        @test has_edge(net, 1, 2)
        @test has_edge(net, 2, 3)
        @test has_edge(net, 3, 1)

        # Undirected from symmetric matrix
        B = [0 1 1; 1 0 0; 1 0 0]
        net_undir = network_from_matrix(B; directed=false)
        @test ne(net_undir) == 2  # (1,2) and (1,3)
    end

    @testset "Network from Edge List" begin
        edges = [(1, 2), (2, 3), (3, 1)]
        net = network_from_edgelist(edges)

        @test nv(net) == 3
        @test ne(net) == 3
        @test has_edge(net, 1, 2)
        @test has_edge(net, 2, 3)
        @test has_edge(net, 3, 1)
    end

    @testset "Common mistakes are actionable: out-of-range IDs, non-edges, undersized n=" begin
        # R's `network` errors on all of these ("vertex IDs out of range",
        # illegal edgelist reference). Until 2026-09 NetworkCore.jl stored the
        # value or dropped the edge silently and handed back plausible, wrong
        # output — the failure class the whole ecosystem was reorganised to stop.
        flo = load_dataset(:florentine_marriage)
        @test nv(flo) == 16

        # --- vertex ID that does not exist: refused, nothing stored ---
        err = try set_vertex_attribute!(flo, :x, 99, 1.0); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("vertex 99 does not exist", err.msg)
        @test occursin("16 vertices", err.msg)
        @test occursin(":x", err.msg)                         # names the attribute
        @test !(:x in list_vertex_attributes(flo))
        @test_throws ArgumentError set_vertex_attribute!(flo, :x, 0, 1.0)
        # The Dict form validates every key before writing any
        @test_throws ArgumentError set_vertex_attribute!(flo, :x, Dict(1 => 1.0, 99 => 2.0))
        @test !(:x in list_vertex_attributes(flo))
        set_vertex_attribute!(flo, :x, 16, 1.0)               # the last vertex is a vertex
        @test get_vertex_attribute(flo, :x, 16) == 1.0

        # --- edge attribute on a dyad that carries no tie ---
        @test !has_edge(flo, 1, 2)
        err = try set_edge_attribute!(flo, :w, 1, 2, 1.0); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("(1, 2) is not an edge", err.msg)
        @test occursin("add_edge!(net, 1, 2)", err.msg)       # says what to do
        @test occursin("require_edge=false", err.msg)         # ... and how to opt out
        @test get_edge_attribute(flo, :w, 1, 2) === nothing
        @test !(:w in list_edge_attributes(flo))
        @test_throws ArgumentError set_edge_attribute!(flo, :w, 1, 99, 1.0)     # endpoint typo
        @test_throws ArgumentError set_edge_attribute!(flo, :w, 99, 1, 1.0)
        # The documented opt-out stores the value ahead of the edge, in writing
        set_edge_attribute!(flo, :w, 1, 2, 1.0; require_edge=false)
        @test get_edge_attribute(flo, :w, 1, 2) == 1.0
        add_edge!(flo, 1, 2)
        @test get_edge_attribute(flo, :w, 1, 2) == 1.0
        # An existing edge is unchanged, in either endpoint order (undirected)
        e = first(edges(flo))
        set_edge_attribute!(flo, :w, dst(e), src(e), 2.0)
        @test get_edge_attribute(flo, :w, src(e), dst(e)) == 2.0
        # Directed: the reverse dyad of an arc is not an edge
        d = network(3; directed=true)
        add_edge!(d, 1, 2)
        set_edge_attribute!(d, :w, 1, 2, 1.0)
        @test_throws ArgumentError set_edge_attribute!(d, :w, 2, 1, 1.0)

        # --- mistyped attribute name in the typed fast path ---
        err = try vertex_attribute_vector(flo, :welth, Float64); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("no attribute :welth", err.msg)
        @test occursin(":wealth", err.msg)                    # the name the user meant
        @test occursin("default=", err.msg)

        # --- undersized n= or a typo'd endpoint in an edge list ---
        err = try network_from_edgelist([(1, 2), (2, 5)]; n=3); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("(2, 5)", err.msg)
        @test occursin("outside 1:3", err.msg)
        @test occursin("n=", err.msg)
        @test_throws ArgumentError network_from_edgelist([1 2; 2 5]; n=3)   # matrix form too
        @test_throws ArgumentError network_from_edgelist([(0, 1)])
        # A self-loop with loops=false used to vanish silently
        err = try network_from_edgelist([(1, 2), (2, 2)]); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("(2, 2) is a self-loop", err.msg)
        @test occursin("loops=true", err.msg)
        l = network_from_edgelist([(1, 2), (2, 2)]; loops=true)
        @test ne(l) == 2 && has_edge(l, 2, 2)
        ok = network_from_edgelist([(1, 2), (2, 5)]; n=5)
        @test nv(ok) == 5 && ne(ok) == 2 && has_edge(ok, 2, 5)

        # --- add_edges!, the bulk adder every example uses, never drops an
        # edge either: R's add.edges "Illegal vertex reference"
        small = network(4)
        err = try add_edges!(small, [(1, 2), (1, 99)]); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("(1, 99)", err.msg)
        @test occursin("outside 1:4", err.msg)
        @test occursin("add_vertices!", err.msg) && occursin("n=", err.msg)
        @test ne(small) == 1                                   # edges before the bad one were added
        @test_throws ArgumentError add_edges!(small, [(0, 1)])
        err = try add_edges!(small, [(2, 2)]); nothing catch e; e end
        @test err isa ArgumentError && occursin("(2, 2) is a self-loop", err.msg) &&
              occursin("loops=true", err.msg)
        @test add_edges!(network(3; loops=true), [(2, 2)]) == 1
        err = try add_edges!(network(5; bipartite=2), [(1, 3), (1, 2)]); nothing catch e; e end
        @test err isa ArgumentError && occursin("(1, 2)", err.msg) && occursin("same mode", err.msg)
        @test add_edges!(small, [(1, 2), (3, 4), (1, 2)]) == 1  # a repeat counts 0, is not an error
        @test ne(small) == 2
        # The single add_edge! keeps Graphs.jl's Bool convention
        @test add_edge!(small, 1, 99) == false

        # --- the Dict form of set_edge_attribute! validates like the single-edge form
        d3 = network(4)
        add_edges!(d3, [(1, 2), (2, 3)])
        err = try set_edge_attribute!(d3, :w, Dict((1, 2) => 1.0, (3, 9) => 2.0)); nothing catch e; e end
        @test err isa ArgumentError && occursin("vertex 9 does not exist", err.msg)
        @test !(:w in list_edge_attributes(d3))                # all-or-nothing: nothing stored
        err = try set_edge_attribute!(d3, :w, Dict((2, 3) => 5.0, (1, 3) => 1.0)); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("(1, 3) is not an edge", err.msg) && occursin("require_edge=false", err.msg)
        @test !(:w in list_edge_attributes(d3))
        @test get_edge_attribute(d3, :w, 2, 3) === nothing     # the good key was not written either
        set_edge_attribute!(d3, :w, Dict((1, 3) => 1.0); require_edge=false)   # the documented opt-out
        @test get_edge_attribute(d3, :w, 1, 3) == 1.0
        set_edge_attribute!(d3, :w, Dict((1, 2) => 1.0, (2, 3) => 2.0))
        @test get_edge_attribute(d3, :w, 2, 3) == 2.0
        # require_edge=false still checks the vertex range
        @test_throws ArgumentError set_edge_attribute!(d3, :w, Dict((1, 42) => 1.0); require_edge=false)

        # --- mistyped edge-attribute name in as_matrix / as_edgelist:
        # used to be an all-zero matrix / a column of `missing`, silently
        set_edge_attribute!(d3, :weight, 1, 2, 2.5)
        err = try as_matrix(d3; attr=:wt); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("no edge attribute :wt", err.msg)
        @test occursin(":weight", err.msg)                     # lists what exists
        @test_throws ArgumentError as_matrix(d3; attr=:wt, sparse=true)
        @test_throws ArgumentError as_matrix(d3; attr=:wt, return_mask=true)
        @test_throws ArgumentError as_matrix(network(5; bipartite=2); attr=:wt)
        @test_throws ArgumentError as_edgelist(d3; attrs=[:weight, :wt])
        @test as_matrix(d3; attr=:weight)[1, 2] == 2.5
        @test as_matrix(d3; attr=:weight)[2, 3] == 0.0         # an edge without a value is 0
        el, vals = as_edgelist(d3; attrs=[:weight])
        @test vals[:weight][1] == 2.5 && ismissing(vals[:weight][2])

        # --- three refusals R makes
        err = try network(-1); nothing catch e; e end          # R: "cannot be of negative order"
        @test err isa ArgumentError && occursin("non-negative", err.msg) && occursin("-1", err.msg)
        @test_throws ArgumentError Network(-3)
        @test_throws ArgumentError Network{Int}(; n=-1)
        @test nv(network(0)) == 0
        err = try get_induced_subgraph(d3, [1, 99]); nothing catch e; e end   # R: "Illegal vertex selection"
        @test err isa ArgumentError && occursin("vertex 99", err.msg)
        err = try get_induced_subgraph(d3, [1, 2, 1]); nothing catch e; e end
        @test err isa ArgumentError && occursin("repeated vertex (1)", err.msg)
        @test_throws ArgumentError get_induced_subgraph(d3, [0, 1])
        @test nv(get_induced_subgraph(d3, [4, 1])) == 2
        err = try permute_vertices(d3, [1, 1, 2, 3]); nothing catch e; e end  # R: "Invalid permutation vector"
        @test err isa ArgumentError && occursin("permutation of 1:4", err.msg)
        @test_throws ArgumentError permute_vertices(d3, [1, 2, 3])          # wrong length
        @test_throws ArgumentError permute_vertices(d3, [0, 1, 2, 3])
        @test_throws ArgumentError permute_vertices(d3, [5, 1, 2, 3])
        @test nv(permute_vertices(d3, [4, 3, 2, 1])) == 4
    end

    @testset "Directed two-mode incidence matrix drops mode-2 → mode-1 arcs, as R" begin
        # R 4.6.1, network 1.20.0:
        #   b <- network.initialize(6, directed=TRUE, bipartite=2)
        #   b <- add.edges(b, tail=c(1,5), head=c(3,1)); b[6,2] <- NA
        #   as.matrix(b)                         # 2×4: only [1,1] == 1, no NA
        #   as.matrix(b, expand.bipartite=TRUE)  # [1,3] == 1, [5,1] == 1, [6,2] == NA
        # The arc 5→1 used to be written at B[1, col(5)] — reported as
        # the arc 1→5 — and the mask on 6→2 as a mask on 2→6. An incidence
        # matrix has no cell for a reverse arc; it is dropped and warned, never
        # mis-oriented.
        b = network(6; directed=true, bipartite=2)
        add_edges!(b, [(1, 3), (5, 1)])
        set_missing_dyad!(b, 6, 2)
        @test ne(b) == 2 && has_edge(b, 5, 1) && is_missing_dyad(b, 6, 2)

        B = @test_logs (:warn, r"dropped 1 arc\(s\) from mode 2 to mode 1") as_matrix(b)
        @test size(B) == (2, 4)
        @test B == [1 0 0 0; 0 0 0 0]                     # R's as.matrix(b)
        Bs = @test_logs (:warn, r"expand_bipartite=true") as_matrix(b; sparse=true)
        @test Matrix(Bs) == B
        # The mask: one warning for the arc (inner face-value pass), one for the mask
        A, mask = @test_logs (:warn, r"1 arc\(s\)") (:warn, r"1 masked dyad\(s\)") as_matrix(b; return_mask=true)
        @test A == B && !any(mask)                        # R: no NA in the incidence matrix
        M = @test_logs (:warn, r"arc") (:warn, r"masked") as_matrix(b; missing_value=missing)
        @test !any(ismissing, M)

        # The square matrix keeps everything, silently (nothing is dropped)
        E = @test_logs as_matrix(b; expand_bipartite=true)
        @test size(E) == (6, 6) && E[1, 3] == 1 && E[5, 1] == 1 && E[1, 5] == 0
        Em, mm = @test_logs as_matrix(b; expand_bipartite=true, return_mask=true)
        @test mm[6, 2] && !mm[2, 6] && count(mm) == 1

        # A directed two-mode network with only mode-1 → mode-2 arcs: no warning, same as before
        ok = network(6; directed=true, bipartite=2)
        add_edges!(ok, [(1, 3), (2, 6)]); set_missing_dyad!(ok, 1, 6)
        Bo, mo = @test_logs as_matrix(ok; return_mask=true)
        @test Bo == [1 0 0 0; 0 0 0 1] && mo == Bool[0 0 0 1; 0 0 0 0]

        # Undirected: an edge listed mode-2-first is the same edge; nothing to drop
        u = network(5; bipartite=2, directed=false)
        add_edge!(u, 4, 1); set_missing_dyad!(u, 5, 2)
        Bu, mu = @test_logs as_matrix(u; return_mask=true)
        @test Bu == [0 1 0; 0 0 0] && mu == Bool[0 0 0; 0 0 1]      # R: as.matrix(u)[1, "4"] == 1
    end

    @testset "show: one line inside containers, R-style block at top level" begin
        # `load_dataset(:s50)` — a NamedTuple of a Vector of networks — used to
        # print a garbled wall of multi-line blocks because the ONLY show
        # method was the 2-arg one. Base uses 2-arg show inside containers.
        net = network(5; directed=false)
        add_edges!(net, [(1, 2), (2, 3), (3, 4)])
        one = repr(net)
        @test one == "Network{Int64,false}(5 vertices, 3 edges)"
        @test !occursin('\n', one)
        @test string(net) == one && sprint(print, net) == one
        @test !occursin('\n', repr(Dict(:a => net)))
        @test !occursin('\n', repr([net, net]))
        s50 = load_dataset(:s50)
        @test !occursin("Vertices:", repr(s50.friendship))
        @test occursin("(50 vertices, 113 edges)", repr(s50))

        # The block: what R's print.network always states, sorted attributes
        set_vertex_attribute!(net, :zeta, collect(1:5)); set_vertex_attribute!(net, :alpha, collect(1:5))
        set_edge_attribute!(net, :weight, 1, 2, 1.0)
        set_network_attribute!(net, :title, "t")
        set_missing_dyad!(net, 4, 5)
        block = repr(MIME"text/plain"(), net)
        @test startswith(block, "Network{Int64}: undirected network\n")
        @test occursin("  Vertices: 5\n  Edges: 3\n", block)
        @test occursin("Directed: false", block)
        @test occursin("Loops: false", block)
        @test occursin("Vertex attributes: alpha, zeta", block)     # sorted, as R lists them
        @test occursin("Edge attributes: weight", block) && occursin("Network attributes: title", block)
        @test occursin("Missing (unobserved) dyads: 1", block)
        @test !occursin("Two-mode", block)
        @test sprint(show, MIME"text/plain"(), net) == block
        @test sprint(io -> show(io, MIME"text/plain"(), net)) == block   # what display() calls

        # loops and the two-mode size are visible (the two classic confusions)
        l = network(3; loops=true)
        @test occursin("Loops: true", repr(MIME"text/plain"(), l))
        @test occursin("loops", repr(l))
        @test repr(MIME"text/plain"(), network(3)) != repr(MIME"text/plain"(), l)
        tm = network(6; bipartite=2)
        tblock = repr(MIME"text/plain"(), tm)
        @test occursin("two-mode directed network", tblock)
        @test occursin("Two-mode: 2 mode-1 vertices (1:2), 4 mode-2 vertices (3:6)", tblock)
        @test occursin("two-mode 2+4", repr(tm))

        # BipartiteNetwork: the same pair, delegating to the wrapped network
        bn = BipartiteNetwork(2, 3; directed=false)
        add_edge!(bn, 1, 3)
        @test repr(bn) == "BipartiteNetwork{Int64}(2+3 vertices, 1 edges)"
        bblock = repr(MIME"text/plain"(), bn)
        @test startswith(bblock, "BipartiteNetwork{Int64}: undirected two-mode network\n")
        @test occursin("Vertices: 5", bblock) && occursin("Edges: 1", bblock)
        @test occursin("Two-mode: 2 mode-1 vertices (1:2), 3 mode-2 vertices (3:5)", bblock)
        @test occursin("Loops: false", bblock)
    end

    @testset "Every exported and public docstring carries a runnable example" begin
        # Every exported and public binding has a docstring with an example,
        # and every ```julia block in those docstrings RUNS. A docs build with
        # checkdocs=:exports checks presence, not content, and a block that
        # is only checked for presence can rot unnoticed (the `as_dataframe`
        # example referred to a `net` it never built). Each block runs in a
        # fresh module that has done `using NetworkCore` (Documenter's
        # DocTestSetup does the same); anything else it needs, it loads.
        meta = Base.Docs.meta(NetworkCore)
        missing_example = String[]
        undocumented = String[]
        blocks = DOC_BLOCKS
        for nm in names(NetworkCore)
            nm === :NetworkCore && continue
            b = Base.Docs.Binding(NetworkCore, nm)
            if !haskey(meta, b)
                push!(undocumented, string(nm))
                continue
            end
            has_example = false
            for (_, ds) in meta[b].docs
                txt = ds.text isa AbstractString ? ds.text : join(string.(ds.text), "\n")
                if occursin("```julia", txt) || occursin("```jldoctest", txt)
                    has_example = true
                end
                for m in eachmatch(r"```julia\n(.*?)```"s, txt)
                    push!(blocks, (nm, m.captures[1]))
                end
            end
            has_example || push!(missing_example, string(nm))
        end
        @test isempty(undocumented)
        @test isempty(missing_example)
        @test length(blocks) >= 120
        # the blocks are EXECUTED by the last testset of this file: some of
        # them add methods to the shared generics (`gof`, `compute`, …),
        # which the "empty generics" testsets below must see empty
        # The Graphs.jl re-exports are documented HERE (NetworkCore-owned docstrings
        # attached to the foreign bindings), so they are covered by the walk above
        for nm in (:degree, :indegree, :outdegree, :nv, :ne, :add_edge!)
            @test haskey(meta, Base.Docs.Binding(NetworkCore, nm))
        end
    end

    @testset "Network Density" begin
        # Complete directed graph
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 1)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 1)
        add_edge!(net, 3, 2)

        @test network_density(net) == 1.0  # 6 edges / 6 possible

        # Partial network
        net2 = network(3)
        add_edge!(net2, 1, 2)
        add_edge!(net2, 2, 3)
        @test network_density(net2) ≈ 2/6

        # --- Missing dyads: refused by default, face value only in writing ---
        masked = copy(net2)
        set_missing_dyad!(masked, 3, 1)
        @test_throws ArgumentError network_density(masked)
        msg = try
            network_density(masked)
        catch e
            sprint(showerror, e)
        end
        @test occursin("network_density", msg)   # names the caller
        @test occursin("missing=:face", msg)     # ...and the opt-in it DOES have
        @test network_density(masked; missing=:face) == network_density(net2)
        @test_throws ArgumentError network_density(masked; missing=:omit)  # R's na.omit is not offered
        @test supports_missing(network_density) == false   # no principled treatment, just the guard

        # --- Two-mode: the wrapper and the flagged Network agree, both ways ---
        bnet = BipartiteNetwork(2, 3; directed=false)
        flagged = network(5; bipartite=2, directed=false)
        for (i, j) in [(1, 3), (1, 4), (2, 5)]
            add_edge!(bnet, i, j)
            add_edge!(flagged, i, j)
        end
        @test network_density(bnet) == network_density(flagged) == 3 / 10   # R default: all dyads
        @test network_density(bnet; discount_bipartite=true) ==
              network_density(flagged; discount_bipartite=true) == 3 / 6     # n1·n2 cross-mode dyads
        # Directed two-mode: one dyad per direction, so the discounted denominator doubles
        bdir = BipartiteNetwork(2, 3; directed=true)
        fdir = network(5; bipartite=2, directed=true)
        for (i, j) in [(1, 3), (1, 4), (4, 2)]
            add_edge!(bdir, i, j)
            add_edge!(fdir, i, j)
        end
        @test network_density(bdir) == network_density(fdir) == 3 / 20
        @test network_density(bdir; discount_bipartite=true) ==
              network_density(fdir; discount_bipartite=true) == 3 / 12
        # The keyword is a no-op on a one-mode network
        @test network_density(net2; discount_bipartite=true) == network_density(net2)
        # The BipartiteNetwork method forwards explicit keywords, not a
        # `kwargs...` splat (which allocated 16–32 B per call): 0 B like the
        # Network method it wraps, for every keyword combination.
        dens_default(b) = network_density(b)
        dens_discount(b) = network_density(b; discount_bipartite=true)
        dens_all(b) = network_density(b; missing=:face, diag=false, discount_bipartite=true)
        for b in (bdir, bnet, fdir, flagged)
            dens_default(b); dens_discount(b); dens_all(b)          # warm up
            @test @allocated(dens_default(b)) == 0
            @test @allocated(dens_discount(b)) == 0
            @test @allocated(dens_all(b)) == 0
        end

        # --- diag: self-loops leave numerator AND denominator under diag=false ---
        lnet = network(3; loops=true)
        add_edge!(lnet, 1, 1)
        add_edge!(lnet, 1, 2)
        @test network_density(lnet) == 2 / 9                # default diag = net.loops = true
        @test network_density(lnet; diag=true) == 2 / 9
        @test network_density(lnet; diag=false) == 1 / 6    # sna::gden(diag=FALSE)
        # Four vertices, not three: on a 3-vertex example both keywords give
        # exactly 1/3 (2/6 and 1/3), so a `diag` that silently no-oped on the
        # undirected branch would have passed. Here the two values differ.
        ulnet = network(4; directed=false, loops=true)
        add_edge!(ulnet, 2, 2)
        add_edge!(ulnet, 1, 2)
        @test network_density(ulnet) == 2 / 10              # n(n+1)/2 = 10 dyads, loop counted
        @test network_density(ulnet; diag=true) == 2 / 10
        @test network_density(ulnet; diag=false) == 1 / 6   # n(n-1)/2 = 6, loop leaves both
        @test network_density(ulnet) != network_density(ulnet; diag=false)
        # diag=true on a network that cannot hold loops only widens the denominator
        @test network_density(net2; diag=true) == 2 / 9

        # The policy vocabulary is declared, not implied
        @test missing_policies(network_density) == (:error, :face)
    end

    @testset "network_density matches R network.density" begin
        # Frozen output of R `network::network.density` (test/fixtures/r/
        # network_density.R) on the same hand-built networks: every branch of
        # the denominator, two-mode with and without discount.bipartite, and an
        # NA edge counted at face value (na.omit=FALSE == missing=:face).
        g = load_golden(joinpath(@__DIR__, "fixtures", "network_density.toml"))
        @test g.provenance["network_version"] == "1.20.0"

        function build(n; directed, loops=false, bipartite=nothing, edges)
            net = network(n; directed=directed, loops=loops, bipartite=bipartite)
            add_edges!(net, edges)
            return net
        end
        check(key, value) = check_golden(g, key, value) || error(golden_report(g, key, value))

        dir_noloops   = build(4; directed=true,  edges=[(1, 2), (2, 3), (3, 4)])
        undir_noloops = build(4; directed=false, edges=[(1, 2), (2, 3), (3, 4)])
        dir_loops     = build(4; directed=true,  loops=true, edges=[(1, 2), (1, 1), (2, 3), (3, 1)])
        undir_loops   = build(4; directed=false, loops=true, edges=[(1, 2), (2, 3), (3, 3)])
        bip_undir     = build(5; directed=false, bipartite=2, edges=[(1, 3), (1, 4), (2, 5)])
        bip_dir       = build(5; directed=true,  bipartite=2, edges=[(1, 3), (1, 4), (4, 2)])

        @test check("directed_noloops",   network_density(dir_noloops))
        @test check("undirected_noloops", network_density(undir_noloops))
        @test check("directed_loops",     network_density(dir_loops))
        @test check("undirected_loops",   network_density(undir_loops))
        @test check("bipartite_undirected",          network_density(bip_undir))
        @test check("bipartite_undirected_discount", network_density(bip_undir; discount_bipartite=true))
        @test check("bipartite_directed",            network_density(bip_dir))
        @test check("bipartite_directed_discount",   network_density(bip_dir; discount_bipartite=true))

        # diag=false is pinned to sna::gden(diag=FALSE) (network.density has no
        # diag argument), on loops networks where the two values differ, so a
        # silent no-op on either branch would fail here; diag=true (the
        # default on a loops network) stays network.density.
        dir_loops_diag   = build(4; directed=true,  loops=true,
                                 edges=[(1, 2), (1, 1), (2, 3), (3, 1), (4, 2)])
        undir_loops_diag = build(4; directed=false, loops=true, edges=[(2, 2), (1, 2)])
        @test check("directed_loops_diag",     network_density(dir_loops_diag))
        @test check("directed_loops_nodiag",   network_density(dir_loops_diag; diag=false))
        @test check("undirected_loops_diag",   network_density(undir_loops_diag))
        @test check("undirected_loops_nodiag", network_density(undir_loops_diag; diag=false))
        @test check("undirected_loops_fixture_nodiag", network_density(undir_loops; diag=false))
        @test g.values["directed_loops_diag"] != g.values["directed_loops_nodiag"]
        @test g.values["undirected_loops_diag"] != g.values["undirected_loops_nodiag"]
        @test g.provenance["sna_version"] == "2.8"

        # The BipartiteNetwork wrapper is the same object to R (one network with
        # a %n% "bipartite" attribute), so it must give the same numbers.
        wb = BipartiteNetwork(2, 3; directed=false)
        foreach(((i, j),) -> add_edge!(wb, i, j), [(1, 3), (1, 4), (2, 5)])
        @test check("bipartite_undirected",          network_density(wb))
        @test check("bipartite_undirected_discount", network_density(wb; discount_bipartite=true))

        # One NA edge, na.omit=FALSE: the unobserved 2->3 tie counts at its stored
        # face value (present) over the full dyad count. That is `missing=:face`;
        # the default refuses, and R's na.omit=TRUE has no Julia counterpart.
        na_net = build(4; directed=true, edges=[(1, 2), (2, 3), (3, 4)])
        set_missing_dyad!(na_net, 2, 3)
        @test_throws ArgumentError network_density(na_net)
        @test check("na_face", network_density(na_net; missing=:face))

        @test g.values["single_vertex_is_nan"] == true
        @test isnan(network_density(network(1)))
    end

    @testset "Induced Subgraph" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        set_vertex_attribute!(net, :name, Dict(1 => "A", 2 => "B", 3 => "C", 4 => "D", 5 => "E"))

        sub = get_induced_subgraph(net, [2, 3, 4])

        @test nv(sub) == 3
        @test ne(sub) == 2  # (2,3) and (3,4) mapped to (1,2) and (2,3)
        @test has_edge(sub, 1, 2)  # was 2→3
        @test has_edge(sub, 2, 3)  # was 3→4
    end

    @testset "Neighborhood" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        # Order 1 neighborhood
        n1 = get_neighborhood(net, 1, 1)
        @test 1 in n1
        @test 2 in n1
        @test !(3 in n1)

        # Order 2 neighborhood
        n2 = get_neighborhood(net, 1, 2)
        @test 1 in n2
        @test 2 in n2
        @test 3 in n2
        @test !(4 in n2)
    end

    @testset "Bipartite Network" begin
        bnet = BipartiteNetwork(3, 4)

        @test nv(bnet) == 7
        @test is_bipartite(bnet) == true
        @test bnet.n_mode1 == 3
        @test bnet.n_mode2 == 4
    end

    @testset "Unsupported multigraph/hypergraph flags" begin
        @test_throws ArgumentError network(3; multiple=true)
        @test_throws ArgumentError network(3; hyper=true)
        @test_throws ArgumentError network(3; bipartite=5)
    end

    @testset "Undirected self-loops" begin
        net = network(3; directed=false, loops=true)
        add_edge!(net, 1, 2)
        @test add_edge!(net, 2, 2)
        @test ne(net) == 2
        @test has_edge(net, 2, 2)

        # Density: 2 edges of n(n+1)/2 = 6 possible
        @test network_density(net) ≈ 2 / 6

        rem_edge!(net, 2, 2)
        @test ne(net) == 1
    end

    @testset "Density edge cases" begin
        @test isnan(network_density(network(1)))
        @test isnan(network_density(network(0)))
    end

    @testset "rem_vertex! remaps swapped attributes" begin
        # Directed: remove middle vertex; last vertex (4) is renumbered to 2
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 1)
        set_vertex_attribute!(net, :name, Dict(1 => "A", 2 => "B", 3 => "C", 4 => "D"))
        set_edge_attribute!(net, :w, Dict((1, 2) => 12.0, (3, 4) => 34.0, (4, 1) => 41.0))

        @test rem_vertex!(net, 2)
        @test nv(net) == 3

        names = get_vertex_attribute(net, :name)
        @test names[1] == "A"
        @test names[2] == "D"  # old vertex 4 swapped into slot 2
        @test names[3] == "C"
        @test !haskey(names, 4)

        # Edges of old vertex 4 now attach to vertex 2, attributes intact
        @test has_edge(net, 3, 2)
        @test has_edge(net, 2, 1)
        @test get_edge_attribute(net, :w, 3, 2) == 34.0
        @test get_edge_attribute(net, :w, 2, 1) == 41.0
        # Attributes of removed edges are gone
        @test isnothing(get_edge_attribute(net, :w, 1, 2)) || !has_edge(net, 1, 2)

        # Removing the last vertex needs no remap
        net2 = network(3)
        add_edge!(net2, 1, 2)
        set_vertex_attribute!(net2, :name, Dict(1 => "A", 2 => "B", 3 => "C"))
        @test rem_vertex!(net2, 3)
        @test get_vertex_attribute(net2, :name) == Dict(1 => "A", 2 => "B")

        # Undirected: canonical (min,max) edge keys are re-canonicalized
        unet = network(4; directed=false)
        add_edge!(unet, 3, 4)
        set_edge_attribute!(unet, :w, 3, 4, 7.5)
        @test rem_vertex!(unet, 1)  # old vertex 4 becomes 1
        @test has_edge(unet, 1, 3)
        @test get_edge_attribute(unet, :w, 1, 3) == 7.5
    end

    @testset "add_edge! with attribute Dict" begin
        net = network(3)
        @test add_edge!(net, 1, 2, Dict(:weight => 1.5, :type => :friend))
        @test get_edge_attribute(net, :weight, 1, 2) == 1.5
        @test get_edge_attribute(net, :type, 1, 2) == :friend
    end

    @testset "Undirected as_matrix symmetry" begin
        net = network(3; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_edge_attribute!(net, :w, 1, 2, 4.0)
        set_edge_attribute!(net, :w, 2, 3, 5.0)

        A = as_matrix(net)
        @test A == transpose(A)
        @test A[1, 2] == A[2, 1] == 1.0

        W = as_matrix(net; attr=:w)
        @test W == transpose(W)
        @test W[3, 2] == 5.0

        S = as_matrix(net; sparse=true)
        @test S == A
    end

    @testset "Weighted network_from_matrix" begin
        A = [0.0 2.5 0.0; 0.0 0.0 1.5; 0.5 0.0 0.0]
        net = network_from_matrix(A)
        @test ne(net) == 3
        @test get_edge_attribute(net, :weight, 1, 2) == 2.5
        @test get_edge_attribute(net, :weight, 3, 1) == 0.5

        # 0/1-valued but genuinely weighted: force storage
        B = [0 1 0; 0 0 1; 1 0 0]
        net_forced = network_from_matrix(B; store_values=true)
        @test get_edge_attribute(net_forced, :weight, 1, 2) == 1

        # Asymmetric input with directed=false warns
        C = [0 1 0; 0 0 0; 0 0 0]
        @test_logs (:warn, r"Asymmetric") network_from_matrix(C; directed=false)
    end

    @testset "Empty edge list" begin
        net = network_from_edgelist(Tuple{Int,Int}[])
        @test nv(net) == 0
        @test ne(net) == 0
        net2 = network_from_edgelist(Tuple{Int,Int}[]; n=4)
        @test nv(net2) == 4
    end

    @testset "Bipartite semantics" begin
        bnet = BipartiteNetwork(2, 3)
        @test add_edge!(bnet, 1, 3)
        @test add_edge!(bnet, 2, 5)
        @test !add_edge!(bnet, 1, 2)  # within mode 1 rejected
        @test !add_edge!(bnet, 3, 4)  # within mode 2 rejected
        @test ne(bnet) == 2

        # Incidence matrix is n_mode1 × n_mode2
        B = as_matrix(bnet)
        @test size(B) == (2, 3)
        @test B[1, 1] == 1.0  # edge (1, 3): column 3 - 2 = 1
        @test B[2, 3] == 1.0  # edge (2, 5): column 5 - 2 = 3
        @test sum(B) == 2.0

        # A two-mode network is measured over ALL dyads by default (R's
        # discount.bipartite=FALSE); the incidence-matrix denominator is opt-in
        # and honours directedness (BipartiteNetwork defaults to directed, so
        # the old hard-coded n1*n2 = 6 was wrong by a factor of two here).
        @test is_directed(bnet)
        @test network_density(bnet) ≈ 2 / 20
        @test network_density(bnet; discount_bipartite=true) ≈ 2 / 12
        ubnet = BipartiteNetwork(2, 3; directed=false)
        add_edge!(ubnet, 1, 3); add_edge!(ubnet, 2, 5)
        @test network_density(ubnet) ≈ 2 / 10
        @test network_density(ubnet; discount_bipartite=true) ≈ 2 / 6

        # expand_bipartite=true: the (n1+n2)-square adjacency matrix a one-mode
        # matrix routine expects, with the incidence block at A[1:k, k+1:end]
        for directed in (true, false)
            b = BipartiteNetwork(2, 3; directed=directed)
            add_edge!(b, 1, 3)
            add_edge!(b, 2, 5)
            add_edge!(b, 1, 4)
            k = 2
            A = as_matrix(b; expand_bipartite=true)
            @test size(A) == (5, 5)
            @test A[1:k, k+1:end] == as_matrix(b)
            @test A == as_adjacency_matrix(b; expand_bipartite=true)
            if directed
                @test all(A[k+1:end, 1:k] .== 0)      # arcs all run mode 1 -> mode 2 here
            else
                @test A == A'                          # mirrored into the lower block
            end
            @test all(A[1:k, 1:k] .== 0) && all(A[k+1:end, k+1:end] .== 0)  # no within-mode ties
            # ...and it round-trips: the square matrix reproduces the edge set
            back = network_from_matrix(A; directed=directed)
            @test nv(back) == 5 && ne(back) == 3
            @test Set(Tuple(e) for e in edges(back)) == Set(Tuple(e) for e in edges(b.network))
            # Every other keyword still applies to the expanded shape
            @test as_matrix(b; expand_bipartite=true, sparse=true) == A
            set_edge_attribute!(b.network, :w, 1, 3, 2.5)
            @test as_matrix(b; expand_bipartite=true, attr=:w)[1, 3] == 2.5
            set_missing_dyad!(b, 2, 4)
            Am, mask = as_matrix(b; expand_bipartite=true, return_mask=true)
            @test size(mask) == (5, 5) && mask[2, 4] && (directed || mask[4, 2])
            @test count(mask) == (directed ? 1 : 2)
            @test ismissing(as_matrix(b; expand_bipartite=true, missing_value=missing)[2, 4])
            # The flagged Network (no wrapper) behaves identically
            @test as_matrix(b.network; expand_bipartite=true, return_mask=true)[1] == Am
        end
        # A one-mode network ignores the keyword
        one = network(3)
        add_edge!(one, 1, 2)
        @test as_matrix(one; expand_bipartite=true) == as_matrix(one)
    end

    @testset "Graphs.jl conversion" begin
        net = network(3)
        add_edge!(net, 1, 2)
        g = SimpleDiGraph(net)
        @test g isa SimpleDiGraph
        @test Graphs.ne(g) == 1
        # A copy, not a view
        add_edge!(g, 2, 3)
        @test !has_edge(net, 2, 3)

        unet = network(3; directed=false)
        add_edge!(unet, 1, 2)
        h = SimpleGraph(unet)
        @test h isa SimpleGraph
        @test Graphs.ne(h) == 1
    end

    @testset "Pajek round-trip" begin
        mktempdir() do dir
            # Weighted directed network; weights of 1.0 must survive because
            # the file carries an explicit weight column
            net = network(3)
            add_edge!(net, 1, 2)
            add_edge!(net, 2, 3)
            set_vertex_attribute!(net, :vertex_names, Dict(1 => "a", 2 => "b", 3 => "c"))
            set_edge_attribute!(net, :weight, 1, 2, 1.0)
            set_edge_attribute!(net, :weight, 2, 3, 2.5)

            path = joinpath(dir, "net.net")
            write_pajek(net, path)
            back = read_pajek(path)

            @test nv(back) == 3
            @test ne(back) == 2
            @test is_directed(back)
            @test has_edge(back, 1, 2) && has_edge(back, 2, 3)
            @test get_edge_attribute(back, :weight, 1, 2) == 1.0
            @test get_edge_attribute(back, :weight, 2, 3) == 2.5
            @test get_vertex_attribute(back, :vertex_names, 1) == "a"

            # Undirected round-trip
            unet = network(3; directed=false)
            add_edge!(unet, 1, 3)
            upath = joinpath(dir, "unet.net")
            write_pajek(unet, upath)
            uback = read_pajek(upath)
            @test !is_directed(uback)
            @test has_edge(uback, 3, 1)
        end
    end

    @testset "CSV edge list quoting" begin
        mktempdir() do dir
            net = network(2)
            add_edge!(net, 1, 2)
            set_edge_attribute!(net, :label, 1, 2, "hello, \"world\"")

            path = joinpath(dir, "edges.csv")
            write_edgelist_csv(net, path; attrs=[:label])
            lines = readlines(path)
            @test lines[1] == "source,target,label"
            @test lines[2] == "1,2,\"hello, \"\"world\"\"\""
        end
    end

    @testset "DataFrames extension" begin
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_vertex_attribute!(net, :name, Dict(1 => "A", 2 => "B", 3 => "C"))
        set_edge_attribute!(net, :w, 1, 2, 9.0)

        edf = as_dataframe(net; attrs=[:w])
        @test nrow(edf) == 2
        @test Set(names(edf)) == Set(["source", "target", "w"])

        vdf = as_dataframe(net; vertices=true)
        @test nrow(vdf) == 3
        @test vdf.name == ["A", "B", "C"]

        df = DataFrame(source=[1, 2], target=[2, 3], weight=[1.0, 2.0])
        net2 = network_from_dataframe(df)
        @test nv(net2) == 3
        @test ne(net2) == 2
        @test get_edge_attribute(net2, :weight, 2, 3) == 2.0
    end

    @testset "Attribute indexing shorthand" begin
        net = network(2)
        add_edge!(net, 1, 2)
        net[:v, :name] = Dict(1 => "A", 2 => "B")
        net[:e, :w] = Dict((1, 2) => 3.0)
        net[:n, :title] = "T"
        @test net[:v, :name][1] == "A"
        @test net[:e, :w][(1, 2)] == 3.0
        @test net[:n, :title] == "T"
    end

    @testset "Vertex Permutation" begin
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_vertex_attribute!(net, :name, Dict(1 => "A", 2 => "B", 3 => "C"))

        # Permute: old 1 -> new 3, old 2 -> new 1, old 3 -> new 2
        perm = [2, 3, 1]  # new position i has old vertex perm[i]
        net_perm = permute_vertices(net, perm)

        @test nv(net_perm) == 3
        @test ne(net_perm) == 2

        # Check vertex attributes were permuted correctly
        names = get_vertex_attribute(net_perm, :name)
        @test names[1] == "B"  # old vertex 2 is now vertex 1
        @test names[2] == "C"  # old vertex 3 is now vertex 2
        @test names[3] == "A"  # old vertex 1 is now vertex 3
    end

    @testset "Copy independence" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_vertex_attribute!(net, :name, Dict(1 => "A", 2 => "B"))
        set_edge_attribute!(net, :weight, Dict((1, 2) => 1.5))
        set_network_attribute!(net, :title, "original")

        net2 = copy(net)
        @test nv(net2) == 4
        @test ne(net2) == 2
        @test has_edge(net2, 1, 2)
        @test has_edge(net2, 2, 3)

        # Mutating the copy's structure does not affect the original
        add_edge!(net2, 3, 4)
        rem_edge!(net2, 1, 2)
        add_vertex!(net2)
        @test !has_edge(net, 3, 4)
        @test has_edge(net, 1, 2)
        @test nv(net) == 4
        @test ne(net) == 2

        # Mutating the copy's attributes does not affect the original
        set_vertex_attribute!(net2, :name, 2, "Z")
        set_vertex_attribute!(net2, :age, Dict(1 => 30))
        set_edge_attribute!(net2, :weight, 2, 3, 9.9)
        set_network_attribute!(net2, :title, "copy")
        delete_vertex_attribute!(net2, :name)
        @test get_vertex_attribute(net, :name, 2) == "B"
        @test :name in list_vertex_attributes(net)
        @test !(:age in list_vertex_attributes(net))
        @test get_edge_attribute(net, :weight, 1, 2) == 1.5
        @test isnothing(get_edge_attribute(net, :weight, 2, 3))
        @test get_network_attribute(net, :title) == "original"

        # And the reverse: mutating the original does not affect the copy
        net3 = copy(net)
        set_vertex_attribute!(net, :name, 1, "X")
        @test get_vertex_attribute(net3, :name, 1) == "A"
    end

    @testset "Copy preserves invariants" begin
        # Directed flag
        dnet = network(3)
        @test is_directed(copy(dnet))
        unet = network(3; directed=false)
        add_edge!(unet, 1, 2)
        unet2 = copy(unet)
        @test !is_directed(unet2)
        @test has_edge(unet2, 2, 1)  # undirected symmetry survives the copy
        @test ne(unet2) == 1

        # Loops flag
        lnet = network(3; loops=true)
        add_edge!(lnet, 1, 1)
        lnet2 = copy(lnet)
        @test lnet2.loops
        @test has_edge(lnet2, 1, 1)
        @test add_edge!(lnet2, 2, 2)
        @test !copy(network(3)).loops

        # Bipartite field on Network
        bip = network(5; bipartite=2)
        bip2 = copy(bip)
        @test bip2.bipartite == 2
        @test !add_edge!(bip2, 1, 2)  # within-mode edge still rejected
        @test add_edge!(bip2, 1, 3)

        # All three attribute namespaces are preserved
        net = network(3)
        add_edge!(net, 1, 2)
        set_vertex_attribute!(net, :name, Dict(1 => "A"))
        set_edge_attribute!(net, :weight, Dict((1, 2) => 2.0))
        set_network_attribute!(net, :title, "t")
        net2 = copy(net)
        @test list_vertex_attributes(net2) == [:name]
        @test list_edge_attributes(net2) == [:weight]
        @test list_network_attributes(net2) == [:title]
        @test get_vertex_attribute(net2, :name, 1) == "A"
        @test get_edge_attribute(net2, :weight, 1, 2) == 2.0
        @test get_network_attribute(net2, :title) == "t"

        # BipartiteNetwork wrapper copy
        bnet = BipartiteNetwork(2, 3)
        add_edge!(bnet, 1, 3)
        set_network_attribute!(bnet.network, :title, "modes")
        bnet2 = copy(bnet)
        @test bnet2 isa BipartiteNetwork{Int}
        @test bnet2.n_mode1 == 2
        @test bnet2.n_mode2 == 3
        @test has_edge(bnet2, 1, 3)
        @test get_network_attribute(bnet2.network, :title) == "modes"
        add_edge!(bnet2, 2, 4)
        @test !has_edge(bnet, 2, 4)

        # Missing-dyad mask survives copy, and copies are independent
        mnet = network(4)
        add_edge!(mnet, 1, 2)
        set_missing_dyad!(mnet, 3, 4)
        set_missing_dyad!(mnet, 1, 2)
        mnet2 = copy(mnet)
        @test n_missing_dyads(mnet2) == 2
        @test is_missing_dyad(mnet2, 3, 4)
        @test is_missing_dyad(mnet2, 1, 2)
        # face values also survive
        @test has_edge(mnet2, 1, 2)
        @test !has_edge(mnet2, 3, 4)
        # independence in both directions
        delete_missing_dyad!(mnet2, 3, 4)
        @test is_missing_dyad(mnet, 3, 4)
        set_missing_dyad!(mnet, 2, 3)
        @test !is_missing_dyad(mnet2, 2, 3)
    end

    @testset "Missing-dyad mask" begin
        net = network(5)
        @test n_missing_dyads(net) == 0
        @test !is_missing_dyad(net, 1, 2)
        @test isempty(collect(missing_dyads(net)))

        # Marking is metadata only: no edge is created or removed
        add_edge!(net, 1, 2)
        @test set_missing_dyad!(net, 1, 2) === net
        set_missing_dyad!(net, 4, 5)
        @test has_edge(net, 1, 2)       # face value unchanged
        @test !has_edge(net, 4, 5)      # face value unchanged
        @test ne(net) == 1
        @test is_missing_dyad(net, 1, 2)
        @test is_missing_dyad(net, 4, 5)
        @test !is_missing_dyad(net, 2, 1)   # directed: arcs masked independently
        @test n_missing_dyads(net) == 2
        @test sort(collect(missing_dyads(net))) == [(1, 2), (4, 5)]

        # Idempotent marking
        set_missing_dyad!(net, 1, 2)
        @test n_missing_dyads(net) == 2

        # Delete and clear
        @test delete_missing_dyad!(net, 1, 2)
        @test !delete_missing_dyad!(net, 1, 2)
        @test !is_missing_dyad(net, 1, 2)
        set_missing_dyad!(net, 2, 3)
        @test clear_missing_dyads!(net) === net
        @test n_missing_dyads(net) == 0

        # Validation
        @test_throws ArgumentError set_missing_dyad!(net, 1, 6)
        @test_throws ArgumentError set_missing_dyad!(net, 0, 2)
        @test_throws ArgumentError set_missing_dyad!(net, 2, 2)  # no loops allowed
        lnet = network(3; loops=true)
        set_missing_dyad!(lnet, 2, 2)
        @test is_missing_dyad(lnet, 2, 2)

        # Undirected: symmetric handling
        unet = network(4; directed=false)
        set_missing_dyad!(unet, 3, 1)
        @test is_missing_dyad(unet, 1, 3)
        @test is_missing_dyad(unet, 3, 1)
        @test n_missing_dyads(unet) == 1
        @test collect(missing_dyads(unet)) == [(1, 3)]   # canonical minmax form
        set_missing_dyad!(unet, 1, 3)                    # same dyad, both orders
        @test n_missing_dyads(unet) == 1
        @test delete_missing_dyad!(unet, 1, 3)
        @test n_missing_dyads(unet) == 0

        # Two-mode metadata: within-mode dyads cannot be masked
        bip = network(5; bipartite=2)
        set_missing_dyad!(bip, 1, 4)
        @test is_missing_dyad(bip, 1, 4)
        @test_throws ArgumentError set_missing_dyad!(bip, 1, 2)

        # BipartiteNetwork forwarding
        bnet = BipartiteNetwork(2, 2)
        @test set_missing_dyad!(bnet, 1, 3) === bnet
        @test is_missing_dyad(bnet, 1, 3)
        @test n_missing_dyads(bnet) == 1
        @test collect(missing_dyads(bnet)) == [(1, 3)]
        @test delete_missing_dyad!(bnet, 1, 3)
        set_missing_dyad!(bnet, 2, 4)
        @test clear_missing_dyads!(bnet) === bnet
        @test n_missing_dyads(bnet) == 0

        # rem_vertex! drops marks of the removed vertex and remaps the
        # swapped-in last vertex
        rnet = network(4)
        set_missing_dyad!(rnet, 1, 2)   # untouched
        set_missing_dyad!(rnet, 2, 3)   # involves removed vertex 3 -> dropped
        set_missing_dyad!(rnet, 1, 4)   # vertex 4 is renumbered to 3
        rem_vertex!(rnet, 3)
        @test nv(rnet) == 3
        @test is_missing_dyad(rnet, 1, 2)
        @test is_missing_dyad(rnet, 1, 3)   # old (1, 4)
        @test n_missing_dyads(rnet) == 2

        # Induced subgraph keeps marks between retained vertices (relabeled)
        snet = network(4)
        set_missing_dyad!(snet, 1, 3)
        set_missing_dyad!(snet, 2, 4)
        sub = get_induced_subgraph(snet, [1, 3])
        @test n_missing_dyads(sub) == 1
        @test is_missing_dyad(sub, 1, 2)    # old (1, 3)

        # Vertex permutation relabels marks
        pnet = network(3)
        set_missing_dyad!(pnet, 1, 2)
        perm = permute_vertices(pnet, [3, 1, 2])   # old 3 -> new 1, old 1 -> new 2, old 2 -> new 3
        @test is_missing_dyad(perm, 2, 3)
        @test n_missing_dyads(perm) == 1
    end

    @testset "zero" begin
        # Type-level zero: empty directed network of the same eltype
        z = zero(Network{Int})
        @test z isa Network{Int}
        @test nv(z) == 0
        @test ne(z) == 0
        @test is_directed(z)
        @test zero(Network) isa Network{Int}

        # Instance-level zero matches directedness (and loops)
        unet = network(5; directed=false, loops=true)
        add_edge!(unet, 1, 2)
        uz = zero(unet)
        @test nv(uz) == 0
        @test ne(uz) == 0
        @test !is_directed(uz)
        @test uz.loops
        @test is_directed(zero(network(2)))
    end

    @testset "add_vertices! extends Graphs" begin
        @test NetworkCore.add_vertices! === Graphs.add_vertices!

        net = network(2)
        @test add_vertices!(net, 3) == 3
        @test nv(net) == 5
        @test add_vertices!(net, 0) == 0
        @test nv(net) == 5

        # The generic Graphs function dispatches to the Network method
        @test Graphs.add_vertices!(net, 2) == 2
        @test nv(net) == 7
    end

    @testset "Graphs.jl interface conformance" begin
        @testset "directed" begin
            net = network(5)  # Network{Int, true}
            add_edge!(net, 1, 2)
            add_edge!(net, 2, 3)

            # Directedness is truthful at both the instance and type level
            @test is_directed(net)
            @test Graphs.is_directed(typeof(net))
            @test net isa Network{Int}          # one-param UnionAll still matches
            @test net isa Network{Int, true}

            @test nv(net) == 5
            @test ne(net) == 2
            @test eltype(net) == Int
            @test edgetype(net) == Graphs.SimpleEdge{Int}

            # Edge iterator: correct count, length, and eltype
            es = edges(net)
            @test length(es) == 2
            @test eltype(es) == Graphs.SimpleEdge{Int}
            @test length(collect(es)) == 2
            @test Graphs.SimpleEdge(1, 2) in collect(es)

            # Directed has_edge is asymmetric
            @test has_edge(net, 1, 2)
            @test !has_edge(net, 2, 1)

            # Degree = indegree + outdegree for directed graphs
            @test Graphs.degree(net, 2) == 2
            @test Graphs.indegree(net, 2) == 1
            @test Graphs.outdegree(net, 2) == 1
            @test Graphs.degree(net) ==
                  Graphs.indegree(net) .+ Graphs.outdegree(net)

            # zero preserves type (including directedness)
            z = zero(typeof(net))
            @test z isa Network{Int, true}
            @test nv(z) == 0 && ne(z) == 0

            # Graphs generics see directed semantics
            @test Graphs.density(net) == 2 / 20
        end

        @testset "undirected" begin
            # The acceptance example: 5 nodes, 2 undirected edges
            net = network(5; directed=false)
            add_edge!(net, 1, 2)
            add_edge!(net, 2, 3)

            # Directedness is truthful at both the instance and type level
            @test !is_directed(net)
            @test !Graphs.is_directed(typeof(net))
            @test net isa Network{Int}          # one-param UnionAll still matches
            @test net isa Network{Int, false}

            @test nv(net) == 5
            @test ne(net) == 2  # each undirected edge counted once

            # Edge iterator yields each undirected edge exactly once,
            # with correct length and eltype
            es = edges(net)
            @test length(es) == ne(net) == 2
            @test eltype(es) == Graphs.SimpleEdge{Int}
            collected = collect(es)
            @test length(collected) == 2
            @test all(e -> src(e) <= dst(e), collected)
            @test Set((src(e), dst(e)) for e in collected) ==
                  Set([(1, 2), (2, 3)])

            # has_edge is symmetric
            @test has_edge(net, 1, 2) && has_edge(net, 2, 1)
            @test has_edge(net, 2, 3) && has_edge(net, 3, 2)
            @test !has_edge(net, 1, 3)

            # Undirected degree counts each incident edge once, and
            # in/out/neighbors all coincide
            @test Graphs.degree(net) == [1, 2, 1, 0, 0]
            @test Graphs.degree(net, 2) == 2
            @test inneighbors(net, 2) == outneighbors(net, 2) == neighbors(net, 2)
            @test sort(neighbors(net, 2)) == [1, 3]
            # Handshake lemma: degrees sum to twice the edge count
            @test sum(Graphs.degree(net)) == 2 * ne(net)

            # zero preserves type (including directedness)
            z = zero(typeof(net))
            @test z isa Network{Int, false}
            @test nv(z) == 0 && ne(z) == 0

            # ACCEPTANCE TEST: Graphs generics now give the correct
            # undirected answers (previously density returned 0.1)
            @test Graphs.density(net) == 0.2
            adj = Matrix(Graphs.adjacency_matrix(net))
            @test adj == adj'  # symmetric
            @test sum(adj) == 2 * ne(net)
            @test Graphs.connected_components(net) ==
                  [[1, 2, 3], [4], [5]]

            # Undirected self-loops are counted once by ne
            lnet = network(3; directed=false, loops=true)
            add_edge!(lnet, 1, 2)
            add_edge!(lnet, 3, 3)
            @test ne(lnet) == 2
            @test length(collect(edges(lnet))) == 2
        end

        @testset "constructor forms with directedness" begin
            # One-parameter form dispatches on the directed kwarg
            u = Network{Int}(; n=4, directed=false)
            @test u isa Network{Int, false}
            @test !is_directed(u)
            d = Network{Int}(; n=4)
            @test d isa Network{Int, true}
            @test is_directed(d)

            # network() picks D from the keyword too
            @test network(3; directed=false) isa Network{Int, false}
            @test network(3) isa Network{Int, true}

            # directedness is a type parameter, not a mutable field:
            # net.directed reads it, assigning it throws
            @test u.directed === false
            @test d.directed === true
            @test_throws ArgumentError (u.directed = true)

            # D must be a Bool
            @test_throws ArgumentError Network{Int, 3}()
        end

        @testset "honest is_bipartite vs is_two_mode" begin
            # Graph-theoretic bipartiteness is computed, not read off metadata
            path = network(4; directed=false)
            add_edge!(path, 1, 2)
            add_edge!(path, 2, 3)
            @test Graphs.is_bipartite(path)   # a path is two-colorable
            @test !is_two_mode(path)          # ... but carries no mode metadata

            add_edge!(path, 1, 3)             # odd cycle
            @test !Graphs.is_bipartite(path)

            # The statnet-style metadata flag lives under is_two_mode
            two_mode = network(5; bipartite=2)
            @test is_two_mode(two_mode)
            @test Graphs.is_bipartite(two_mode)  # no edges: trivially bipartite

            bnet = BipartiteNetwork(2, 3)
            @test is_two_mode(bnet)
            @test Graphs.is_bipartite(bnet)
        end

        @testset "BipartiteNetwork neighbor forwarding" begin
            bnet = BipartiteNetwork(2, 3)
            add_edge!(bnet, 1, 3)
            add_edge!(bnet, 2, 3)
            @test sort(outneighbors(bnet, 1)) == [3]
            @test sort(inneighbors(bnet, 3)) == [1, 2]
            @test sort(neighbors(bnet, 1)) == [3]
            # degree forwards through the wrapped network (directed default)
            @test Graphs.degree(bnet.network, 3) == 2
        end
    end

    @testset "Typed attribute accessors" begin
        net = network(4; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_vertex_attribute!(net, :wealth, Dict(1 => 10, 2 => 36, 3 => 55, 4 => 44))
        set_edge_attribute!(net, :weight, Dict((1, 2) => 1.5, (2, 3) => 2.0))

        # Typed vertex accessor: concrete return type, single conversion pass
        w = @inferred Dict{Int, Int} get_vertex_attribute(net, :wealth, Int)
        @test w isa Dict{Int, Int}
        @test w == Dict(1 => 10, 2 => 36, 3 => 55, 4 => 44)
        wf = @inferred Dict{Int, Float64} get_vertex_attribute(net, :wealth, Float64)
        @test wf isa Dict{Int, Float64}
        @test wf[3] === 55.0
        # Missing attribute -> empty typed dict, not an error
        @test @inferred(Dict{Int, Int}, get_vertex_attribute(net, :nope, Int)) ==
              Dict{Int, Int}()
        # Unconvertible values throw rather than silently corrupting
        set_vertex_attribute!(net, :label, Dict(1 => "a"))
        @test_throws MethodError get_vertex_attribute(net, :label, Int)

        # Vector form: dense, vertex-ordered, concretely typed
        wv = @inferred vertex_attribute_vector(net, :wealth, Int)
        @test wv isa Vector{Int}
        @test wv == [10, 36, 55, 44]
        wvf = @inferred vertex_attribute_vector(net, :wealth, Float64)
        @test wvf == [10.0, 36.0, 55.0, 44.0]
        # Missing entries: default fills, no default throws
        set_vertex_attribute!(net, :partial, Dict(1 => 7, 3 => 9))
        @test vertex_attribute_vector(net, :partial, Int; default=0) == [7, 0, 9, 0]
        # No default: an ArgumentError naming attribute and vertex (was a bare
        # `KeyError(i)`, which did not even say which attribute was missing)
        err = try vertex_attribute_vector(net, :partial, Int); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin(":partial is not set for vertex 2", err.msg)
        @test occursin("default=", err.msg)
        @test !occursin("existing:", err.msg)          # the attribute exists; only a gap
        err = try vertex_attribute_vector(net, :nope, Int); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("no attribute :nope", err.msg)
        @test occursin(":wealth", err.msg)             # lists what the network does have
        @test vertex_attribute_vector(net, :nope, Int; default=-1) == fill(-1, 4)

        # Typed edge accessor with canonical (minmax) keys on undirected nets
        ew = @inferred Dict{Tuple{Int, Int}, Float64} get_edge_attribute(net, :weight, Float64)
        @test ew isa Dict{Tuple{Int, Int}, Float64}
        @test ew == Dict((1, 2) => 1.5, (2, 3) => 2.0)
        @test @inferred(Dict{Tuple{Int, Int}, Int},
                        get_edge_attribute(net, :nope, Int)) ==
              Dict{Tuple{Int, Int}, Int}()

        # Untyped accessors are unchanged
        @test get_vertex_attribute(net, :wealth) isa Dict{Int, Any}
        # The whole-Dict getters allocate nothing on a hit: the eager default
        # `get(d, k, Dict())` built an 80 B empty Dict on EVERY call, paid per
        # sampler step by a term reading the getter (an ERGMUserterms term).
        get_vertex_attribute(net, :wealth); get_edge_attribute(net, :wealth)
        @test @allocated(get_vertex_attribute(net, :wealth)) == 0
        set_edge_attribute!(net, :w, first(edges(net)).src, first(edges(net)).dst, 1.0)
        get_edge_attribute(net, :w)
        @test @allocated(get_edge_attribute(net, :w)) == 0
        # A miss still returns an empty, correctly typed Dict
        @test get_vertex_attribute(net, :absent) == Dict{Int, Any}()
        @test get_edge_attribute(net, :absent) == Dict{Tuple{Int,Int}, Any}()
        @test get_vertex_attribute(net, :wealth, 2) == 36
        @test get_edge_attribute(net, :weight) isa Dict{Tuple{Int, Int}, Any}
    end

    @testset "Datasets" begin
        @testset "florentine_marriage" begin
            flo = load_dataset(:florentine_marriage)
            @test nv(flo) == 16
            @test ne(flo) == 20
            @test !is_directed(flo)
            # Raw statistics for the standard statnet flomarriage data:
            # 20 undirected ties among 16 families, density 20/120.
            # (ERGM.jl's Florentine golden masters can adopt this loader.)
            @test network_density(flo) ≈ 20 / 120
            @test has_edge(flo, 1, 9)     # Acciaiuoli - Medici
            @test has_edge(flo, 9, 1)     # undirected: both directions
            @test !has_edge(flo, 1, 2)
            # Pucci (12) is the isolate
            @test isempty(neighbors(flo, 12))
            # Pinned vertex attributes (statnet flomarriage)
            @test get_vertex_attribute(flo, :name, 9) == "Medici"
            @test get_vertex_attribute(flo, :name, 15) == "Strozzi"
            @test get_vertex_attribute(flo, :wealth, 9) == 103
            @test get_vertex_attribute(flo, :wealth, 15) == 146
            @test get_vertex_attribute(flo, :priorates, 9) == 53
            @test get_vertex_attribute(flo, :totalties, 9) == 54
            @test get_vertex_attribute(flo, :priorates, 12) == 0
            @test vertex_attribute_vector(flo, :wealth, Int) ==
                  [10, 36, 55, 44, 20, 32, 8, 42, 103, 48, 49, 3, 27, 10, 146, 48]
        end

        @testset "florentine_business" begin
            flob = load_dataset(:florentine_business)
            @test nv(flob) == 16
            @test ne(flob) == 15
            @test !is_directed(flob)
            @test has_edge(flob, 3, 5)    # Barbadori - Castellani
            @test has_edge(flob, 9, 16)   # Medici - Tornabuoni
            @test !has_edge(flob, 1, 9)   # marriage tie, not a business tie
            # Same family attributes as the marriage network
            @test get_vertex_attribute(flob, :name, 3) == "Barbadori"
            @test get_vertex_attribute(flob, :wealth, 3) == 55
            @test get_vertex_attribute(flob, :totalties, 5) == 18
        end

        @testset "sampson" begin
            samp = load_dataset(:sampson)
            @test nv(samp) == 18
            @test ne(samp) == 88
            @test is_directed(samp)
            @test has_edge(samp, 1, 2)
            @test !has_edge(samp, 2, 3)   # liking is not symmetric
            # Pinned vertex attributes (statnet samplike)
            @test get_vertex_attribute(samp, :name, 1) == "John Bosco"
            @test get_vertex_attribute(samp, :name, 18) == "Simplicius"
            @test get_vertex_attribute(samp, :group, 1) == "Turks"
            @test get_vertex_attribute(samp, :group, 4) == "Loyal"
            @test get_vertex_attribute(samp, :group, 18) == "Outcasts"
            @test get_vertex_attribute(samp, :cloisterville, 1) === true
            @test get_vertex_attribute(samp, :cloisterville, 18) === false
            groups = vertex_attribute_vector(samp, :group, String)
            @test count(==("Turks"), groups) == 7
            @test count(==("Loyal"), groups) == 7
            @test count(==("Outcasts"), groups) == 4
            # Edge attribute: number of time points the tie was named (1-3)
            noms = get_edge_attribute(samp, :nominations, Int)
            @test length(noms) == 88
            @test noms[(2, 1)] == 3
            @test noms[(1, 2)] == 1
            @test all(1 <= v <= 3 for v in values(noms))
        end

        @testset "faux_mesa_high" begin
            fmh = load_dataset(:faux_mesa_high)
            @test fmh isa Network{Int,false}
            @test nv(fmh) == 205
            @test ne(fmh) == 203
            @test !is_directed(fmh)
            @test Set(list_vertex_attributes(fmh)) == Set([:Grade, :Race, :Sex])
            grade = vertex_attribute_vector(fmh, :Grade, Int)
            @test all(in(7:12), grade)
            # statnet's own grade distribution (table(faux.mesa.high %v% "Grade"))
            @test [count(==(g), grade) for g in 7:12] == [62, 40, 42, 25, 24, 12]
            race = vertex_attribute_vector(fmh, :Race, String)
            @test Set(race) == Set(["Black", "Hisp", "NatAm", "Other", "White"])
            @test count(==("Hisp"), race) == 109
            sex = vertex_attribute_vector(fmh, :Sex, String)
            @test Set(sex) == Set(["F", "M"])
            @test count(==("F"), sex) == 99
            # Same-grade ties: the nodematch("Grade") statistic. Cross-checked
            # against the provenanced ERGMEgo.jl golden fixture
            # test/fixtures/fauxmesa_ego_census.toml (ergm.ego 1.1.4 / ergm
            # 4.12.0, a census of faux.mesa.high), whose target statistics are
            # `targets = [203, 163]` for edges + nodematch("Grade").
            same_grade = count(grade[src(e)] == grade[dst(e)] for e in edges(fmh))
            @test same_grade == 163
            # No self-loops, every edge stored once
            @test all(src(e) != dst(e) for e in edges(fmh))
            @test length(collect(edges(fmh))) == 203
        end

        @testset "s50" begin
            s50 = load_dataset(:s50)
            @test s50 isa NamedTuple
            @test propertynames(s50) == (:friendship, :alcohol, :smoke)
            @test s50.friendship isa Vector{Network{Int,true}}
            @test length(s50.friendship) == 3
            @test all(net -> nv(net) == 50 && is_directed(net), s50.friendship)
            # Tie counts equal the row sums of the source matrices. Computed
            # from the files here and frozen as literals below: the RSiena
            # s501/s502/s503 matrices carry 113, 116 and 122 ties (verified
            # against Siena.jl's test/data copies and data/README.md).
            for (w, net) in enumerate(s50.friendship)
                A = NetworkCore._read_int_matrix_tsv(joinpath("s50", "friendship_wave$w.tsv"))
                @test size(A) == (50, 50)
                @test ne(net) == sum(A)
                @test as_matrix(net) == A
                @test all(A[i, i] == 0 for i in 1:50)
            end
            @test ne.(s50.friendship) == [113, 116, 122]
            # The binary waves carry no edge attribute (store_values=false)
            @test all(isempty(list_edge_attributes(net)) for net in s50.friendship)
            @test s50.alcohol isa Matrix{Int}
            @test size(s50.alcohol) == (50, 3)
            @test all(in(1:5), s50.alcohol)
            @test s50.smoke isa Matrix{Int}
            @test size(s50.smoke) == (50, 3)
            @test all(in(1:3), s50.smoke)
            # First rows of RSiena::s50a / s50s
            @test s50.alcohol[1, :] == [3, 1, 3]
            @test s50.smoke[1, :] == [2, 1, 1]
        end

        @testset "wtc_police_calls" begin
            wtc = load_dataset(:wtc_police_calls)
            @test wtc isa NamedTuple
            @test propertynames(wtc) == (:events, :n_actors, :is_icr)
            @test wtc.events isa Matrix{Int}
            @test size(wtc.events) == (481, 3)
            @test wtc.events[:, 1] == 1:481            # ordinal event numbers
            @test all(in(1:37), wtc.events[:, 2])
            @test all(in(1:37), wtc.events[:, 3])
            @test all(wtc.events[:, 2] .!= wtc.events[:, 3])   # no self-calls
            @test wtc.events[1, :] == [1, 16, 32]      # head(WTCPoliceCalls)
            @test wtc.events[481, :] == [481, 16, 32]  # tail(WTCPoliceCalls)
            @test wtc.n_actors == 37
            @test wtc.is_icr isa Vector{Bool}
            @test length(wtc.is_icr) == 37
            # Frozen from the R object WTCPoliceIsICR: three coordinators
            @test count(wtc.is_icr) == 3
            @test findall(wtc.is_icr) == [13, 14, 32]
            # 35 of the 37 actors ever participate; the universe is still 37
            @test length(union(wtc.events[:, 2], wtc.events[:, 3])) == 35
        end

        @testset "independence and errors" begin
            # Each call returns a fresh, independent network
            a = load_dataset(:florentine_marriage)
            b = load_dataset(:florentine_marriage)
            add_edge!(a, 1, 2)
            @test !has_edge(b, 1, 2)
            set_vertex_attribute!(a, :wealth, 9, 0)
            @test get_vertex_attribute(b, :wealth, 9) == 103

            # ...and so do the NamedTuple datasets (no cached, shared arrays)
            s1 = load_dataset(:s50)
            s2 = load_dataset(:s50)
            @test s1.friendship[1] !== s2.friendship[1]
            @test s1.alcohol !== s2.alcohol
            add_edge!(s1.friendship[1], 1, 2)
            s1.alcohol[1, 1] = 99
            @test !has_edge(s2.friendship[1], 1, 2)
            @test s2.alcohol[1, 1] == 3
            w1 = load_dataset(:wtc_police_calls)
            w2 = load_dataset(:wtc_police_calls)
            @test w1.events !== w2.events
            w1.events[1, 2] = 0
            @test w2.events[1, 2] == 16

            @test_throws ArgumentError load_dataset(:no_such_dataset)
            # The error lists every dataset by name
            msg = sprint(showerror, try
                load_dataset(:no_such_dataset)
            catch e
                e
            end)
            @test occursin(":no_such_dataset", msg)
            for name in (:florentine_marriage, :florentine_business, :sampson,
                         :faux_mesa_high, :s50, :wtc_police_calls)
                @test occursin(":$name", msg)
            end
        end
    end

    @testset "Result presentation" begin
        @testset "format_pvalue" begin
            # Fixed-point display
            @test format_pvalue(0.03511) == "0.0351"
            @test format_pvalue(1.0) == "1.0000"
            @test format_pvalue(0.5; digits=2) == "0.50"
            # Scientific notation below the fixed-point resolution
            @test format_pvalue(3.2e-7) == "3.2e-07"
            @test format_pvalue(0.00005) == "5.0e-05"
            @test format_pvalue(0.00005; digits=5) == "0.00005"
            # Display floor: never prints 0.0
            @test format_pvalue(0.0) == "<1e-16"
            @test format_pvalue(1e-300) == "<1e-16"
            @test format_pvalue(5e-17) == "<1e-16"
            @test format_pvalue(1e-16) == "1.0e-16"   # at the floor, not below
            @test format_pvalue(0.0; floor=1e-10) == "<1e-10"
            # NaN passes through (e.g. undefined standard error)
            @test format_pvalue(NaN) == "NaN"
            # Invalid input
            @test_throws ArgumentError format_pvalue(-0.1)
            @test_throws ArgumentError format_pvalue(1.5)
            @test_throws ArgumentError format_pvalue(0.5; digits=0)
        end

        @testset "z_pvalues" begin
            # The ONE z → two-sided-normal p helper.
            # Reference: Distributions' ccdf, which NetworkCore.jl deliberately does
            # NOT depend on — erfc from SpecialFunctions is the implementation.
            grid = vcat(collect(-37.0:0.25:37.0), [-8.3, 8.3, 1.96, -1.96, 2.5758, 1e-3, -1e-3])
            p = z_pvalues(grid)
            @test p isa Vector{Float64}
            @test length(p) == length(grid)
            @test all(isfinite, p) && all(0 .< p .<= 1)
            @test p ≈ 2 .* ccdf.(Normal(), abs.(grid)) rtol = 1e-12
            @test z_pvalues([0.0]) == [1.0]
            # Where the naive 2(1 − Φ(|z|)) cancels to nonsense, erfc does not
            @test z_pvalues([8.3])[1] ≈ 2 * ccdf(Normal(), 8.3) rtol = 1e-12
            @test !(z_pvalues([8.3])[1] ≈ 2 * (1 - cdf(Normal(), 8.3)))
            @test z_pvalues([37.0])[1] ≈ 2 * ccdf(Normal(), 37.0) rtol = 1e-12
            @test z_pvalues([37.0])[1] > floatmin(Float64)
            # Floor: a finite statistic never has p exactly 0
            @test z_pvalues([50.0]) == [floatmin(Float64)]
            @test z_pvalues([-50.0]) == [floatmin(Float64)]
            @test z_pvalues([40.0])[1] >= floatmin(Float64)
            @test z_pvalues([Inf])[1] == 0.0          # an infinite statistic may
            @test z_pvalues([-Inf])[1] == 0.0
            # NaN stays NaN
            @test isnan(z_pvalues([NaN])[1])
            @test isequal(z_pvalues([NaN, 0.0, NaN]), [NaN, 1.0, NaN])
            # Symmetric and monotone decreasing in |z|
            @test z_pvalues([-2.5]) == z_pvalues([2.5])
            zs = collect(0.0:0.1:45.0)
            pz = z_pvalues(zs)
            @test issorted(pz; rev=true)                   # incl. through the underflow
            @test all(pz[k] > pz[k + 1] for k in 1:370)   # strictly, up to |z| = 37
            @test all(pz[391:end] .== floatmin(Float64))  # floored from |z| = 39, never subnormal
            @test all(p -> p >= floatmin(Float64), pz)
            # Any Real input, Float64 output
            @test z_pvalues([0, 1, 2]) ≈ z_pvalues([0.0, 1.0, 2.0])
            @test z_pvalues(Float32[1.5])[1] ≈ z_pvalues([1.5])[1] rtol = 1e-6
            @test z_pvalues(Float64[]) == Float64[]

            # Two-argument form: z = est/se, NaN wherever se is not positive
            zp = z_pvalues([1.2, -0.3, 0.8, 2.0, 1.0], [0.4, 0.1, 0.0, -1.0, NaN])
            @test zp isa NamedTuple{(:z, :p)}
            @test zp.z[1:2] == [1.2 / 0.4, -0.3 / 0.1]
            @test zp.z[1:2] ≈ [3.0, -3.0]
            @test all(isnan, zp.z[3:5]) && all(isnan, zp.p[3:5])
            @test zp.p[1:2] == z_pvalues(zp.z[1:2])
            @test zp.p[1:2] ≈ z_pvalues([3.0, -3.0])
            @test isnan(z_pvalues([1.0], [0.0]).p[1])
            @test isnan(z_pvalues([1.0], [0.0]).z[1])
            @test_throws ArgumentError z_pvalues([1.0, 2.0], [0.5])
            # Feeds print_coeftable directly
            io = IOBuffer()
            print_coeftable(io, ["a", "b"], [1.2, 0.8], [0.4, 0.0], zp.p[[1, 3]]; z_values=zp.z[[1, 3]])
            out = String(take!(io))
            @test occursin("0.0027", out) && occursin("NaN", out)
        end

        @testset "signif_code" begin
            @test signif_code(0.0005) == "***"
            @test signif_code(0.005) == "**"
            @test signif_code(0.03) == "*"
            @test signif_code(0.07) == "."
            @test signif_code(0.2) == ""
            @test signif_code(NaN) == ""
            # Strict boundaries (matching the ecosystem's existing show methods)
            @test signif_code(0.001) == "**"
            @test signif_code(0.01) == "*"
            @test signif_code(0.05) == "."
            @test signif_code(0.1) == ""
        end

        @testset "print_coeftable" begin
            names = ["edges", "gwesp.fixed.0.5"]
            est = [-3.2145, 0.8113]
            se = [0.1936, 0.2854]
            p = [1.2e-62, 0.0045]
            out = sprint(io -> print_coeftable(io, names, est, se, p))
            lines = split(out, '\n')

            # Header row with the four R-style column titles
            header = lines[1]
            for h in ("Estimate", "Std.Error", "z value", "Pr(>|z|)")
                @test occursin(h, header)
            end

            # Column alignment: each numeric value is right-aligned so it
            # ends where its header ends, on every row
            for (title, values) in (("Estimate", ["-3.2145", "0.8113"]),
                                    ("Std.Error", ["0.1936", "0.2854"]),
                                    ("z value", ["-16.6038", "2.8427"]),
                                    ("Pr(>|z|)", ["<1e-16", "0.0045"]))
                stop = findfirst(title, header).stop
                for (i, val) in enumerate(values)
                    row = lines[1 + i]
                    r = findfirst(val, row)
                    @test r !== nothing && r.stop == stop
                end
            end

            # Name column is left-aligned from column 1
            @test startswith(lines[2], "edges ")
            @test startswith(lines[3], "gwesp.fixed.0.5")

            # z values are computed from est/se; p-value floor applies
            @test occursin("<1e-16", lines[2])
            # Significance codes per row and the legend below
            @test endswith(lines[2], "***")
            @test endswith(lines[3], "**")
            @test occursin("---", out)
            @test occursin(SIGNIF_LEGEND, out)

            # legend=false suppresses the legend but keeps the codes
            out2 = sprint(io -> print_coeftable(io, names, est, se, p; legend=false))
            @test !occursin("Signif. codes", out2)
            @test endswith(split(out2, '\n')[2], "***")

            # signif_codes=false suppresses codes and legend
            out3 = sprint(io -> print_coeftable(io, names, est, se, p;
                                                signif_codes=false))
            @test !occursin("***", out3)
            @test !occursin("Signif. codes", out3)

            # Symbols as names; custom header and explicit z (t) values
            out4 = sprint(io -> print_coeftable(io, [:density, :reciprocity],
                                                [0.1, 2.0], [0.05, 0.5],
                                                [0.0455, 6.3e-5];
                                                z_values=[2.0, 4.0],
                                                header=("Estimate", "Std.Error",
                                                        "t value", "Pr(>|t|)")))
            @test occursin("t value", out4)
            @test occursin("density", out4)
            @test occursin("2.0000", out4)
            @test occursin("6.3e-05", out4)

            # Length mismatch throws
            @test_throws ArgumentError print_coeftable(devnull, ["a"], [1.0, 2.0],
                                                       [0.1, 0.2], [0.5, 0.5])
            @test_throws ArgumentError print_coeftable(devnull, ["a", "b"],
                                                       [1.0, 2.0], [0.1], [0.5, 0.5])
            # --- optional Std.Error column + p_floor ---
            # A permutation test (SNA.netlm) has a statistic and a p-value but no
            # standard error: std_errors=nothing drops the column, z_values is
            # mandatory, and the header is a 3-tuple.
            out5 = sprint(io -> print_coeftable(io, ["intercept", "x"], [0.12, 0.85],
                                                nothing, [0.31, 0.0];
                                                z_values=[1.02, 6.7],
                                                header=("Estimate", "t value",
                                                        "Pr(>=|t|)"),
                                                p_floor=1e-3))
            lines5 = split(out5, '\n')
            @test occursin("t value", lines5[1]) && occursin("Pr(>=|t|)", lines5[1])
            @test !occursin("Std.Error", out5)
            @test !occursin("z value", out5)
            # p_floor = 1/reps: a zero permutation count prints "<0.001", never 0.0
            @test occursin("<0.001", lines5[3])
            @test !occursin("0.0000", lines5[3])
            @test endswith(lines5[3], "***")
            @test occursin("0.3100", lines5[2])
            # Three numeric columns, right-aligned under their titles
            for (title, values) in (("Estimate", ["0.1200", "0.8500"]),
                                    ("t value", ["1.0200", "6.7000"]),
                                    ("Pr(>=|t|)", ["0.3100", "<0.001"]))
                stop = findfirst(title, lines5[1]).stop
                for (i, val) in enumerate(values)
                    r = findfirst(val, lines5[1 + i])
                    @test r !== nothing && r.stop == stop
                end
            end
            # Default 3-column header
            out6 = sprint(io -> print_coeftable(io, ["a"], [1.0], nothing, [0.5];
                                                z_values=[0.7]))
            @test occursin("Estimate", out6) && occursin("z value", out6) &&
                  occursin("Pr(>|z|)", out6) && !occursin("Std.Error", out6)
            # Without std_errors, z_values must be given
            @test_throws ArgumentError print_coeftable(devnull, ["a"], [1.0], nothing, [0.5])
            # A 4-tuple header on a 3-column table (and vice versa) throws
            @test_throws ArgumentError print_coeftable(devnull, ["a"], [1.0], nothing, [0.5];
                                                       z_values=[0.7],
                                                       header=("Estimate", "Std.Error",
                                                               "z value", "Pr(>|z|)"))
            @test_throws ArgumentError print_coeftable(devnull, ["a"], [1.0], [0.1], [0.5];
                                                       header=("Estimate", "z value",
                                                               "Pr(>|z|)"))
            # Length mismatch without std_errors still throws
            @test_throws ArgumentError print_coeftable(devnull, ["a", "b"], [1.0, 2.0],
                                                       nothing, [0.5, 0.5]; z_values=[1.0])
            # p_floor applies to the 4-column table too, and the default floor is
            # unchanged for the existing callers (1e-16)
            out7 = sprint(io -> print_coeftable(io, ["a"], [1.0], [0.1], [0.0]; p_floor=1e-3))
            @test occursin("<0.001", out7) && endswith(split(out7, '\n')[2], "***")
            out8 = sprint(io -> print_coeftable(io, ["a"], [1.0], [0.1], [0.0]))
            @test occursin("<1e-16", out8)
            # Explicit header=nothing and z_values=nothing reproduce the defaults
            @test sprint(io -> print_coeftable(io, names, est, se, p;
                                               z_values=nothing, header=nothing)) == out
        end

        @testset "CoefficientTable / coeftable" begin
            names = ["edges", "gwesp.fixed.0.5", "mutual"]
            est = [-3.2145, 0.8113, 0.42]
            se = [0.1936, 0.2854, 0.0]          # the third SE is undefined
            tbl = CoefficientTable(names, est, se)
            @test tbl isa CoefficientTable
            @test tbl.names == names && tbl.estimates == est && tbl.std_errors == se
            @test length(tbl) == 3
            @test tbl.header == ("Estimate", "Std.Error", "z value", "Pr(>|z|)")
            @test tbl.p_floor == 1e-16
            # Defaults come from the ONE z → p helper (NaN where se ≤ 0)
            zp = z_pvalues(est, se)
            @test isequal(tbl.z_values, zp.z) && isequal(tbl.p_values, zp.p)
            @test isnan(tbl.z_values[3]) && isnan(tbl.p_values[3])
            @test tbl.z_values[1] == est[1] / se[1]
            # show delegates to print_coeftable on the same vectors: byte-identical
            expected = sprint(io -> print_coeftable(io, names, est, se, zp.p;
                                                    z_values=zp.z))
            @test sprint(show, tbl) == expected
            @test sprint(show, MIME("text/plain"), tbl) == expected
            @test repr(tbl) == expected
            @test occursin("<1e-16", expected) && occursin("NaN", expected)
            # Row access by index and by name (Symbol or String)
            row = tbl[2]
            @test row.name == "gwesp.fixed.0.5" && row.estimate == 0.8113 &&
                  row.std_error == 0.2854 && row.z_value == zp.z[2] &&
                  row.p_value == zp.p[2]
            @test tbl["mutual"].estimate == 0.42
            @test tbl[:edges] == tbl[1]
            @test isequal(tbl[end], tbl[3])
            @test_throws KeyError tbl["nope"]
            @test_throws KeyError tbl[:nope]
            @test_throws BoundsError tbl[0]
            @test_throws BoundsError tbl[4]
            # Iteration yields the rows in order
            @test isequal(collect(tbl), [tbl[1], tbl[2], tbl[3]])
            @test [r.name for r in tbl] == names
            @test all(r -> r isa NamedTuple && propertynames(r) ==
                          (:name, :estimate, :std_error, :z_value, :p_value), tbl)
            # Symbols as names are stringified
            @test CoefficientTable([:a, :b], [1.0, 2.0], [0.1, 0.2]).names == ["a", "b"]
            # Explicit z/p override the defaults; p only → z still from z_pvalues
            t2 = CoefficientTable(["a", "b"], [1.0, 2.0], [0.5, 0.5];
                                  z_values=[3.0, 4.0], p_values=[0.01, 0.02])
            @test t2.z_values == [3.0, 4.0] && t2.p_values == [0.01, 0.02]
            t3 = CoefficientTable(["a", "b"], [1.0, 2.0], [0.5, 0.5]; p_values=[0.01, 0.02])
            @test t3.z_values == [2.0, 4.0] && t3.p_values == [0.01, 0.02]
            t4 = CoefficientTable(["a", "b"], [1.0, 2.0], [0.5, 0.5]; z_values=[3.0, 4.0])
            @test t4.p_values == z_pvalues([3.0, 4.0])
            # Custom header (a t-based model) is stored and printed
            t5 = CoefficientTable(["a"], [1.0], [0.5];
                                  header=("Estimate", "Std.Error", "t value", "Pr(>|t|)"))
            @test occursin("t value", sprint(show, t5))
            # Without standard errors: 3 columns, z required, p defaults to z_pvalues(z)
            t6 = CoefficientTable(["intercept", "x"], [0.12, 0.85], nothing;
                                  z_values=[1.02, 6.7], p_values=[0.31, 0.0],
                                  header=("Estimate", "t value", "Pr(>=|t|)"),
                                  p_floor=1/1000)
            @test t6.std_errors === nothing && length(t6.header) == 3
            @test t6[1].std_error === nothing
            s6 = sprint(show, t6)
            @test s6 == sprint(io -> print_coeftable(io, t6.names, t6.estimates, nothing,
                                                     t6.p_values; z_values=t6.z_values,
                                                     header=t6.header, p_floor=1e-3))
            @test !occursin("Std.Error", s6) && occursin("<0.001", s6) &&
                  !occursin("0.0000 ", s6)
            t7 = CoefficientTable(["a"], [1.0], nothing; z_values=[2.0])
            @test t7.p_values == z_pvalues([2.0])
            @test t7.header == ("Estimate", "z value", "Pr(>|z|)")
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], nothing)
            # Length and header validation
            @test_throws ArgumentError CoefficientTable(["a", "b"], [1.0], [0.1, 0.2])
            @test_throws ArgumentError CoefficientTable(["a", "b"], [1.0, 2.0], [0.1])
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], [0.1]; z_values=[1.0, 2.0])
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], [0.1]; p_values=[0.1, 0.2])
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], [0.1];
                                                        header=("Estimate", "z", "p"))
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], nothing; z_values=[1.0],
                                                        header=("E", "S", "z", "p"))
            @test_throws ArgumentError CoefficientTable(["a"], [1.0], [0.1]; p_floor=2.0)

            # The StatsAPI verb: ONE binding, no NetworkCore method on it. A package
            # adds `StatsAPI.coeftable(::ItsType)` and gets `coeftable(fit)` for free.
            @test NetworkCore.coeftable === StatsAPI.coeftable
            @test :coeftable in Base.names(NetworkCore) && Base.isexported(NetworkCore, :coeftable)
            @test !any(m -> m.module === NetworkCore, methods(StatsAPI.coeftable))
            @test !hasmethod(coeftable, Tuple{CoefficientTable})
            fit = DummyFullFit([1.2, -0.3], [0.16 0.0; 0.0 0.01])
            ct = coeftable(fit)
            @test ct isa CoefficientTable
            @test ct.names == ["a", "b"] && ct.estimates == [1.2, -0.3]
            @test ct.std_errors == [0.4, 0.1]
            @test ct["a"].p_value ≈ z_pvalues([3.0])[1]
            @test sprint(show, ct) == sprint(show, fit)   # show IS the table
        end

        @testset "check_statsapi" begin
            full = DummyFullFit([1.2, -0.3], [0.16 0.0; 0.0 0.01])
            r = check_statsapi(full)
            @test r isa NamedTuple
            @test keys(r) == (:coef, :stderror, :vcov, :confint, :loglikelihood,
                              :nobs, :dof, :aic, :bic, :coeftable)
            @test all(values(r))
            @test all(check_statsapi(full))
            @test check_statsapi(full; strict=true) === r
            @test check_statsapi(full; strict=true) !== nothing
            # The subset form keeps the caller's order
            @test check_statsapi(full; required=(:coeftable, :coef)) ==
                  (coeftable=true, coef=true)
            @test check_statsapi(full; required=[:aic, :bic]) == (aic=true, bic=true)

            # Three of ten: the NamedTuple names exactly the gaps ...
            partial = DummyPartialFit([1.2, -0.3], [0.4, 0.1])
            rp = check_statsapi(partial)
            @test rp.coef && rp.stderror && rp.coeftable
            @test [k for (k, v) in pairs(rp) if !v] ==
                  [:vcov, :confint, :loglikelihood, :nobs, :dof, :aic, :bic]
            # ... and strict=true throws naming them, and nothing else
            err = try
                check_statsapi(partial; strict=true)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("DummyPartialFit", err.msg)
            for v in (:vcov, :confint, :loglikelihood, :nobs, :dof, :aic, :bic)
                @test occursin("no method StatsAPI.$v(::DummyPartialFit)", err.msg)
            end
            @test !occursin("StatsAPI.coef(", err.msg)
            @test !occursin("StatsAPI.stderror(", err.msg)
            @test !occursin("StatsAPI.coeftable(", err.msg)
            # Restricted to what it has, the partial fit passes strict
            @test check_statsapi(partial; required=(:coef, :stderror, :coeftable),
                                 strict=true) == (coef=true, stderror=true, coeftable=true)

            # Consistency checks on methods that exist: a stderror of the wrong
            # length, a vcov of the wrong shape, a coeftable that is not a
            # CoefficientTable, a coeftable with the wrong row count, and a method
            # that throws are all reported false with the reason in the message.
            bad_se = DummyPartialFit([1.2, -0.3], [0.4])
            rb = check_statsapi(bad_se; required=(:coef, :stderror))
            @test rb == (coef=true, stderror=false)
            eb = try check_statsapi(bad_se; required=(:coef, :stderror), strict=true)
                 catch e; e end
            @test eb isa ArgumentError && occursin("length 1", eb.msg) &&
                  occursin("length 2", eb.msg)
            bad_v = DummyFullFit([1.2, -0.3], [0.16 0.0 0.0; 0.0 0.01 0.0])
            @test check_statsapi(bad_v; required=(:coef, :vcov)) == (coef=true, vcov=false)
            bad_v2 = DummyFullFit([1.2, -0.3, 0.5], [0.16 0.0; 0.0 0.01])
            @test check_statsapi(bad_v2; required=(:coef, :vcov)) == (coef=true, vcov=false)
            @test check_statsapi(bad_v2; required=(:coef, :coeftable)) == (coef=true, coeftable=false)
            @test occursin("rows", (try check_statsapi(bad_v2; required=(:coef, :coeftable), strict=true)
                                    catch e; e end).msg)
            df_fit = DummyDataFrameCoeftable()
            @test check_statsapi(df_fit; required=(:coeftable,)) == (coeftable=false,)
            @test occursin("not a CoefficientTable",
                           (try check_statsapi(df_fit; required=(:coeftable,), strict=true)
                            catch e; e end).msg)
            throwing = DummyThrowingFit()
            @test check_statsapi(throwing; required=(:coef,)) == (coef=false,)
            @test occursin("threw", (try check_statsapi(throwing; required=(:coef,), strict=true)
                                     catch e; e end).msg)
            # A verb that is not a StatsAPI function, or an empty list, is refused
            @test_throws ArgumentError check_statsapi(full; required=(:coef, :nonsense))
            @test_throws ArgumentError check_statsapi(full; required=())
            # The surface has ten verbs required by default, and `coefnames`
            @test NetworkCore.STATSAPI_VERBS == (:coef, :stderror, :vcov, :confint,
                                              :loglikelihood, :nobs, :dof, :aic, :bic,
                                              :coeftable)
            @test NetworkCore.STATSAPI_OPTIONAL_VERBS == (:coefnames,)
            @test_throws ArgumentError check_statsapi(full; required=(:predict,))
        end

        @testset "check_statsapi: coefnames" begin
            # The verb is the StatsAPI binding (StatsBase and GLM export the same one)
            @test NetworkCore.coefnames === StatsAPI.coefnames
            @test Base.isexported(NetworkCore, :coefnames)
            @test !any(m -> m.module === NetworkCore, methods(StatsAPI.coefnames))
            full = DummyFullFit([1.2, -0.3], [0.16 0.0; 0.0 0.01])
            @test coefnames(full) == coeftable(full).names == ["a", "b"]
            @test check_statsapi(full; required=(:coef, :coeftable, :coefnames)) ==
                  (coef=true, coeftable=true, coefnames=true)
            @test check_statsapi(full; required=(:coefnames,), strict=true) ==
                  (coefnames=true,)
            # Not required by default: a fit without it still passes strict ...
            partial = DummyPartialFit([1.2, -0.3], [0.4, 0.1])
            @test !hasmethod(StatsAPI.coefnames, Tuple{DummyPartialFit})
            @test check_statsapi(partial; required=(:coef, :stderror, :coeftable),
                                 strict=true) == (coef=true, stderror=true, coeftable=true)
            # ... but fails when it is listed
            @test check_statsapi(partial; required=(:coef, :coefnames)) ==
                  (coef=true, coefnames=false)
            e0 = try check_statsapi(partial; required=(:coefnames,), strict=true)
                 catch e; e end
            @test e0 isa ArgumentError &&
                  occursin("no method StatsAPI.coefnames(::DummyPartialFit)", e0.msg)
            # Names that disagree with the table, or with the number of
            # coefficients, or that are not strings, are reported false
            swapped = DummyNamedFit([1.0, 2.0], ["b", "a"])
            @test check_statsapi(swapped; required=(:coef, :coefnames)) ==
                  (coef=true, coefnames=false)
            short = DummyNamedFit([1.0, 2.0], ["a"])
            @test check_statsapi(short; required=(:coefnames,)) == (coefnames=false,)
            es = try check_statsapi(short; required=(:coefnames,), strict=true)
                 catch e; e end
            @test es isa ArgumentError && occursin("length 1", es.msg) &&
                  occursin("coefnames", es.msg)
            syms = DummyNamedFit([1.0, 2.0], [:a, :b])
            @test check_statsapi(syms; required=(:coefnames,)) == (coefnames=false,)
            good = DummyNamedFit([1.0, 2.0], ["a", "b"])
            @test check_statsapi(good; required=(:coef, :coeftable, :coefnames),
                                 strict=true) == (coef=true, coeftable=true, coefnames=true)
            # Implemented but not listed: still checked, and strict throws naming it
            @test check_statsapi(swapped; required=(:coef,)) == (coef=true,)
            ew = try check_statsapi(swapped; required=(:coef,), strict=true)
                 catch e; e end
            @test ew isa ArgumentError && occursin("missing or inconsistent: coefnames", ew.msg) &&
                  occursin("differs from coeftable", ew.msg)
            @test check_statsapi(good; required=(:coef,), strict=true) == (coef=true,)
        end

        @testset "mc_pvalue" begin
            # Never exactly zero: (1 + k)/(N + 1) estimator
            @test mc_pvalue(zeros(99), 100.0) == 2 * 1 / 100
            @test mc_pvalue([1.0, 2.0, 3.0, 4.0], 5.0) == 2 * (1 + 0) / (4 + 1)
            # Observed in the middle of the distribution: capped at 1
            @test mc_pvalue([1.0, 2.0, 3.0], 2.0) == 1.0
            @test_throws ArgumentError mc_pvalue(Float64[], 1.0)
            # An undefined statistic is never "significant": a NaN
            # observed value, or no defined draw, gives NaN — it used to give
            # 2/(N+1) = 0.02 for N = 99.
            @test isnan(mc_pvalue(randn(Xoshiro(1), 99), NaN))
            @test isnan(mc_pvalue(fill(NaN, 99), 1.0))
            # NaN draws are dropped and N counts the defined ones
            @test mc_pvalue([1.0, NaN, 3.0, 4.0], 5.0) == 2 * (1 + 0) / (3 + 1)
            @test mc_pvalue([1.0, NaN, 2.0, NaN, 3.0], 2.0) == mc_pvalue([1.0, 2.0, 3.0], 2.0)
            # brute force: NaN-padded draws ≡ the defined draws alone
            rng = Xoshiro(7)
            for _ in 1:200
                sim = randn(rng, rand(rng, 1:30))
                obs = randn(rng)
                padded = vcat(sim, fill(NaN, rand(rng, 0:10)))
                @test mc_pvalue(shuffle(rng, padded), obs) == mc_pvalue(sim, obs)
            end
            # GOFStatistic defaults and display inherit it: NaN p-value, NaN
            # envelope, no "quantiles undefined in presence of NaNs" error
            st = GOFStatistic("x", ["a", "b"], [1.0, NaN], [1.0 NaN; 2.0 NaN; NaN NaN])
            @test st.p_values[1] == mc_pvalue([1.0, 2.0], 1.0)
            @test isnan(st.p_values[2])
            out = sprint(show, st)
            @test occursin("NaN", out) && occursin("1.0000", out)
        end

        @testset "GOFStatistic" begin
            sims = [1.0 2.0; 3.0 4.0; 5.0 6.0]   # 3 sims × 2 levels
            stat = GOFStatistic("degree", ["0", "1"], [3, 10], sims)
            @test stat.name == "degree"
            @test stat.labels == ["0", "1"]
            @test stat.observed == [3.0, 10.0]
            # Default p-values via mc_pvalue
            @test stat.p_values[1] == mc_pvalue(sims[:, 1], 3.0)
            @test stat.p_values[2] == mc_pvalue(sims[:, 2], 10.0)
            # Explicit p-values are taken as given
            stat2 = GOFStatistic("degree", [0, 1], [3, 10], sims;
                                 p_values=[0.5, 0.001])
            @test stat2.labels == ["0", "1"]   # non-string labels stringified
            @test stat2.p_values == [0.5, 0.001]
            # Dimension validation
            @test_throws ArgumentError GOFStatistic("d", ["0"], [3, 10], sims)
            @test_throws ArgumentError GOFStatistic("d", ["0", "1"], [3, 10],
                                                    zeros(3, 1))
            @test_throws ArgumentError GOFStatistic("d", ["0", "1"], [3, 10],
                                                    zeros(0, 2))
            @test_throws ArgumentError GOFStatistic("d", ["0", "1"], [3, 10], sims;
                                                    p_values=[0.5])
        end

        @testset "GOFResult show" begin
            sims = [0.0 10.0; 2.0 12.0; 4.0 14.0; 6.0 16.0]   # 4 sims × 2 levels
            deg = GOFStatistic("degree", ["deg0", "deg1"], [3, 100], sims)
            esp = GOFStatistic("esp", ["esp0"], [2],
                               reshape([1.0, 2.0, 3.0, 4.0], 4, 1))
            result = GOFResult([deg, esp]; model="ERGM", p_overall=0.42)
            @test n_simulations(result) == 4

            out = sprint(show, result)
            lines = split(out, '\n')
            @test occursin("Goodness-of-fit assessment: ERGM", out)
            @test occursin("Based on 4 simulated networks", out)
            @test occursin("Overall p-value: 0.4200", out)
            # One table heading per statistic category
            @test occursin("Goodness-of-fit for degree", out)
            @test occursin("Goodness-of-fit for esp", out)
            # Observed vs simulated quantile columns
            for h in ("Obs", "2.5%", "Mean", "97.5%", "MC p-value")
                @test occursin(h, out)
            end
            # One row per level, observed values shown compactly
            degrow = lines[findfirst(l -> startswith(l, "deg0"), lines)]
            @test occursin(r"deg0\s+3\s", degrow * " ")
            @test occursin("3.00", degrow)   # simulated mean of level 1
            # An observed value far outside the envelope gets the smallest
            # attainable MC p-value, never 0
            deg1row = lines[findfirst(l -> startswith(l, "deg1"), lines)]
            @test occursin("0.4000", deg1row)   # 2*(1+0)/(4+1)
            @test !occursin(" 0.0000", deg1row)

            # Column alignment within a GOF table
            hline = lines[findfirst(l -> occursin("MC p-value", l), lines)]
            stop = findfirst("Mean", hline).stop
            r = findfirst("3.00", degrow)
            @test r !== nothing && r.stop == stop

            # Single-statistic convenience constructor, no model/p_overall
            r1 = GOFResult(deg)
            out1 = sprint(show, r1)
            @test occursin("Goodness-of-fit assessment", out1)
            @test !occursin("Overall p-value", out1)
            # GOFStatistic has its own show
            @test occursin("Goodness-of-fit for degree", sprint(show, deg))

            @test_throws ArgumentError GOFResult(GOFStatistic[])
        end

        @testset "gof generic" begin
            # Empty generic: exported, no methods; packages add their own
            @test length(methods(gof)) == 0
            @test_throws MethodError gof(1)
        end
    end

    @testset "Shared statistic protocol" begin
        # `compute` / `name` / `compute_all` live here for the same reason `gof`
        # does. ERGM.jl and REM.jl each used to define and export their OWN
        # `compute` and `name`; two different functions with one name, so
        # `using ERGM, REM` — cross-sections plus dynamics, the core statnet
        # workflow — left both verbs UNDEFINED by Julia's conflicting-export
        # rule. One generic, extended by everyone, cannot collide.

        @testset "empty generics, exported" begin
            for f in (compute, name, compute_all)
                @test f isa Function
                @test length(methods(f)) == 0    # NetworkCore defines no methods
                @test_throws MethodError f(1)
            end
            for s in (:compute, :name, :compute_all)
                @test s in names(NetworkCore)
            end
        end

        @testset "one generic carries incompatible signatures" begin
            # The property the fix rests on: the model packages' methods differ
            # in signature, so they coexist on one function. Simulated here with
            # two term types standing in for an ERGM term and a REM statistic.
            struct FakeTerm end
            struct FakeStat end
            NetworkCore.compute(::FakeTerm, net) = 1.0
            NetworkCore.compute(::FakeStat, state, sender::Int, receiver::Int) = 2.0
            NetworkCore.name(::FakeTerm) = "term"
            NetworkCore.name(::FakeStat) = "stat"

            @test compute(FakeTerm(), nothing) == 1.0
            @test compute(FakeStat(), nothing, 1, 2) == 2.0
            @test name(FakeTerm()) == "term"
            @test name(FakeStat()) == "stat"
        end
    end

    @testset "Module/type names" begin
        # The module is `NetworkCore`, the type is `Network`. Until v0.2 they shared
        # a name, so the exported type shadowed the module and all three of these
        # were broken. They are the whole point of the rename — pin them.

        @testset "the documented constructor is callable" begin
            # Used to throw `MethodError: objects of type Module are not callable`
            net = Network(5)
            @test net isa Network{Int,true}
            @test nv(net) == 5
            @test Network(4; directed=false) isa Network{Int,false}
            # The lowercase alias keeps working — it was the only thing that did
            @test network(5) isa Network{Int,true}
            @test typeof(network(5)) === typeof(Network(5))
        end

        @testset "qualified module access resolves to the module" begin
            # `Network.supports_missing` used to be parsed as FIELD ACCESS on the
            # exported type (a UnionAll) and failed to precompile in every
            # dependent package. Qualified access must now reach the module.
            @test NetworkCore.supports_missing(sum) == false
            @test NetworkCore.MISSING_POLICIES == (:error, :face)
            @test NetworkCore.is_directed(Network{Int,false}) == false
            @test isa(NetworkCore, Module)
            @test NetworkCore.Network === Network
        end

        @testset "no Base.Docs workaround is needed" begin
            # The collision forced a custom `Base.Docs.doc!(::Type{Network}, ...)`
            # method that reached into Base internals, and precompilation warned
            # that docs for `Network.Network` were being replaced. With the module
            # renamed, neither exists: the type's docstring is its own again.
            @test isempty(methods(Base.Docs.doc!, (Type{Network}, Base.Docs.Binding,
                                                   Base.Docs.DocStr)))
            @test !isempty(string(@doc Network))
        end
    end

    @testset "Golden-fixture harness" begin
        # A frozen reference value is only evidence if it can be REGENERATED.
        # The harness enforces that: a fixture without provenance won't load.

        fixdir = joinpath(@__DIR__, "fixtures")

        @testset "loads a provenanced fixture" begin
            g = load_golden(joinpath(fixdir, "florentine_sna.toml"))
            @test g isa GoldenFixture
            @test g.name == "florentine_sna"
            # Provenance is what makes the number auditable rather than folklore
            @test haskey(g.provenance, "r_version")
            @test haskey(g.provenance, "seed")
            @test haskey(g.provenance, "script")
            @test g.provenance["sna_version"] == "2.8"
            # The generating script was resolved on disk, to an absolute path
            @test g.script_path == joinpath(pkgdir(NetworkCore), "test", "fixtures",
                                            "r", "harness_selftest.R")
            @test isabspath(g.script_path)
            @test isfile(g.script_path)
            @test occursin(g.script_path, sprint(show, g))
        end

        @testset "refuses a fixture it cannot audit" begin
            mktempdir() do dir
                # No [provenance] at all
                bare = joinpath(dir, "bare.toml")
                write(bare, "[values]\nx = 1.0\n")
                @test_throws ArgumentError load_golden(bare)

                # Provenance present but incomplete — no seed, so a Monte-Carlo
                # fixture could never be reproduced
                partial = joinpath(dir, "partial.toml")
                write(partial, """
                    [provenance]
                    r_version = "4.6.1"
                    script = "x.R"

                    [values]
                    x = 1.0
                    """)
                @test_throws ArgumentError load_golden(partial)

                # Values missing entirely
                novals = joinpath(dir, "novals.toml")
                write(novals, """
                    [provenance]
                    r_version = "4.6.1"
                    seed = 1
                    script = "x.R"
                    """)
                @test_throws ArgumentError load_golden(novals)

                @test_throws ArgumentError load_golden(joinpath(dir, "absent.toml"))

                # Complete provenance naming a script that does not exist: the
                # fixture cannot be regenerated, so it does not load either.
                ghost = joinpath(dir, "ghost.toml")
                write(ghost, """
                    [provenance]
                    r_version = "4.6.1"
                    seed = 1
                    script = "r/ghost.R"

                    [values]
                    x = 1.0
                    """)
                @test_throws ArgumentError load_golden(ghost)
                msg = try
                    load_golden(ghost)
                catch e
                    sprint(showerror, e)
                end
                # The message names the fixture, the script string and the
                # directories that were searched (the fixture's own directory
                # is the first ancestor tried)
                @test occursin(ghost, msg)
                @test occursin("\"r/ghost.R\"", msg)
                @test occursin(dir, msg)
                @test occursin("searched:", msg)
                @test occursin("root=", msg)

                # ...and loads as soon as the script is where it says it is
                # (resolved against the fixture's directory, the first
                # ancestor searched)
                mkpath(joinpath(dir, "r"))
                write(joinpath(dir, "r", "ghost.R"), "# placeholder\n")
                g = load_golden(ghost)
                @test g isa GoldenFixture
                @test g.script_path == joinpath(dir, "r", "ghost.R")
                @test occursin(g.script_path, golden_report(g, "x", 1.0))
            end
        end

        @testset "resolves provenance.script by walking up the ancestors" begin
            mktempdir() do dir
                # Script two levels above the fixture, under the conventional
                # package-relative path — found by the walk, no Project.toml
                # needed
                mkpath(joinpath(dir, "test", "fixtures", "r"))
                write(joinpath(dir, "test", "fixtures", "r", "deep.R"), "# r\n")
                fixture = joinpath(dir, "test", "fixtures", "deep.toml")
                write(fixture, """
                    [provenance]
                    r_version = "4.6.1"
                    seed = 1
                    script = "test/fixtures/r/deep.R"

                    [values]
                    x = 1.0
                    """)
                g = load_golden(fixture)
                @test g.script_path == joinpath(dir, "test", "fixtures", "r", "deep.R")

                # A fixture copied elsewhere cannot find it from its own
                # ancestors, but `root=` says where the path is relative to
                mkpath(joinpath(dir, "elsewhere"))
                moved = joinpath(dir, "elsewhere", "deep.toml")
                cp(fixture, moved)
                @test_throws ArgumentError cd(homedir()) do
                    load_golden(moved; root=homedir())
                end
                g2 = load_golden(moved; root=dir)
                @test g2.script_path == joinpath(dir, "test", "fixtures", "r", "deep.R")
                # Ancestors of `root` are searched too
                g3 = load_golden(moved; root=joinpath(dir, "elsewhere"))
                @test g3.script_path == g2.script_path

                # An absolute script path is used as given
                absfix = joinpath(dir, "abs.toml")
                # TOML must escape the backslashes in a native Windows path.
                absolute_fixture = TOML.parsefile(fixture)
                absolute_fixture["provenance"]["script"] =
                    joinpath(dir, "test", "fixtures", "r", "deep.R")
                open(absfix, "w") do io
                    TOML.print(io, absolute_fixture)
                end
                @test load_golden(absfix).script_path ==
                      joinpath(dir, "test", "fixtures", "r", "deep.R")
            end
        end

        @testset "resolves provenance.script against the package root" begin
            # The conventional `script = "test/fixtures/r/<name>.R"` is relative
            # to the PACKAGE root, not to the working directory, so a fixture
            # must load from anywhere.
            g = cd(homedir()) do
                load_golden(joinpath(fixdir, "florentine_sna.toml"))
            end
            @test g isa GoldenFixture
            @test isfile(joinpath(pkgdir(NetworkCore), g.provenance["script"]))
            @test g.script_path == normpath(joinpath(pkgdir(NetworkCore), g.provenance["script"]))
        end

        @testset "checks values against the declared tolerance" begin
            g = load_golden(joinpath(fixdir, "florentine_sna.toml"))

            # The tolerance is read from the fixture, not invented at the call site
            @test golden_tolerance(g, "density") == 1e-10
            @test golden_tolerance(g, "anything_else") == 1e-10   # the default

            @test check_golden(g, "density", 1 / 6)
            @test !check_golden(g, "density", 0.2)
            @test check_golden(g, "degree", g.values["degree"])
            @test !check_golden(g, "degree", [0, 0, 0])           # length mismatch

            @test_throws ArgumentError check_golden(g, "no_such_value", 1.0)

            # A failing comparison must explain itself, provenance included
            rep = golden_report(g, "density", 0.2)
            @test occursin("reference (R)", rep)
            @test occursin("observed (jl)", rep)
            @test occursin("atol", rep)
            @test occursin("r_version=4.6.1", rep)
            @test occursin("harness_selftest.R", rep)
        end

        @testset "NetworkCore reproduces R sna on the Florentine data" begin
            # The harness proving itself on real frozen R output.
            g = load_golden(joinpath(fixdir, "florentine_sna.toml"))
            flo = load_dataset(:florentine_marriage)

            @test nv(flo) == length(g.values["degree"])
            # NetworkCore's own density must match sna::gden on the same graph
            @test check_golden(g, "density", Graphs.density(flo)) ||
                  error(golden_report(g, "density", Graphs.density(flo)))
        end
    end

    @testset "Result-metadata protocol" begin
        # Six issues across the ecosystem ask the same question in different
        # words: what did this fit ACTUALLY do? The protocol makes the answer
        # inspectable instead of a sentence in a `show` method.

        @testset "defaults are conservative, not optimistic" begin
            struct UnopedInResult end
            r = UnopedInResult()

            # A package that has not declared must NOT be assumed exact
            @test is_exact(r) == false
            @test estimand(r) == :unspecified
            @test objective(r) == :unspecified
            @test se_method(r) == :unspecified
            @test missing_method(r) == :unspecified
            @test tie_method(r) == :not_applicable
            @test approximations(r) == String[]
        end

        @testset "metadata collects the protocol" begin
            struct FakeFit end
            NetworkCore.estimand(::FakeFit) = :stergm_formation
            NetworkCore.objective(::FakeFit) = :conditional_pseudolikelihood
            NetworkCore.is_exact(::FakeFit) = false
            NetworkCore.se_method(::FakeFit) = :bootstrap
            NetworkCore.missing_method(::FakeFit) = :available_case
            NetworkCore.approximations(::FakeFit) = ["CMPLE of a dyad-dependent formula"]

            md = fit_metadata(FakeFit())
            @test md isa ResultMetadata
            @test md.estimand == :stergm_formation
            @test md.objective == :conditional_pseudolikelihood
            @test md.is_exact == false
            @test md.se_method == :bootstrap
            @test md.missing_method == :available_case
            @test md.approximations == ["CMPLE of a dyad-dependent formula"]

            out = sprint(show, md)
            @test occursin("conditional_pseudolikelihood", out)
            @test occursin("approximate", out)          # is_exact == false
            @test occursin("bootstrap", out)
            @test occursin("CMPLE of a dyad-dependent formula", out)
        end

        @testset "is_exact is a property of the fit, not the estimator" begin
            # The whole point: the SAME estimator is exact for one model and an
            # approximation for another. CMPLE of a dyad-independent formula IS
            # the CMLE; of a dependent one it is not.
            struct IndepFit end
            struct DepFit end
            NetworkCore.objective(::Union{IndepFit,DepFit}) = :conditional_pseudolikelihood
            NetworkCore.is_exact(::IndepFit) = true
            NetworkCore.is_exact(::DepFit) = false

            @test objective(IndepFit()) == objective(DepFit())
            @test is_exact(IndepFit()) && !is_exact(DepFit())
            @test occursin("exact for this model", sprint(show, fit_metadata(IndepFit())))
            @test occursin("approximate", sprint(show, fit_metadata(DepFit())))
        end
    end

    @testset "Tied-event contract: TIE_POLICIES and check_tie_policy" begin
        # The sibling of the missing-data contract, for the models with an event
        # clock (REM.jl, Revel.jl). ONE vocabulary, defined
        # here, so `ties=:efron` cannot mean two things in two packages — and ONE
        # guard, so a policy a model cannot honour FAILS rather than no-ops.
        @test TIE_POLICIES == (:error, :ordered, :breslow, :efron, :batch)

        # A supported policy passes straight through
        for t in (:error, :ordered, :breslow)
            @test check_tie_policy(t, (:error, :ordered, :breslow);
                                   model="a model") === t
        end

        # A policy that IS in the vocabulary but NOT defined for this model is
        # refused — and the refusal says why, in the model's own terms
        e = try
            check_tie_policy(:batch, (:error, :ordered, :breslow);
                             model="the ordinal likelihood",
                             reasons=Dict(:batch => "a frozen batch IS Breslow here"))
            nothing
        catch err
            err
        end
        @test e isa ArgumentError
        m = sprint(showerror, e)
        @test occursin("`:batch` is not defined for the ordinal likelihood", m)
        @test occursin("a frozen batch IS Breslow here", m)   # the REASON
        @test occursin("`:breslow`", m)                       # ... and the way out

        # Without a reason it still refuses (silence is never the fallback)
        @test_throws ArgumentError check_tie_policy(:efron, (:error, :ordered);
                                                    model="a timing model")

        # A symbol outside the vocabulary altogether names the vocabulary
        e2 = try
            check_tie_policy(:sort_it_out, (:error, :ordered); model="a model")
            nothing
        catch err
            err
        end
        @test e2 isa ArgumentError
        m2 = sprint(showerror, e2)
        @test occursin("unknown tie policy `:sort_it_out`", m2)
        @test occursin("TIE_POLICIES", m2)
        @test all(occursin("`:$p`", m2) for p in TIE_POLICIES)

        # `tie_method` reports a policy that RAN, so `:error` can never be one:
        # under it a tie throws instead of producing a fit
        struct TiedFit end
        NetworkCore.tie_method(::TiedFit) = :efron
        @test fit_metadata(TiedFit()).tie_method === :efron
        @test occursin("tied events:    efron", sprint(show, fit_metadata(TiedFit())))
        # ... and a model without an event clock says so, rather than :none
        @test fit_metadata(IndepFit()).tie_method === :not_applicable
        @test !occursin("tied events", sprint(show, fit_metadata(IndepFit())))
    end

    @testset "Missing-data contract" begin
        # A masked dyad is UNOBSERVED, not absent. These tests pin the two
        # halves of the ecosystem contract: the trait routines opt into, and
        # the guard that stops face values being read by accident.

        @testset "supports_missing trait" begin
            # Default false: a routine that never declared a missing-data
            # method is one that must not be handed a masked network.
            @test supports_missing(sum) == false
            @test supports_missing(nothing) == false
        end

        @testset "missing_policies vocabulary" begin
            # The routine-level vocabulary: what a `missing=` keyword ACTUALLY
            # takes, so tooling prints the truth instead of assuming `:face`.
            @test missing_policies(sum) == (:error,)          # no keyword at all
            @test missing_policies(nothing) == (:error,)
            @test missing_policies(network_density) == (:error, :face)
            @test missing_policies(network_density) isa Tuple{Vararg{Symbol}}
            @test first(missing_policies(network_density)) === :error
            @test all(p in MISSING_POLICIES for p in missing_policies(network_density))

            # Introspection: every exported routine with a `missing=` keyword
            # must declare more than the bare default — otherwise the keyword
            # exists but the vocabulary claims it does not (or vice versa).
            declared = Symbol[]
            for nm in names(NetworkCore)
                f = getfield(NetworkCore, nm)
                f isa Function || continue
                has_kw = any(:missing in Base.kwarg_decl(m) for m in methods(f, NetworkCore))
                has_kw || continue
                push!(declared, nm)
                @test length(missing_policies(f)) >= 2
            end
            @test :network_density in declared

            # The per-result-type form for SHARED generics: `gof` has one method
            # per package, so a package declares `missing_policies(gof,
            # ::Type{<:ItsResult})` and tooling queries the two-argument form;
            # the fallback answers the one-argument form for everything else.
            @test missing_policies(network_density, Network) == (:error, :face)
            @test missing_policies(sum, Int) == (:error,)
            @test missing_policies(gof) == (:error,)
            @test missing_policies(gof, Nothing) == (:error,)
            struct MPResult end
            NetworkCore.missing_policies(::typeof(gof), ::Type{MPResult}) =
                (:error, :condition_on_face)
            @test missing_policies(gof, MPResult) == (:error, :condition_on_face)
            @test missing_policies(gof) == (:error,)     # the generic itself unchanged
        end

        @testset "require_observed guard" begin
            net = network(4)
            add_edge!(net, 1, 2)

            # No mask: passes under every policy, returns the network
            @test require_observed(net) === net
            @test require_observed(net, :face) === net

            set_missing_dyad!(net, 3, 4)
            @test_throws ArgumentError require_observed(net)
            @test_throws ArgumentError require_observed(net, :error)

            # :face is the explicit, auditable opt-in to face values
            @test require_observed(net, :face) === net

            # Unknown policies are rejected rather than silently ignored
            @test_throws ArgumentError require_observed(net, :ignore)

            # The message names the caller and the number of masked dyads
            msg = try
                require_observed(net; context="my_measure")
            catch e
                sprint(showerror, e)
            end
            @test occursin("my_measure", msg)
            @test occursin("1 masked dyad", msg)
            @test occursin("missing=:face", msg)     # this caller HAS the keyword

            # ...but most estimators do NOT expose a `missing=` keyword, and must
            # not: for TERGM/ERGMCount/ERGMMulti a `:face` fit would enumerate
            # unobserved dyads as observed — the very thing this contract stops.
            # `face_ok=false` drops the suggestion, because telling a user to pass
            # a keyword that does not exist is its own small lie.
            msg2 = try
                require_observed(net; context="fit_ergm_count", face_ok=false)
            catch e
                sprint(showerror, e)
            end
            @test occursin("fit_ergm_count", msg2)
            @test occursin("1 masked dyad", msg2)
            @test !occursin("missing=:face", msg2)   # not advertised
            @test occursin("clear_missing_dyads!", msg2)   # the escape that DOES exist

            # `hint=`: a routine with its OWN opt-in vocabulary names it in the
            # same message, after the shared bullets (ERGM's `:condition_on_face`
            # / `:mle` used to be glued on by catch-and-rethrow).
            msg3 = try
                require_observed(net; context="mcmle", face_ok=false,
                                 hint="  • pass `missing=:condition_on_face` to hold masked dyads fixed.\n")
            catch e
                sprint(showerror, e)
            end
            @test occursin("mcmle does not support", msg3)
            @test occursin("supports_missing(f) == true`).\n  • pass `missing=:condition_on_face`", msg3)
            @test !occursin("missing=:face", msg3)
            # No hint, no trailing newline: the message is unchanged
            @test !endswith(msg2, "\n")
            # A clean network passes regardless of the hint
            @test require_observed(network(3); hint="  • whatever\n") isa Network
        end

        @testset "as_matrix missing_value" begin
            net = network(4)
            add_edge!(net, 1, 2)
            set_missing_dyad!(net, 3, 4)

            # Default: face values only, mask invisible (the lossy path)
            @test as_matrix(net) isa Matrix{Float64}

            A = as_matrix(net; missing_value=missing)
            @test ismissing(A[3, 4])
            @test A[1, 2] == 1.0
            @test count(ismissing, A) == 1

            # Any other sentinel is written literally
            @test as_matrix(net; missing_value=-1.0)[3, 4] == -1.0

            # Sparse + missing_value cannot be represented; say so, don't guess
            @test_throws ArgumentError as_matrix(net; sparse=true, missing_value=missing)
        end

        @testset "as_matrix return_mask" begin
            net = network(4)
            add_edge!(net, 1, 2)
            set_missing_dyad!(net, 3, 4)

            A, mask = as_matrix(net; return_mask=true)
            @test A isa Matrix{Float64}
            @test mask isa BitMatrix
            @test mask[3, 4]
            @test count(mask) == 1
            @test A[1, 2] == 1.0        # face values untouched

            # An undirected mask is symmetric in the matrix, as it is in the net
            u = network(4; directed=false)
            set_missing_dyad!(u, 3, 4)
            _, umask = as_matrix(u; return_mask=true)
            @test umask[3, 4] && umask[4, 3]
        end

        @testset "round-trip preserves the mask" begin
            # Absent-face AND present-face masked dyads: the present-face case
            # is the one that distinguishes the two conversion paths.
            net = network(4)
            add_edge!(net, 1, 2)
            add_edge!(net, 2, 3)
            set_missing_dyad!(net, 1, 2)     # present face
            set_missing_dyad!(net, 3, 4)     # absent face

            # missing_value: mask survives, face value does not (documented)
            m = network_from_matrix(as_matrix(net; missing_value=missing))
            @test n_missing_dyads(m) == 2
            @test is_missing_dyad(m, 1, 2) && is_missing_dyad(m, 3, 4)
            @test !has_edge(m, 1, 2)         # face value was overwritten

            # missing_mask: fully lossless — mask AND face value survive
            A, mask = as_matrix(net; return_mask=true)
            r = network_from_matrix(A; missing_mask=mask, directed=true)
            @test n_missing_dyads(r) == 2
            @test is_missing_dyad(r, 1, 2) && is_missing_dyad(r, 3, 4)
            @test has_edge(r, 1, 2)          # present face preserved
            @test ne(r) == ne(net)

            @test_throws ArgumentError network_from_matrix(A; missing_mask=falses(2, 2))
        end

        @testset "round-trip: undirected" begin
            u = network(4; directed=false)
            add_edge!(u, 1, 2)
            set_missing_dyad!(u, 1, 2)      # present face
            set_missing_dyad!(u, 3, 4)      # absent face

            A, mask = as_matrix(u; return_mask=true)
            r = network_from_matrix(A; missing_mask=mask, directed=false)
            @test n_missing_dyads(r) == 2
            @test is_missing_dyad(r, 4, 3)  # symmetric lookup
            @test has_edge(r, 1, 2)
            @test ne(r) == ne(u)
        end

        @testset "round-trip: no mask is a no-op" begin
            net = network(4)
            add_edge!(net, 1, 2)
            r = network_from_matrix(as_matrix(net; missing_value=missing))
            @test n_missing_dyads(r) == 0
            @test ne(r) == 1
        end
    end

    # The ecosystem conversion contract (src/conversion.jl). NetworkCore.jl owns
    # the vocabulary; the adapters (DynamicNetworks, TSNA, REM, Siena) build the
    # reports. The per-path invariant table is docs/src/guide/conversion_invariants.md.
    @testset "Conversion contract: ConversionReport" begin
        lossless = ConversionReport(:Network, :DynamicNetwork)
        @test is_lossless(lossless)
        @test isempty(dropped_fields(lossless))
        @test occursin("lossless", sprint(show, lossless))
        @test occursin("lossless", sprint(show, MIME"text/plain"(), lossless))

        lossy = ConversionReport(:DynamicNetwork, :Network)
        record_drop!(lossy, :spells, "a static network has no time axis")
        record_drop!(lossy, :observation_period, "no static counterpart")
        @test !is_lossless(lossy)
        @test dropped_fields(lossy) == [:spells, :observation_period]
        @test lossy.source === :DynamicNetwork && lossy.target === :Network
        # record_drop! returns the report, so it chains
        @test record_drop!(lossy, :x, "y") === lossy
        @test length(lossy.dropped) == 3

        plain = sprint(show, MIME"text/plain"(), lossy)
        @test occursin("DynamicNetwork → Network", plain)
        @test occursin("spells", plain)
        @test occursin("a static network has no time axis", plain)
        @test occursin("3 field(s)", plain)
        @test occursin("3 dropped", sprint(show, lossy))
    end
    @testset "Exported API inventory" begin
        # The export list is FROZEN. Adding or removing an export is a deliberate
        # act: ERGM.jl and SNA.jl re-export a curated copy of this list
        # (minus the golden harness, `bootstrap_cov`, `record_drop!` and the
        # module name), and downstream `using ERGM; Network(5)` breaks silently
        # when the two drift. Update the literal below together with CHANGELOG.
        #
        # `names(M)` lists exported AND `public` names (Julia ≥ 1.11), and the
        # dependants' re-export loops iterate it, so both lists are frozen:
        # `frozen` is the exported surface, `frozen_public` the supported-but-
        # unexported one: the shared optimizer (re-exported by ERGM.jl), the
        # separation verdict and its policy helpers, and the task helpers,
        # all tooling for package authors rather than for a user's namespace.
        frozen = [
            :AbstractNetwork, :BipartiteNetwork, :CoefficientTable, :ConversionReport, :GOFResult, :GOFStatistic, :GoldenFixture,
            :MISSING_POLICIES, :Network, :NetworkCore, :ResultMetadata, :SIGNIF_LEGEND, :TIE_POLICIES,
            :add_edge!, :add_edges!, :add_vertex!, :add_vertices!, :approximations, :as_adjacency_matrix,
            :as_dataframe, :as_edgelist, :as_matrix, :bootstrap_cov, :check_golden, :check_se, :check_statsapi, :check_tie_policy,
            :clear_missing_dyads!, :coefnames, :coeftable, :compute, :compute_all, :degree, :delete_edge_attribute!, :delete_missing_dyad!,
            :delete_network_attribute!, :delete_vertex_attribute!, :dropped_fields, :dst, :edges, :estimand, :fit_metadata,
            :format_pvalue, :get_edge_attribute, :get_induced_subgraph, :get_neighborhood, :get_network_attribute, :get_vertex_attribute,
            :gof, :golden_report, :golden_tolerance, :has_edge, :has_vertex, :indegree,
            :inneighbors, :is_directed, :is_exact, :is_lossless, :is_missing_dyad, :is_two_mode,
            :list_edge_attributes, :list_network_attributes, :list_vertex_attributes, :load_dataset, :load_golden, :mc_pvalue,
            :missing_dyads, :missing_method, :missing_policies, :n_missing_dyads, :n_simulations, :name,
            :ne, :neighbors, :network, :network_density, :network_edgecount, :network_from_dataframe,
            :network_from_edgelist, :network_from_matrix, :network_initialize, :network_size, :nv, :objective,
            :outdegree, :outneighbors, :permute_vertices, :print_coeftable, :read_pajek, :record_drop!,
            :rem_edge!, :rem_vertex!, :require_observed, :se_method, :set_edge_attribute!, :set_missing_dyad!,
            :set_network_attribute!, :set_vertex_attribute!, :signif_code, :src, :supports_missing, :tie_method, :vertex_attribute_vector,
            :vertices, :write_edgelist_csv, :write_graphml, :write_pajek, :z_pvalues,
        ]
        frozen_public = [:SeparationVerdict, :clogit_separation, :logistic_derivatives,
                         :logistic_separation, :newton_fit, :poisson_separation,
                         :separation_caveat, :separation_from_margins, :separation_message,
                         :spawn_all, :unwrap_task_exception, :warn_separation]
        @test issorted(frozen) && issorted(frozen_public)
        all_names = names(NetworkCore)
        actual = sort(filter(nm -> Base.isexported(NetworkCore, nm), all_names))
        actual_public = sort(filter(nm -> !Base.isexported(NetworkCore, nm), all_names))
        @test actual == frozen
        # Name the drift when it happens, rather than dumping two 100-element vectors
        @test isempty(setdiff(actual, frozen))   # an unplanned NEW export
        @test isempty(setdiff(frozen, actual))   # an unplanned REMOVED export
        @test actual_public == frozen_public
        @test all(nm -> Base.ispublic(NetworkCore, nm) && !Base.isexported(NetworkCore, nm), frozen_public)
        @test sort(all_names) == sort(vcat(frozen, frozen_public))

        # The degree family is re-exported from Graphs.jl:
        # the same binding, so `using NetworkCore, Graphs` is not ambiguous.
        @test :degree in names(NetworkCore)
        @test :indegree in names(NetworkCore) && :outdegree in names(NetworkCore)
        @test NetworkCore.degree === Graphs.degree
        @test NetworkCore.indegree === Graphs.indegree && NetworkCore.outdegree === Graphs.outdegree
        # `edges(net)` yields Graphs.jl edges: the endpoint accessors travel
        # with the iterator (a term author's `for e in edges(net); src(e) ...`
        # was an UndefVarError after `using ERGM`).
        @test NetworkCore.src === Graphs.src && NetworkCore.dst === Graphs.dst
        @test :src in names(NetworkCore) && :dst in names(NetworkCore)
        let net = Network(3)
            add_edge!(net, 1, 3)
            @test [(src(e), dst(e)) for e in edges(net)] == [(1, 3)]
        end

        # The golden harness, bootstrap_cov and record_drop! stay EXPORTED here:
        # ten dependant test suites call `load_golden` unqualified after
        # `using NetworkCore`. A namespace leak would be in the
        # dependants' re-export loops, which curate this list.
        for nm in (:load_golden, :check_golden, :golden_report, :golden_tolerance,
                   :GoldenFixture, :bootstrap_cov, :record_drop!)
            @test nm in names(NetworkCore)
        end

        # Fresh process: a plain `using NetworkCore` answers the first questions a
        # user asks — construct a network and take a degree — with nothing else
        # loaded.
        code = "using NetworkCore; net = Network(5); add_edge!(net, 1, 2); " *
               "degree(net, 1) == 1 || error(\"degree\"); " *
               "src(first(edges(net))) == 1 || error(\"src\"); " *
               "dst(first(edges(net))) == 2 || error(\"dst\"); " *
               "indegree(net, 2) == 1 || error(\"indegree\"); " *
               "outdegree(net, 1) == 1 || error(\"outdegree\"); " *
               "Network(3) isa Network || error(\"Network\")"
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`
        @test success(cmd)
    end

    @testset "rem_edge!/add_edge! are allocation-free on attribute-free networks" begin
        # The sampler toggle: `rem_edge!` used to walk
        # every edge-attribute Dict per removal, halving MH throughput on any
        # attribute-bearing network. Now it never touches the attribute store on
        # an attribute-free network, walks it by value otherwise, and returns
        # early on an absent edge.
        for directed in (true, false)
            net = network(6; directed=directed)
            add_edges!(net, [(1, 3), (3, 2), (1, 4), (2, 4), (5, 6)])
            # warm up
            add_edge!(net, 1, 2); rem_edge!(net, 1, 2); add_edge!(net, 1, 2); rem_edge!(net, 1, 2)
            @test !has_edge(net, 1, 2)
            @test @allocated(rem_edge!(net, 1, 2)) == 0          # absent: early return
            add_edge!(net, 1, 2)
            @test @allocated(rem_edge!(net, 1, 2)) == 0          # present: real removal
            @test !has_edge(net, 1, 2)

            # `add_edge!` adds NOTHING on top of Graphs.jl's own `add_edge!` on
            # the bare SimpleDiGraph. (Base.insert! re-grows a Memory-backed
            # adjacency vector after deleteat! on alternate calls under Julia
            # 1.12, so a bare SimpleDiGraph itself is not always at 0 B; what
            # NetworkCore.jl guarantees is the wrapper costs zero bytes more.)
            g = copy(net.graph)
            toggle_net!(net) = (add_edge!(net, 1, 2); rem_edge!(net, 1, 2); nothing)
            function toggle_graph!(g)
                Graphs.add_edge!(g, Graphs.SimpleEdge{Int}(1, 2))
                directed || Graphs.add_edge!(g, Graphs.SimpleEdge{Int}(2, 1))
                Graphs.rem_edge!(g, Graphs.SimpleEdge{Int}(1, 2))
                directed || Graphs.rem_edge!(g, Graphs.SimpleEdge{Int}(2, 1))
                nothing
            end
            toggle_net!(net); toggle_graph!(g)
            net_bytes = @allocated(for _ in 1:64; toggle_net!(net); end)
            graph_bytes = @allocated(for _ in 1:64; toggle_graph!(g); end)
            @test net_bytes == graph_bytes
            @test @allocated(rem_edge!(net, 1, 2)) == 0

            # Correctness: attributes are still dropped when present ...
            set_edge_attribute!(net, :w, 1, 3, 1.5)
            set_edge_attribute!(net, :w, 3, 2, 2.5)
            set_edge_attribute!(net, :label, 3, 2, "x")
            @test rem_edge!(net, 3, 2)
            @test get_edge_attribute(net, :w, 3, 2) === nothing
            @test get_edge_attribute(net, :label, 3, 2) === nothing
            @test get_edge_attribute(net, :w, 1, 3) == 1.5          # neighbours untouched
            @test !has_edge(net, 3, 2) && (directed || !has_edge(net, 2, 3))
            # ... and the attribute-bearing removal is allocation-free too
            add_edge!(net, 3, 2); set_edge_attribute!(net, :w, 3, 2, 2.5); rem_edge!(net, 3, 2)
            add_edge!(net, 3, 2); set_edge_attribute!(net, :w, 3, 2, 2.5)
            @test @allocated(rem_edge!(net, 3, 2)) == 0
            @test get_edge_attribute(net, :w, 3, 2) === nothing
            # Removing an absent edge does not disturb attributes of other edges
            @test !rem_edge!(net, 5, 1)
            @test get_edge_attribute(net, :w, 1, 3) == 1.5
        end
    end

    @testset "Shared bootstrap: bootstrap_cov" begin
        # The ONE resampling loop of the ecosystem (src/bootstrap.jl). Every
        # `se=:bootstrap` in ERGM/ERGMCount/ERGMRank/ERGMMulti and REM's repeated
        # control sampling runs on it, so its contract is pinned here: replicates
        # are drawn once through the caller's rng, refits are collected, and the
        # empirical covariance comes back.
        θ̂ = [3.0]
        sim(rng, B) = [randn(rng, 200) .+ 3.0 for _ in 1:B]
        refit(x) = [sum(x) / length(x)]

        boot = bootstrap_cov(refit, sim, θ̂; n_boot=300, rng=MersenneTwister(1))

        # Recovers the known sampling SE of a mean: 1/sqrt(200)
        @test isapprox(boot.se[1], 1 / sqrt(200); rtol=0.15)
        @test size(boot.vcov) == (1, 1)
        @test size(boot.replicates) == (300, 1)
        @test boot.vcov[1, 1] ≈ boot.se[1]^2

        # Reproducible under a fixed rng — the whole point of routing the
        # randomness through `simulate`
        boot2 = bootstrap_cov(refit, sim, θ̂; n_boot=300, rng=MersenneTwister(1))
        @test boot.se == boot2.se
        @test boot.replicates == boot2.replicates
        @test bootstrap_cov(refit, sim, θ̂; n_boot=300,
                            rng=MersenneTwister(2)).se != boot.se

        # Threaded and serial refits agree exactly (the refits are deterministic
        # given their replicate, so the result is thread-count-independent)
        @test bootstrap_cov(refit, sim, θ̂; n_boot=50, rng=MersenneTwister(3),
                            threaded=false).replicates ==
              bootstrap_cov(refit, sim, θ̂; n_boot=50, rng=MersenneTwister(3),
                            threaded=true).replicates

        # A covariance needs at least two replicates
        @test_throws ArgumentError bootstrap_cov(refit, sim, θ̂; n_boot=1)
        # `simulate` must honour the requested number of replicates
        @test_throws ArgumentError bootstrap_cov(refit, (rng, B) -> sim(rng, B - 1),
                                                 θ̂; n_boot=10)

        # Multivariate: off-diagonal covariance is estimated, not just the
        # diagonal (two perfectly correlated coordinates)
        mv = bootstrap_cov(x -> [sum(x) / length(x), 2 * sum(x) / length(x)],
                           sim, [3.0, 6.0]; n_boot=200, rng=MersenneTwister(4))
        @test size(mv.vcov) == (2, 2)
        @test mv.vcov[1, 2] ≈ 2 * mv.vcov[1, 1]
        @test mv.se[2] ≈ 2 * mv.se[1]
    end

    @testset "bootstrap_cov: a refit error propagates unwrapped" begin
        sim(rng, B) = [randn(rng, 5) for _ in 1:B]
        bad(x) = throw(ArgumentError("refit refused"))
        for threaded in (true, false)
            err = try
                bootstrap_cov(bad, sim, [0.0]; n_boot=8, rng=Xoshiro(1), threaded=threaded)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError && err.msg == "refit refused"
            @test_throws ArgumentError bootstrap_cov(bad, sim, [0.0]; n_boot=8,
                                                     rng=Xoshiro(1), threaded=threaded)
        end
        # only some replicates fail: still the underlying exception
        k = Threads.Atomic{Int}(0)
        some(x) = (Threads.atomic_add!(k, 1); x[1] > 0 ? throw(DomainError(x[1], "positive")) : [x[1]])
        @test_throws DomainError bootstrap_cov(some, sim, [0.0]; n_boot=16, rng=Xoshiro(2))
        # results and rng order are untouched by the try/catch
        good(x) = [sum(x) / length(x)]
        a = bootstrap_cov(good, sim, [0.0]; n_boot=50, rng=Xoshiro(3), threaded=true)
        b = bootstrap_cov(good, sim, [0.0]; n_boot=50, rng=Xoshiro(3), threaded=false)
        @test a.replicates == b.replicates && a.se == b.se
        @test NetworkCore.unwrap_task_exception(ErrorException("x")) isa ErrorException
    end

    @testset "Shared se= validator: check_se" begin
        # ONE validator for the `se=` keyword, replacing
        # the seven `se in (...) || throw(...)` lines in ERGM/ERGMCount/ERGMMulti/
        # ERGMRank/TERGM/REM. Accepted symbol returned by identity; rejected
        # symbol's message names the context and the allowed tuple.
        @test check_se(:hessian, (:hessian, :bootstrap); context="mple") === :hessian
        @test check_se(:bootstrap, (:hessian, :bootstrap); context="mple") === :bootstrap
        # REM's two vocabularies (`fit_rem_from_data` has no bootstrap: the
        # controls are fixed in the DataFrame)
        @test check_se(:sandwich, (:hessian, :sandwich); context="fit_rem_from_data") === :sandwich
        @test check_se(:bootstrap, (:hessian, :sandwich, :bootstrap); context="fit_rem") === :bootstrap
        err = try
            check_se(:bootstrap, (:hessian, :sandwich); context="fit_rem_from_data")
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("fit_rem_from_data: se must be one of (:hessian, :sandwich) (got :bootstrap)",
                       err.msg)
        err2 = try
            check_se(:sandwich, (:hessian, :bootstrap); context="mple")
        catch e
            e
        end
        @test err2 isa ArgumentError && occursin("mple:", err2.msg) &&
              occursin("(:hessian, :bootstrap)", err2.msg) && occursin("(got :sandwich)", err2.msg)
        @test_throws ArgumentError check_se(:jackknife, (:hessian, :sandwich, :bootstrap); context="fit_rem")
        # A vector of allowed symbols prints as the same tuple
        err3 = try
            check_se(:x, [:hessian, :bootstrap]; context="c")
        catch e
            e
        end
        @test occursin("(:hessian, :bootstrap)", err3.msg)
        @test check_se(:hessian, [:hessian, :bootstrap]; context="c") === :hessian
    end

    @testset "Shared Newton optimizer" begin
        # `newton_fit` and `logistic_derivatives` are hosted here, `public` but
        # not exported until ERGM.jl aliases them (its
        # test suite calls `newton_fit` unqualified after `using ERGM, NetworkCore`).
        @test Base.ispublic(NetworkCore, :newton_fit) && !Base.isexported(NetworkCore, :newton_fit)
        @test Base.ispublic(NetworkCore, :logistic_derivatives) &&
              !Base.isexported(NetworkCore, :logistic_derivatives)
        newton_fit = NetworkCore.newton_fit
        logistic_derivatives = NetworkCore.logistic_derivatives

        @testset "quadratic: one exact step, exact vcov" begin
            A = [2.0 0.5; 0.5 1.0]
            target = [1.0, -0.5]
            quad(θ) = (-0.5 * dot(θ - target, A * (θ - target)), -A * (θ - target), -A)
            fit = newton_fit(quad, [10.0, -10.0])
            @test fit.converged
            @test fit.θ ≈ target atol = 1e-8
            @test fit.vcov ≈ inv(A) atol = 1e-8
            @test fit.se ≈ sqrt.([inv(A)[1, 1], inv(A)[2, 2]]) atol = 1e-8
            @test fit.loglik ≈ 0.0 atol = 1e-12
            @test fit.vcov isa Matrix{Float64} && fit.se isa Vector{Float64}
            @test fit.θ isa Vector{Float64}
            @test keys(fit) == (:θ, :se, :vcov, :loglik, :converged, :iterations, :stop)
            # θ0 is copied, not mutated
            θ0 = [10.0, -10.0]
            newton_fit(quad, θ0)
            @test θ0 == [10.0, -10.0]
            # Integer starting values are accepted
            @test newton_fit(quad, [10, -10]).θ ≈ target atol = 1e-8
        end

        @testset "Poisson closed form (step halving engages)" begin
            k = 7.0
            pois(θ) = (k * θ[1] - exp(θ[1]), [k - exp(θ[1])], hcat(-exp(θ[1])))
            pfit = newton_fit(pois, [8.0])
            @test pfit.converged
            @test pfit.θ[1] ≈ log(k) atol = 1e-6
            @test pfit.se[1] ≈ 1 / sqrt(k) atol = 1e-6
            @test pfit.loglik ≈ k * log(k) - k atol = 1e-8
            @test pfit.iterations >= 1
            # The docstring example
            dfit = newton_fit(pois, [0.0])
            @test dfit.θ[1] ≈ log(k) && dfit.se[1] ≈ 1 / sqrt(k) && dfit.converged
            # From θ = -10 the Hessian is ~e^-10, so the first Newton step is
            # ~1.5e5 too long and the first improving step is 2^-14 of it. Ten
            # halvings stop short of it; the default does not.
            far = newton_fit(pois, [-10.0])
            @test far.converged
            @test far.θ[1] ≈ log(k) atol = 1e-6
            short = newton_fit(pois, [-10.0]; max_halvings=10)
            @test !short.converged && short.θ[1] == -10.0
            # maxiter=1 from far away: honestly unconverged
            one = newton_fit(pois, [8.0]; maxiter=1)
            @test one.converged == false
            @test one.iterations == 1
            @test isfinite(one.loglik) && isfinite(one.θ[1])
        end

        @testset "finiteness guard: −Inf on the full step" begin
            # ll(θ) = k log θ − θ on θ > 0, −Inf outside. From θ0 = 20 the full
            # Newton step overshoots to θ < 0; the halving must reject it and
            # still converge to θ̂ = k.
            k = 7.0
            barrier(θ) = θ[1] <= 0 ? (-Inf, [NaN], hcat(NaN)) :
                         (k * log(θ[1]) - θ[1], [k / θ[1] - 1], hcat(-k / θ[1]^2))
            bfit = newton_fit(barrier, [20.0])
            @test bfit.converged
            @test bfit.θ[1] ≈ k atol = 1e-6
            @test bfit.se[1] ≈ sqrt(k) atol = 1e-6
            @test all(isfinite, bfit.θ) && isfinite(bfit.loglik)
            # NaN on the full step is treated the same way
            nanwall(θ) = θ[1] <= 0 ? (NaN, [NaN], hcat(NaN)) : barrier(θ)
            @test newton_fit(nanwall, [20.0]).θ[1] ≈ k atol = 1e-6
            # Every halving non-finite: stop at the last finite iterate, loudly
            # unconverged, never a NaN θ
            # (a Newton step of 1e15, so that even the smallest of the 30 halved
            # steps, 2^-29 of it, lands beyond the wall)
            wall(θ) = θ[1] > 1.0 ? (-Inf, [NaN], hcat(NaN)) : (θ[1], [1.0], hcat(-1e-15))
            wfit = newton_fit(wall, [0.5])
            @test !wfit.converged
            @test wfit.θ == [0.5] && wfit.loglik == 0.5
            # A non-finite objective at the start is a caller error
            @test_throws ArgumentError newton_fit(barrier, [-1.0])
        end

        @testset "indefinite Hessian → NaN standard errors, warned" begin
            saddle(θ) = (-(θ[1]^2) / 2 + (θ[2]^2) / 2, [-θ[1], θ[2]], [-1.0 0.0; 0.0 1.0])
            sfit = @test_logs (:warn, r"not negative definite") newton_fit(saddle, zeros(2))
            @test all(isnan, sfit.se) && all(isnan, sfit.vcov)
            @test sfit.θ == zeros(2)
            # Positive semi-definite (a flat direction) is not enough either:
            # at the maximum the singular Hessian gives no Newton direction
            # (unconverged) and no covariance
            flat(θ) = (-(θ[1]^2) / 2, [-θ[1], 0.0], [-1.0 0.0; 0.0 0.0])
            ffit = @test_logs (:warn, r"not negative definite") newton_fit(flat, [0.0, 0.0])
            @test all(isnan, ffit.se) && all(isnan, ffit.vcov)
            @test ffit.θ == [0.0, 0.0] && !ffit.converged
        end

        @testset "singular Hessian: no Newton direction, unconverged" begin
            sing(θ) = (θ[1], [1.0], hcat(0.0))
            sfit = @test_logs (:warn, r"not negative definite") newton_fit(sing, [0.0])
            @test !sfit.converged
            @test sfit.iterations == 1
        end

        @testset "information rank is invariant to coefficient units" begin
            # A constant timing covariate and the fitted baseline describe
            # the same direction. Cholesky can accept these rank-one matrices
            # on some BLAS platforms; no finite marginal SE is identified.
            for multiplier in (0.1, 4.0, 17.3), scales in ([1.0, 1.0], [1e-100, 1e100])
                H = multiplier .* (scales * scales')
                V, se = @test_logs (:warn, r"numerically identifiable") NetworkCore._observed_information_cov(-H, 2)
                @test all(isnan, V) && all(isnan, se)
            end
            # Mere unit changes must not reject a well identified model.
            C = [1.0 0.25; 0.25 1.0]
            for scales in ([1.0, 1.0], [1e-100, 1e100], [1e100, 1e-100])
                H = C .* (scales * scales')
                V, se = NetworkCore._observed_information_cov(-H, 2)
                @test all(isfinite, V) && all(isfinite, se)
                @test V .* (scales * scales') ≈ inv(C) rtol=1e-14
                @test se .* scales ≈ sqrt.(diag(inv(C))) rtol=1e-14
            end
            # Resolve a difficult but identifiable direction above the
            # numerical-rank boundary; reject one at rounding precision.
            rho = 1.0 - 1e-10
            near = [1.0 rho; rho 1.0]
            @test all(isfinite, first(NetworkCore._observed_information_cov(-near, 2)))
            unresolved = [1.0 prevfloat(1.0); prevfloat(1.0) 1.0]
            V, se = @test_logs (:warn, r"numerically identifiable") NetworkCore._observed_information_cov(-unresolved, 2)
            @test all(isnan, V) && all(isnan, se)
            saddle(θ) = (-θ[1]^2 / 2 + θ[2]^2 / 2, [-θ[1], θ[2]], [-1.0 0.0; 0.0 1.0])
            result = @test_logs (:warn, r"not negative definite") newton_fit(saddle, zeros(2))
            @test !result.converged
        end

        @testset "information rank tolerance reflects Hessian evaluation precision" begin
            # An assembled Hessian can have more rounding noise than the
            # p-dimensional eigensolve itself. Let its owner declare that
            # precision without changing the optimizer's default policy.
            H = [1.0 1.0-1e-14; 1.0-1e-14 1.0]
            quadratic(θ) = (-dot(θ, H * θ) / 2, -H * θ, -H)
            default = newton_fit(quadratic, zeros(2))
            @test default.converged && all(isfinite, default.vcov)
            guarded = @test_logs (:warn, r"numerically identifiable") newton_fit(
                quadratic, zeros(2); information_rtol=1e-12)
            @test !guarded.converged
            @test all(isnan, guarded.vcov) && all(isnan, guarded.se)
            # Equilibration still protects well-identified parameters with
            # very different units; absolute eigenvalues must not decide.
            scaled(θ) = (-θ[1]^2 / 2 - 1e-20 * θ[2]^2 / 2,
                         [-θ[1], -1e-20 * θ[2]], [-1.0 0.0; 0.0 -1e-20])
            healthy = newton_fit(scaled, zeros(2); information_rtol=1e-12)
            @test healthy.converged && healthy.se ≈ [1.0, 1e10]
            for bad in (-1.0, 1.0, Inf, NaN)
                @test_throws ArgumentError newton_fit(quadratic, zeros(2); information_rtol=bad)
            end
        end

        @testset "scale-free verdict: Newton decrement and rounding noise" begin
            # On a large design (TERGM's CMPLE, ERGMRank's swap MPLE, REM) the
            # gradient's own units put `sqrt(tol)` below the objective's
            # rounding floor. One Newton step short of the maximum of
            # ℓ = −½·10⁶·(θ−1)² − 1000 the gradient is 1e-3 (> sqrt(tol)) while
            # the step could gain only 5e-13 (< tol), and the evaluated ℓ is
            # flat to a few ulps (eps(1000) ≈ 1.1e-13). The rounding error is
            # made deterministic here — a penalty of `k` ulps inside |θ−1| <
            # 1e-9, i.e. on the full step and on every halving — so that every
            # halving "decreases" the objective exactly as BLAS rounding did
            # on the TERGM golden panel. The gradient and Hessian are exact.
            h = 1e6
            function noisy(k)
                return θ -> begin
                    d = θ[1] - 1
                    ll = -0.5 * h * d^2 - 1000.0 - k * eps(1000.0) * (abs(d) < 1e-9)
                    (ll, [-h * d], hcat(-h))
                end
            end
            # (a) noise beyond the 8-ulp acceptance band: no halving is taken,
            # but the decrement 5e-13 < tol says the full step is the polish —
            # taken, converged, at the maximum (R's glm reports converged here)
            big = newton_fit(noisy(20), [1 + 1e-9])
            @test big.converged
            @test abs(big.θ[1] - 1) < 1e-15
            @test big.se[1] ≈ 1 / sqrt(h)
            # (b) noise within the band and the gradient shrinking with it: the
            # full step is accepted as a non-decrease (REM's 1-ulp case) and
            # the next iteration converges; θ moved to the optimum
            small = newton_fit(noisy(6), [1 + 1e-9])
            @test small.converged
            @test abs(small.θ[1] - 1) < 1e-15
            # Same problem without any noise: one exact step, unchanged verdict
            clean = newton_fit(noisy(0), [1 + 1e-9])
            @test clean.converged && abs(clean.θ[1] - 1) < 1e-15
            # The decrement is scale-free: the SAME stall on a unit-Hessian
            # problem is a gradient of 1e-9, converged under both criteria
            unit = newton_fit(θ -> (-0.5 * (θ[1] - 1)^2 - 1000.0, [-(θ[1] - 1)], hcat(-1.0)),
                              [1 + 1e-9])
            @test unit.converged
            # (c) Negative control: a genuinely unconverged stall (the step
            # predicts a gain of 50, every halving loses, the gradient never
            # shrinks) stays put and stays unconverged — the decrement rule
            # never rescues a large predicted gain
            liar(θ) = (-(θ[1]^2), [10.0], hcat(-1.0))
            stall = newton_fit(liar, [0.0])
            @test !stall.converged
            @test stall.θ == [0.0] && stall.loglik == 0.0
            # (d) The gain the step could make is what decides, never the
            # noise: under a `tol` below that gain (5e-13) the same stall is
            # no polish — no halving is taken, the gradient (1e-3) is not
            # below sqrt(tol), and the fit is honestly unconverged, θ unchanged
            tight = newton_fit(noisy(20), [1 + 1e-9]; tol=1e-14)
            @test !tight.converged && tight.θ == [1 + 1e-9]
        end

        @testset "a bug in the closure surfaces; only a failed factorisation is 'no direction'" begin
            # The Newton-step `try` used to be a bare `catch` that turned ANY
            # error into `converged=false, iterations=1` with NaN standard
            # errors — a swallowed programming error dressed up as
            # non-convergence. Now only Singular/LAPACK/ZeroPivot/PosDef
            # exceptions mean "no Newton direction".
            bad_grad(θ) = (-sum(abs2, θ), [1.0, 2.0, 3.0], [-2.0 0.0; 0.0 -2.0])   # p = 2, ∇ has 3
            @test_throws DimensionMismatch newton_fit(bad_grad, [1.0, 1.0])
            bad_hess(θ) = (-sum(abs2, θ), -2 .* θ, [-2.0 0.0])                    # 1×2 Hessian
            @test_throws DimensionMismatch newton_fit(bad_hess, [1.0, 1.0])
            # A non-matrix Hessian fails inside the try (`-hess`): rethrown, not caught
            not_a_matrix(θ) = (-θ[1]^2, [-2θ[1]], "hess")
            @test_throws MethodError newton_fit(not_a_matrix, [1.0])
            # ... and an error raised by `\` itself that is not a factorisation
            # failure propagates too
            exploding(θ) = (-θ[1]^2, [-2θ[1]], ExplodingHessian())   # declared at top level
            @test_throws ErrorException newton_fit(exploding, [1.0])
        end

        @testset "logistic_derivatives: closed-form edges-only MPLE" begin
            # Florentine marriage: 16 actors, 120 undirected dyads, 20 ties.
            # θ̂ = logit(20/120) = log(20/100), se = 1/sqrt(120 p (1−p)). This is
            # what ERGM.jl's phase re-runs against flomarriage_ergm.toml at 1e-6
            # before dropping Optim.
            d = logistic_derivatives(ones(1, 1), [120.0], [20.0])
            fit = newton_fit(d, [0.0])
            @test fit.converged
            @test fit.θ[1] ≈ log(20 / 100) atol = 1e-10
            @test fit.se[1] ≈ sqrt(1 / (120 * (1 / 6) * (5 / 6))) atol = 1e-10
            @test fit.loglik ≈ 20 * log(1 / 6) + 100 * log(5 / 6) atol = 1e-10
            # The bundled dataset agrees with the literal
            flo = load_dataset(:florentine_marriage)
            @test nv(flo) == 16 && ne(flo) == 20
            # The same fit from 120 binary rows and from 120 unit weights
            X = ones(120, 1)
            y = [i <= 20 for i in 1:120]
            @test newton_fit(logistic_derivatives(X, y), [0.0]).θ[1] ≈ log(20 / 100) atol = 1e-10
            @test newton_fit(logistic_derivatives(X, y; weights=ones(120)), [0.0]).se[1] ≈
                  fit.se[1] atol = 1e-10
            # Integer design / counts are accepted
            @test newton_fit(logistic_derivatives(ones(Int, 1, 1), [120], [20]), [0.0]).θ[1] ≈
                  log(20 / 100) atol = 1e-10
        end

        @testset "logistic_derivatives: weighted ≡ duplicated rows, binomial ≡ both" begin
            rng = Xoshiro(11)
            n, p = 40, 3
            Xu = randn(rng, n, p)
            yu = rand(rng, Bool, n)
            w = rand(rng, 1:4, n)
            off = 0.3 .* randn(rng, n)
            # Duplicate every row w_r times
            idx = vcat([fill(r, w[r]) for r in 1:n]...)
            Xd, yd, offd = Xu[idx, :], yu[idx], off[idx]
            β = [0.4, -0.2, 0.1]
            ll_d, g_d, h_d = logistic_derivatives(Xd, yd; offset=offd)(β)
            ll_w, g_w, h_w = logistic_derivatives(Xu, yu; weights=Float64.(w), offset=off)(β)
            @test ll_w ≈ ll_d atol = 1e-12 * abs(ll_d)
            @test g_w ≈ g_d atol = 1e-12 * norm(g_d)
            @test h_w ≈ h_d atol = 1e-12 * norm(h_d)
            # Binomial rows: n_tot = w, n_one = w·y
            ll_b, g_b, h_b = logistic_derivatives(Xu, Float64.(w), Float64.(w .* yu); offset=off)(β)
            @test ll_b ≈ ll_d atol = 1e-12 * abs(ll_d)
            @test g_b ≈ g_d atol = 1e-12 * norm(g_d)
            @test h_b ≈ h_d atol = 1e-12 * norm(h_d)
            # Unit weights ≡ no weights, exactly
            @test logistic_derivatives(Xu, yu)(β) == logistic_derivatives(Xu, yu; weights=ones(n))(β)
            # A zero-weight row contributes nothing
            w0 = ones(n); w0[3] = 0.0
            keep = setdiff(1:n, 3)
            @test logistic_derivatives(Xu, yu; weights=w0)(β)[1] ≈
                  logistic_derivatives(Xu[keep, :], yu[keep])(β)[1]
            # Fractional successes (a binomial row need not be integer)
            @test isfinite(logistic_derivatives(Xu, fill(2.5, n), fill(1.25, n))(β)[1])
            # The two fits agree on θ̂ and se
            fd = newton_fit(logistic_derivatives(Xd, yd), zeros(p))
            fw = newton_fit(logistic_derivatives(Xu, yu; weights=Float64.(w)), zeros(p))
            fb = newton_fit(logistic_derivatives(Xu, Float64.(w), Float64.(w .* yu)), zeros(p))
            @test fd.converged && fw.converged && fb.converged
            @test fw.θ ≈ fd.θ atol = 1e-8
            @test fb.θ ≈ fd.θ atol = 1e-8
            @test fw.se ≈ fd.se atol = 1e-8
            @test fb.se ≈ fd.se atol = 1e-8
        end

        @testset "logistic_derivatives: analytic gradient vs central differences" begin
            rng = Xoshiro(7)
            n, p = 60, 4
            X = randn(rng, n, p)
            y = rand(rng, Bool, n)
            off = randn(rng, n)
            w = rand(rng, n) .+ 0.5
            d = logistic_derivatives(X, y; offset=off, weights=w)
            β = [0.3, -0.7, 0.2, 0.05]
            ll, grad, hess = d(β)
            h = 1e-5
            for k in 1:p
                e = zeros(p); e[k] = h
                fd = (d(β .+ e)[1] - d(β .- e)[1]) / (2h)
                @test grad[k] ≈ fd atol = 1e-6
                # Hessian column from the gradient
                fdh = (d(β .+ e)[2] .- d(β .- e)[2]) ./ (2h)
                @test hess[:, k] ≈ fdh atol = 1e-6
            end
            @test hess ≈ hess'                                # symmetric
            @test all(eigvals(Symmetric(hess)) .< 0)          # negative definite
            # Large |η| is numerically stable (no NaN/Inf from exp overflow)
            big = logistic_derivatives(1000.0 .* X, y)([1.0, 1.0, 1.0, 1.0])
            @test isfinite(big[1]) && all(isfinite, big[2]) && all(isfinite, big[3])
        end

        @testset "logistic_derivatives: validation" begin
            X = ones(4, 2)
            @test_throws ArgumentError logistic_derivatives(X, [true, false])
            @test_throws ArgumentError logistic_derivatives(X, trues(4); offset=[1.0])
            @test_throws ArgumentError logistic_derivatives(X, trues(4); weights=[1.0])
            @test_throws ArgumentError logistic_derivatives(X, trues(4); weights=[1.0, -1.0, 1.0, 1.0])
            @test_throws ArgumentError logistic_derivatives(X, trues(4); weights=[1.0, NaN, 1.0, 1.0])
            @test_throws ArgumentError logistic_derivatives(X, [1.0, 1.0], [1.0, 1.0, 1.0, 1.0])
            @test_throws ArgumentError logistic_derivatives(X, ones(4), [1.0, 1.0])
            @test_throws ArgumentError logistic_derivatives(X, ones(4), [1.0, 2.0, 0.0, 0.0])   # n_one > n_tot
            @test_throws ArgumentError logistic_derivatives(X, ones(4), [1.0, -0.5, 0.0, 0.0])
            @test_throws ArgumentError logistic_derivatives(X, [-1.0, 1.0, 1.0, 1.0], zeros(4))
            @test_throws ArgumentError logistic_derivatives(X, [Inf, 1.0, 1.0, 1.0], zeros(4))
            @test_throws ArgumentError logistic_derivatives(X, ones(4), zeros(4); offset=[1.0])
            # Empty design: zero objective, zero derivatives
            e = logistic_derivatives(zeros(0, 2), Bool[])(zeros(2))
            @test e == (0.0, zeros(2), zeros(2, 2))
        end

        @testset "logistic_derivatives allocates one gradient and one Hessian per call" begin
            # The workspaces (η, residual, weighted design) are
            # built once in the closure; an evaluation allocates exactly the
            # length-p gradient and p×p Hessian it returns — not O(n·p²). Pinned
            # by equality with the bytes of those two arrays.
            refalloc(p) = (Vector{Float64}(undef, p), Matrix{Float64}(undef, p, p))
            for p in (1, 3, 6)
                rng = Xoshiro(p)
                nrows = 800
                X = randn(rng, nrows, p)
                y = rand(rng, Bool, nrows)
                β = fill(0.1, p)
                d = logistic_derivatives(X, y)
                dw = logistic_derivatives(X, y; offset=randn(rng, nrows), weights=rand(rng, nrows) .+ 1)
                db = logistic_derivatives(X, fill(3.0, nrows), [y[i] ? 2.0 : 1.0 for i in 1:nrows])
                d(β); dw(β); db(β); refalloc(p)      # warm up
                ref = @allocated refalloc(p)
                @test @allocated(d(β)) == ref
                @test @allocated(dw(β)) == ref
                @test @allocated(db(β)) == ref
                # and it does not grow with the number of rows
                Xb = randn(rng, 20 * nrows, p)
                dbig = logistic_derivatives(Xb, rand(rng, Bool, 20 * nrows))
                dbig(β)
                @test @allocated(dbig(β)) == ref
            end
        end
    end

    # ------------------------------------------------------------------
    # Integer vertex types, Graphs.jl generics, the DataFrames hint, the
    # data directory, and Graphs.jl algorithms on the backing store.
    # ------------------------------------------------------------------

    @testset "Two-mode structural edits keep the partition (fuzz vs brute force)" begin
        # Brute-force model: a dense adjacency matrix, the mode-1 size, one
        # attribute per vertex (its ORIGINAL id) and the mask as a Set of
        # original-id pairs. Every edit is applied to both and compared.
        function check(net, A, k, ids, mask, directed)
            @test nv(net) == size(A, 1) == length(ids)
            @test net.bipartite == k
            M = as_matrix(net; expand_bipartite=true)
            @test M == A
            @test ne(net) == (directed ? count(A) : count(A) ÷ 2)
            @test vertex_attribute_vector(net, :id, Int) == ids
            # edge attributes are the original endpoint pair
            for e in edges(net)
                @test get_edge_attribute(net, :orig, src(e), dst(e)) ==
                      (directed ? (ids[src(e)], ids[dst(e)]) : minmax(ids[src(e)], ids[dst(e)]))
            end
            got = Set((ids[i], ids[j]) for (i, j) in missing_dyads(net))
            @test got == Set(directed ? d : minmax(d...) for d in mask if d[1] in ids && d[2] in ids)
            # the partition is still a partition: every edge is cross-mode,
            # and the incidence form agrees with the square one
            @test all(e -> (src(e) <= k) ⊻ (dst(e) <= k), edges(net))
            directed || @test as_matrix(net) == M[1:k, k+1:end]   # directed incidence drops (and warns about) mode-2 → mode-1 arcs
        end
        rng = Xoshiro(20261002)
        for trial in 1:60
            directed = isodd(trial)
            n1, n2 = rand(rng, 1:5), rand(rng, 1:6)
            n = n1 + n2
            net = network(n; bipartite=n1, directed=directed)
            A = falses(n, n); k = n1; ids = collect(1:n); mask = Set{Tuple{Int,Int}}()
            set_vertex_attribute!(net, :id, collect(1:n))
            for _ in 1:3n
                i, j = rand(rng, 1:n1), rand(rng, n1+1:n)
                rand(rng, Bool) && directed && ((i, j) = (j, i))
                if add_edge!(net, i, j)
                    set_edge_attribute!(net, :orig, i, j, directed ? (i, j) : minmax(i, j))
                    A[i, j] = true; directed || (A[j, i] = true)
                end
            end
            for _ in 1:2
                i, j = rand(rng, 1:n1), rand(rng, n1+1:n)
                set_missing_dyad!(net, i, j); push!(mask, directed ? (i, j) : minmax(i, j))
            end
            check(net, A, k, ids, mask, directed)

            # a within-mode permutation; a cross-mode one is refused
            perm = vcat(shuffle(rng, 1:k), shuffle(rng, k+1:n))
            p = permute_vertices(net, perm)
            check(p, A[perm, perm], k, ids[perm], mask, directed)
            if k < n
                bad = copy(perm); bad[1], bad[end] = bad[end], bad[1]
                @test_throws ArgumentError permute_vertices(net, bad)
            end

            # a sorted induced subgraph stays two-mode with the selected mode-1 count
            sel = sort(shuffle(rng, 1:n)[1:rand(rng, 1:n)])
            sub = get_induced_subgraph(net, sel)
            check(sub, A[sel, sel], count(<=(k), sel), ids[sel], mask, directed)

            # remove vertices one at a time (stable, R's delete.vertices)
            while nv(net) > 1 && rand(rng) < 0.7
                v = rand(rng, 1:nv(net))
                @test rem_vertex!(net, v)
                keep = setdiff(1:size(A, 1), v)
                A = A[keep, keep]; k -= (v <= k); ids = ids[keep]
                check(net, A, k, ids, mask, directed)
            end
            # and several at once
            if nv(net) > 2
                vs = unique(rand(rng, 1:nv(net), 2))
                vmap = Graphs.rem_vertices!(net, vs)
                keep = setdiff(1:size(A, 1), vs)
                @test vmap == keep
                A = A[keep, keep]; k -= count(<=(k), vs); ids = ids[keep]
                check(net, A, k, ids, mask, directed)
            end
        end

        # The original failing cases, verbatim
        b = network(5; bipartite=2, directed=false)
        add_edges!(b, [(1, 3), (1, 4), (2, 5)])
        @test_throws ArgumentError permute_vertices(b, [3, 1, 2, 4, 5])   # used to drop 2 of 3 edges
        b2 = network(5; bipartite=2, directed=false)
        add_edges!(b2, [(1, 3), (2, 5), (1, 5)])
        rem_vertex!(b2, 1)
        @test b2.bipartite == 1                                    # R: delete.vertices → bipartite 1
        @test [(src(e), dst(e)) for e in edges(b2)] == [(1, 4)]   # old (2, 5)
        @test as_matrix(b2) == [0.0 0.0 1.0]                       # used to throw BoundsError
        @test !add_edge!(b2, 1, 1) && !has_edge(b2, 1, 2)
        sub = get_induced_subgraph(b, [1, 3, 4])
        @test sub.bipartite == 1 && ne(sub) == 2                   # R: get.inducedSubgraph keeps bipartite = 1
        @test_throws ArgumentError get_induced_subgraph(b, [3, 1, 4])  # interleaved modes
        @test_throws ArgumentError Graphs.rem_vertices!(b, [1, 9])
        @test nv(b) == 5                                           # refused before touching anything
        # one-mode rem_vertex! keeps its Graphs.jl swap-and-pop convention
        o = network(3); set_vertex_attribute!(o, :name, ["a", "b", "c"])
        rem_vertex!(o, 1)
        @test vertex_attribute_vector(o, :name, String) == ["c", "b"]
        # ... while rem_vertices! is stable on one-mode networks too
        o2 = network(4); set_vertex_attribute!(o2, :name, ["a", "b", "c", "d"])
        add_edges!(o2, [(1, 4), (2, 3)])
        @test Graphs.rem_vertices!(o2, [1]) == [2, 3, 4]
        @test vertex_attribute_vector(o2, :name, String) == ["b", "c", "d"]
        @test has_edge(o2, 1, 2) && ne(o2) == 1
    end

    @testset "Non-Int64 vertex types: Network{Int32}" begin
        for directed in (true, false), loops in (false, true)
            net = Network{Int32}(; n=Int32(6), directed=directed, loops=loops)
            @test net isa Network{Int32, directed}
            add_edges!(net, [(Int32(1), Int32(2)), (Int32(2), Int32(3)), (Int32(4), Int32(5))])
            loops && add_edge!(net, Int32(6), Int32(6))
            set_vertex_attribute!(net, :g, Dict(Int32(v) => v % 2 for v in 1:6))
            set_edge_attribute!(net, :w, Int32(1), Int32(2), 2.5)
            set_missing_dyad!(net, Int32(1), Int32(3))
            c = copy(net)                                  # used to throw a TypeError
            @test c isa Network{Int32, directed}
            @test ne(c) == ne(net) == 3 + loops
            @test as_matrix(c) == as_matrix(net)
            @test get_edge_attribute(c, :w, 1, 2) == 2.5
            @test is_missing_dyad(c, 1, 3)
            @test add_edge!(c, Int32(5), Int32(6))
            @test rem_edge!(c, Int32(1), Int32(2))
            @test ne(c) == ne(net)                          # one added, one removed
            @test ne(net) == 3 + loops                      # the original is untouched
            @test size(as_matrix(c)) == (6, 6)
            @test rem_vertex!(c, 2) && nv(c) == 5
            @test eltype(collect(vertices(c))) == Int32
            p = permute_vertices(net, [6, 5, 4, 3, 2, 1])
            @test p isa Network{Int32, directed} && has_edge(p, 5, 4)
            s = get_induced_subgraph(net, [1, 2, 3])
            @test s isa Network{Int32, directed} && ne(s) == 2
            @test zero(net) isa Network{Int32, directed}
            @test complement(net) isa Network{Int32, directed}
        end
        # The wrapper carries D, so its `network` field is concrete
        for directed in (true, false)
            bw = BipartiteNetwork(2, 3; directed=directed)
            @test bw isa BipartiteNetwork{Int, directed} && bw isa BipartiteNetwork{Int}
            @test isconcretetype(fieldtype(typeof(bw), :network))
            @test bw.network isa Network{Int, directed}
            @test typeof(copy(bw)) === typeof(bw)
            @test is_directed(bw) == is_directed(typeof(bw)) == directed
            @test (@inferred ne(bw)) == 0
        end
        @test BipartiteNetwork{Int, false}(2, 3; loops=true).network.loops
        b = BipartiteNetwork{Int32}(Int32(2), Int32(3))
        add_edge!(b, Int32(1), Int32(4))
        @test ne(copy(b)) == 1 && copy(b).network isa Network{Int32}
        @test network(Int32(4)) isa Network{Int, true}       # any Integer n
        @test nv(Network(UInt8(3))) == 3
    end

    @testset "Graphs.jl induced_subgraph and complement on a Network" begin
        rng = Xoshiro(16)
        for directed in (true, false), loops in (false, true)
            net = network(12; directed=directed, loops=loops)
            for _ in 1:30
                i, j = rand(rng, 1:12), rand(rng, 1:12)
                (i != j || loops) && add_edge!(net, i, j)
            end
            set_vertex_attribute!(net, :id, collect(1:12))
            g = directed ? SimpleDiGraph(12) : SimpleGraph(12)
            for e in edges(net); add_edge!(g, src(e), dst(e)); end

            vl = [3, 1, 7, 12, 5]
            sub, vmap = induced_subgraph(net, vl)                 # used to be a MethodError
            gsub, gmap = induced_subgraph(g, vl)
            @test vmap == gmap == vl
            @test Matrix(adjacency_matrix(sub)) == Matrix(adjacency_matrix(gsub))
            @test vertex_attribute_vector(sub, :id, Int) == vl     # attributes carried
            @test sub.loops == loops
            mask = falses(12); mask[vl] .= true
            @test induced_subgraph(net, mask)[2] == sort(vl)
            @test net[vl] isa Network                             # Graphs' getindex sugar
            @test_throws ArgumentError induced_subgraph(net, [1, 1])

            c = complement(net)                                   # used to be a MethodError
            @test c isa typeof(net) && c.loops == loops
            for i in 1:12, j in 1:12
                (i == j && !loops) && continue
                @test has_edge(c, i, j) == !has_edge(net, i, j)
            end
            loops || @test Matrix(adjacency_matrix(c)) == Matrix(adjacency_matrix(complement(g)))
            @test vertex_attribute_vector(c, :id, Int) == collect(1:12)
            @test complement(c) |> x -> Matrix(adjacency_matrix(x)) == Matrix(adjacency_matrix(net))
        end
        # two-mode: only cross-mode dyads, partition kept; the mask is kept
        b = network(5; bipartite=2, directed=false)
        add_edges!(b, [(1, 3), (1, 4), (2, 5)])
        set_missing_dyad!(b, 2, 3)
        cb = complement(b)
        @test cb.bipartite == 2
        @test sort([(src(e), dst(e)) for e in edges(cb)]) == [(1, 5), (2, 3), (2, 4)]   # R: !net
        @test is_missing_dyad(cb, 2, 3)
        # an edge attribute does not survive (the complement's edges are non-edges)
        set_edge_attribute!(b, :w, 1, 3, 1.0)
        @test isempty(list_edge_attributes(complement(b)))
        # self-loops are part of the complement exactly when loops=true (R: !net)
        l = network(3; loops=true)
        add_edges!(l, [(1, 1), (1, 2)])
        @test sort([(src(e), dst(e)) for e in edges(complement(l))]) ==
              [(1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
    end

    @testset "ne is O(1) with loops: the cached self-loop count stays exact" begin
        brute(net) = net.directed ? Graphs.ne(net.graph) :
            (Graphs.ne(net.graph) + count(v -> Graphs.has_edge(net.graph, v, v), 1:nv(net))) ÷ 2
        rng = Xoshiro(19)
        for directed in (false, true), bip in (nothing, 3)
            net = network(8; directed=directed, loops=bip === nothing, bipartite=bip)
            for step in 1:400
                r = rand(rng)
                i, j = rand(rng, 1:nv(net)), rand(rng, 1:nv(net))
                if r < 0.45
                    add_edge!(net, i, j)
                elseif r < 0.85
                    rem_edge!(net, i, j)
                elseif r < 0.9 && nv(net) > 2
                    rem_vertex!(net, i)
                elseif r < 0.93
                    add_vertex!(net)
                elseif r < 0.96
                    net = copy(net)
                elseif nv(net) > 3
                    net = permute_vertices(net, bip === nothing ? shuffle(rng, 1:nv(net)) :
                        vcat(shuffle(rng, 1:net.bipartite), shuffle(rng, net.bipartite+1:nv(net))))
                end
                @test ne(net) == brute(net)
                @test ne(net) == length(collect(edges(net)))
            end
        end
        # O(1): no allocation, and no dependence on n (it used to count the
        # loops over every vertex, 5.6 μs at n = 2000)
        big = network(20_000; directed=false, loops=true)
        add_edges!(big, [(1, 1), (1, 2)])
        ne(big)
        @test @allocated(ne(big)) == 0
        @test ne(big) == 2
        # assigning the backing graph re-derives the count
        g = SimpleDiGraph(3); add_edge!(g, 2, 2); add_edge!(g, 1, 3); add_edge!(g, 3, 1)
        u = network(3; directed=false, loops=true)
        u.graph = g
        @test ne(u) == 2
    end

    @testset "Separation verdict: every family against the brute-force oracle" begin
        # The verdict (a linear programme plus an exact certificate) against
        # the exhaustive extreme-ray enumeration above, on random small
        # integer designs of each likelihood family; about half are separated.
        # For each separated verdict: it is certified, its direction is a
        # direction of recession supported on `terms`, and `terms` is
        # inclusion-minimal (no proper subset carries one).
        rng = Xoshiro(20261006)
        tally = Dict(f => [0, 0, 0] for f in (:logit, :binomial, :poisson, :clogit))
        for _ in 1:1200
            q = rand(rng, 1:3); n = rand(rng, 2:8)
            X = rand(rng, -2:2, n, q)
            rand(rng) < 0.5 && (X[:, 1] .= 1)
            fam = rand(rng, (:logit, :binomial, :poisson, :clogit))
            if fam === :logit
                y = rand(rng, Bool, n)
                v = NetworkCore.logistic_separation(Float64.(X), y)
                M = oracle_logit_rows(X, y)
            elseif fam === :binomial
                nt = rand(rng, 0:3, n); n1 = [rand(rng, 0:t) for t in nt]
                v = NetworkCore.logistic_separation(Float64.(X), Float64.(nt), Float64.(n1))
                M = oracle_binomial_rows(X, nt, n1)
            elseif fam === :poisson
                y = rand(rng, 0:2, n) .* rand(rng, Bool, n)
                v = NetworkCore.poisson_separation(Float64.(X), y)
                M = oracle_poisson_rows(X, y)
            else
                ns = rand(rng, 1:3)
                st = sort(vcat(collect(1:ns), rand(rng, 1:ns, max(n - ns, 0))))
                ch = falses(length(st))
                for s in 1:ns
                    idx = findall(==(s), st)
                    ch[rand(rng, idx)] = true
                    rand(rng) < 0.3 && (ch[rand(rng, idx)] = true)   # a multi-case stratum
                end
                X = rand(rng, -2:2, length(st), q)
                v = NetworkCore.clogit_separation(Float64.(X), ch, st)
                M = oracle_clogit_rows(X, ch, st)
            end
            truth = oracle_separated(M)
            t = tally[fam]
            t[1] += 1; t[2] += truth
            ok = v.separated == truth
            if ok && truth
                d = v.direction
                margins = M * d
                ok &= v.certified
                ok &= all(>=(-1e-9), margins) && maximum(margins) > 1e-9
                ok &= all(j -> d[j] == 0, setdiff(1:q, v.terms))
                ok &= length(v.terms) == 1 ||
                      !any(j -> oracle_separated(M[:, setdiff(v.terms, [j])]), v.terms)
            end
            ok &= !truth || !isempty(v.units)
            t[3] += !ok
        end
        for (fam, (total, separated, wrong)) in tally
            @test wrong == 0
            @test 0.25 * total < separated < 0.75 * total    # both classes exercised
        end
    end

    @testset "Separation verdict: the likelihood confirms it" begin
        # A separated verdict's direction is one along which the family's own
        # log-likelihood increases without decreasing; a design without
        # separation has a finite maximiser that the shared Newton optimizer
        # reaches with a vanishing gradient.
        newton_fit = NetworkCore.newton_fit
        logistic_derivatives = NetworkCore.logistic_derivatives
        rng = Xoshiro(7)
        n_sep = 0
        for _ in 1:150
            n = rand(rng, 6:20)
            X = hcat(ones(n), rand(rng, -2:2, n), rand(rng, 0:1, n))
            LinearAlgebra.rank(X) == 3 || continue    # collinear: not identified, not separated
            y = rand(rng, Bool, n)
            v = NetworkCore.logistic_separation(X, y)
            f = logistic_derivatives(X, y)
            if v.separated
                n_sep += 1
                ll = [first(f(t .* v.direction)) for t in (0.0, 4.0, 16.0, 64.0)]
                @test issorted(ll) && ll[end] > ll[1]
            else
                fit = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                    newton_fit(f, zeros(3))
                end
                @test fit.converged
                @test norm(f(fit.θ)[2]) < 1e-6
                # a finite maximum: the likelihood falls along every direction
                ll̂ = first(f(fit.θ))
                dirs = vcat([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]], [randn(rng, 3) for _ in 1:10])
                @test all(u -> first(f(fit.θ .+ 50 .* u ./ norm(u))) < ll̂ - 1e-3 &&
                               first(f(fit.θ .- 50 .* u ./ norm(u))) < ll̂ - 1e-3, dirs)
            end
        end
        @test 20 < n_sep < 130
        # Poisson: the rate on the zero-count rows goes to 0 along the direction
        Xp = [1.0 0.0; 1.0 0.0; 1.0 1.0; 1.0 1.0]
        pll(θ, y) = sum(y .* (Xp * θ) .- exp.(Xp * θ))
        v = NetworkCore.poisson_separation(Xp, [3, 1, 0, 0])
        @test v.separated && v.terms == [2] && v.units == [3, 4]
        @test pll(10 .* v.direction, [3, 1, 0, 0]) > pll(zeros(2), [3, 1, 0, 0])
        # conditional logit: the case's share of every stratum rises along d
        Xc = reshape([1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0], :, 1)
        ch = Bool[1, 0, 0, 1, 0, 0, 1, 0, 0]; st = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        cll(θ) = sum(let idx = findall(==(s), st)
                         c = idx[findfirst(i -> ch[i], idx)]
                         (Xc[c, :]' * θ) - log(sum(exp.(Xc[idx, :] * θ)))
                     end for s in 1:3)
        v = NetworkCore.clogit_separation(Xc, ch, st)
        @test v.separated && v.units == [1, 2, 3] && v.family === :clogit
        @test cll(8 .* v.direction) > cll(2 .* v.direction) > cll([0.0])
    end

    @testset "Separation verdict: designs every fitter meets" begin
        lsep = NetworkCore.logistic_separation
        x = collect(-5.0:1.0:5.0)
        # complete and quasi-complete separation; the slope alone carries it
        v = lsep(hcat(ones(11), x), x .> 0)
        @test v.separated && v.certified && v.terms == [2] && v.direction == [0.0, 1.0]
        @test lsep(hcat(ones(6), [-2.0, -1, 0, 0, 1, 2]), Bool[0, 0, 0, 1, 1, 1]).terms == [2]
        # a zero cell (a nodematch with no within-group tie), binary and binomial
        @test lsep(hcat(ones(7), [0, 0, 0, 0, 1, 1, 1]), Bool[0, 1, 0, 1, 1, 1, 1]).terms == [2]
        @test lsep([1.0 0.0; 1.0 1.0], [40.0, 10.0], [8.0, 0.0]).units == [2]
        @test !lsep([1.0 0.0; 1.0 1.0], [40.0, 10.0], [8.0, 3.0]).separated
        # a zero cell next to a covariate on a scale 10⁴ larger
        rng = Xoshiro(2)
        z = 1e4 .* randn(rng, 500); xb = rand(rng, 500) .< 0.2
        yb = (rand(rng, 500) .< 1 ./ (1 .+ exp.(-(z ./ 1e4)))) .| xb
        v = lsep(hcat(ones(500), z, xb), yb)
        @test v.separated && v.terms == [3] && v.certified
        # the same design at any column scale gives the same verdict
        for s in (1e-6, 1e6)
            @test lsep(hcat(ones(500), s .* z, s .* xb), yb).terms == [3]
        end
        # finite maxima are never flagged, however skewed
        @test !lsep(ones(1, 1), [1e9], [1.0]).separated
        y3 = x .> 0; y3[5] = true                       # one overlap: finite
        @test !lsep(hcat(ones(11), x), y3).separated
        for n in (20, 200, 2000), scale in (1.0, 1e3)
            X = hcat(ones(n), scale .* randn(rng, n), randn(rng, n))
            yy = rand(rng, n) .< 1 ./ (1 .+ exp.(-(X * [-1, 1.5 / scale, 0.5])))
            @test !lsep(X, yy).separated
        end
        # collinear columns are non-identification, not separation
        @test !lsep(hcat(ones(6), [1.0, 2, 3, 1, 2, 3], [2.0, 4, 6, 2, 4, 6]),
                    Bool[1, 0, 1, 0, 1, 0]).separated
        # weights: a zero-weight row carries no information
        @test lsep(hcat(ones(4), [0.0, 0, 1, 1]), Bool[0, 1, 1, 0]; weights=[1, 1, 1, 0]).terms == [2]
        @test !lsep(hcat(ones(4), [0.0, 0, 1, 1]), Bool[0, 1, 1, 0]).separated
        # Poisson: a zero count at one level of a factor
        psep = NetworkCore.poisson_separation
        Xp = [1.0 0.0; 1.0 0.0; 1.0 1.0; 1.0 1.0]
        @test psep(Xp, [3, 1, 0, 0]).terms == [2]
        @test !psep(Xp, [3, 1, 0, 2]).separated
        @test psep(ones(3, 1), [0, 0, 0]).separated           # a rate of zero: log-rate −Inf
        # conditional logit: one case per stratum (a relational event), and a
        # stratum with two observed values that must tie (a count dyad class)
        csep = NetworkCore.clogit_separation
        X = reshape([1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0], :, 1)
        @test !csep(X, Bool[1, 0, 0, 1, 0, 0, 1, 0, 0], [1, 1, 1, 2, 2, 2, 3, 3, 3]).separated
        Xs = reshape([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], :, 1); sts = [1, 1, 1, 2, 2, 2]
        @test csep(Xs, Bool[1, 0, 0, 1, 0, 0], sts).terms == [1]   # every observed count is 0
        @test !csep(Xs, Bool[1, 1, 0, 1, 0, 0], sts).separated      # two observed values must tie
        # weights remove an alternative from the choice set (strata labels kept)
        @test csep(reshape([1.0, 0.0, 5.0], :, 1), Bool[1, 0, 0], [7, 7, 7];
                   weights=[1, 1, 0]).units == [7]
        @test_throws ArgumentError csep(ones(2, 1), Bool[0, 0], [1, 1])   # no chosen alternative
        # a REM-sized design (481 events × 1332 alternatives) stays fast
        R = 1332; ne_ = 481
        Xb = randn(Xoshiro(3), ne_ * R, 4); chb = falses(ne_ * R); chb[1:R:end] .= true
        stb = repeat(1:ne_, inner=R)
        @test !csep(Xb, chb, stb).separated
        Xb[:, 2] .= chb                               # an indicator that is 1 on the case only
        v = csep(Xb, chb, stb)
        @test v.separated && v.terms == [2] && length(v.units) == ne_ && v.certified
        # the core criterion and its validation
        @test NetworkCore.separation_from_margins([1.0 0.0; -1.0 1.0]).terms == [2]
        @test !NetworkCore.separation_from_margins(reshape([1.0, -1.0], 2, 1)).separated
        @test !NetworkCore.separation_from_margins(zeros(3, 2)).separated
        @test !NetworkCore.separation_from_margins(zeros(0, 2)).separated
        @test_throws ArgumentError NetworkCore.separation_from_margins([1.0 Inf])
        @test_throws ArgumentError lsep(ones(2, 1), [1.0, 1.0], [2.0, 0.0])
        @test_throws ArgumentError lsep(ones(3, 1), Bool[1, 0])
        @test_throws ArgumentError psep(ones(2, 1), [-1, 2])
        # show
        @test occursin("separated, family=:logistic, terms=[2]", sprint(show, lsep(hcat(ones(11), x), x .> 0)))
        @test sprint(show, lsep(ones(1, 1), [2.0], [1.0])) == "SeparationVerdict(not separated, family=:logistic)"
    end

    @testset "Separation verdict: the certificate stays cheap and stays a proof" begin
        # The verdict reduces the direction to a minimal support in floating
        # point and certifies only that, so a quasi-separated design (most
        # margins exactly zero) costs about what the linear programme costs.
        # Before this, the exact certificate eliminated every near-zero row of
        # the full design in Rational{BigInt}: 13–51 s and 11–34 GiB on the
        # REM-sized designs below, and 6 s / 5.8 GiB on the logistic one.
        csep = NetworkCore.clogit_separation
        lsep = NetworkCore.logistic_separation
        csep(randn(Xoshiro(1), 30, 2), repeat(Bool[1, 0, 0], 10), repeat(1:10, inner=3))
        lsep(hcat(ones(40), [zeros(30); ones(10)]), [rand(Xoshiro(1), Bool, 30); trues(10)])
        # 481 events × 1332 alternatives, 8 statistics; the last is a
        # receiver indicator that only controls carry (officers never called)
        R = 1332; ne_ = 481
        X = randn(Xoshiro(3), ne_ * R, 8); ch = falses(ne_ * R); ch[1:R:end] .= true
        st = repeat(1:ne_, inner=R)
        X[:, 8] .= 0
        for e in 1:ne_, j in 2:4
            X[(e - 1) * R + j, 8] = 1.0
        end
        t = @timed csep(X, ch, st)
        v = t.value
        @test v.separated && v.certified && v.terms == [8] && length(v.units) == ne_
        @test v.direction == [0, 0, 0, 0, 0, 0, 0, -1.0]
        @test t.time < 6.0
        @test t.bytes < 500 * 2^20
        # the same with the indicator on ten cases only (quasi-complete)
        X[:, 8] .= 0
        X[(0:9) .* R .+ 1, 8] .= 1.0
        t = @timed csep(X, ch, st)
        @test t.value.separated && t.value.certified && t.value.terms == [8]
        @test t.value.units == collect(1:10)
        @test t.time < 6.0 && t.bytes < 500 * 2^20
        # logistic, integer ERGM-like statistics plus a zero cell
        rng = Xoshiro(4)
        Xl = Float64.(hcat(ones(10_000), rand(rng, 0:6, 10_000, 28), zeros(10_000)))
        yl = rand(rng, 10_000) .< 0.3
        Xl[findall(yl)[1:20], end] .= 1.0
        t = @timed lsep(Xl, yl)
        @test t.value.separated && t.value.certified && t.value.terms == [30]
        @test t.time < 3.0 && t.bytes < 200 * 2^20

        # Still a proof: the brute-force oracle on designs the small fuzz does
        # not reach: many duplicated rows (so a large near-zero set with few
        # distinct rows), dyadic fractions, power-of-two column scales and an
        # exactly collinear column. The oracle runs on the distinct rows.
        rng = Xoshiro(20261007)
        wrong = 0; n_sep = 0
        for _ in 1:300
            q = rand(rng, 1:4); nu = rand(rng, 2:7)
            U = rand(rng, -3:3, nu, q) .// rand(rng, (1, 2, 8), nu, q)
            rand(rng) < 0.5 && (U[:, 1] .= 1)
            q >= 3 && rand(rng) < 0.3 && (U[:, q] = U[:, 1] .+ 2 .* U[:, 2])
            scale = [2.0^rand(rng, -40:40) for _ in 1:q]
            y = rand(rng, Bool, nu)
            Mu = [(y[r] ? 1 : -1) * U[r, j] for r in 1:nu, j in 1:q]
            reps = rand(rng, 1:40, nu)                   # each distinct row repeated
            rows = reduce(vcat, [fill(r, reps[r]) for r in 1:nu])
            X = Float64.(U[rows, :]) .* scale'
            v = lsep(X, y[rows])
            truth = oracle_separated(Mu)
            n_sep += truth
            ok = v.separated == truth
            if ok && truth
                d = v.direction .* scale                 # in the units of U: X d = U (scale .* d)
                Mf = Float64.(Mu)
                margins = Mf * d
                mags = abs.(Mf) * abs.(d)                # relative: the scales span 2^±40
                ok &= v.certified && all(k -> margins[k] >= -1e-9 * mags[k], eachindex(margins)) &&
                      any(k -> margins[k] > 1e-9 * mags[k], eachindex(margins))
                ok &= length(v.terms) == 1 ||
                      !any(j -> oracle_separated(Mu[:, setdiff(v.terms, [j])]), v.terms)
            end
            wrong += !ok
        end
        @test wrong == 0
        @test 60 < n_sep < 240

        # The exact margin sign: Int128 when the scaled ray fits, BigInt
        # otherwise, both equal to Rational{BigInt} arithmetic, including
        # sums that cancel exactly and entries 2^±600 apart.
        rng = Xoshiro(9)
        for _ in 1:500
            s = rand(rng, 1:5)
            row = [rand(rng) < 0.3 ? 0.0 : ldexp(Float64(rand(rng, -9:9)), rand(rng, -60:60))
                   for _ in 1:s]
            N = BigInt[rand(rng, -10^6:10^6) for _ in 1:s]
            if s >= 2 && rand(rng) < 0.5                  # force an exact cancellation
                row[1] = ldexp(3.0, 5); row[2] = ldexp(-3.0, 4); N[1] = 1; N[2] = 2
                for j in 3:s; N[j] = 0; end
            end
            A = reshape(row, 1, s)
            exact = sign(sum(Rational{BigInt}(row[j]) * N[j] for j in 1:s))
            @test NetworkCore._exact_sign_bigint(A, 1, collect(1:s), N) == exact
            r128 = NetworkCore._exact_sign_int128(A, 1, collect(1:s), Int64.(N))
            @test r128 === nothing || r128 == exact
        end
        A = [ldexp(1.0, 600) -ldexp(1.0, -600) ldexp(1.0, 600)]
        N = BigInt[1, big(2)^1200, -1]
        @test NetworkCore._exact_sign_int128(A, 1, [1, 3], Int64[1, -1]) == 0
        @test NetworkCore._exact_sign_int128(A, 1, [1, 2], Int64[0, 1]) == -1
        @test NetworkCore._exact_sign_int128(A, 1, [1, 2], Int64[1, typemax(Int64)]) === nothing
        @test NetworkCore._exact_sign_bigint(A, 1, [1, 2, 3], N) == -1
        @test all(x -> ldexp(Float64(first(NetworkCore._dyadic(x))), last(NetworkCore._dyadic(x))) == x,
                  [0.1, -3.0, 5e-324, floatmax(Float64), 1.5e-300, -7.25])
        # the big-integer budget: past it the certificate is abandoned, so the
        # verdict says certified=false rather than running on
        uf = Float64[1.0, 0.0, -1.0]
        @test NetworkCore._margin_sign(A, 1, [1, 2, 3], uf, false, N, nothing, Ref(0)) === nothing
        @test NetworkCore._margin_sign(A, 1, [1, 2, 3], uf, false, N, nothing, Ref(1)) == -1
        @test NetworkCore._EXACT_BUDGET == 200_000
    end

    @testset "Separation policy: one warning, one caveat, inference withheld" begin
        x = collect(-3.0:3.0)
        v = NetworkCore.logistic_separation(hcat(ones(7), x), x .> 0)
        none = NetworkCore.logistic_separation(ones(1, 1), [2.0], [1.0])
        msg = NetworkCore.separation_message("my_fit", v, ["(Intercept)", "x"]; estimate="MPLE",
                                             note="R warns too.")
        @test startswith(msg, "my_fit: the MPLE does not exist (separation).")
        @test occursin("`x`", msg) && !occursin("`(Intercept)`", msg)
        @test occursin("6 observations perfectly", msg)
        @test occursin("`converged == false`", msg) && occursin("withheld", msg)
        @test endswith(msg, " R warns too.")
        # The objective is named after the estimate: a pseudo-likelihood is not a likelihood.
        @test occursin("The log-pseudo-likelihood keeps increasing", msg)
        @test occursin("The log-pseudo-likelihood", NetworkCore.separation_message("f", v; estimate="swap-MPLE"))
        @test occursin("The log-likelihood keeps increasing", NetworkCore.separation_message("f", v))
        @test occursin("coefficient 2", NetworkCore.separation_message("f", v))
        @test occursin("strata", NetworkCore.separation_message("f",
            NetworkCore.clogit_separation(reshape([1.0, 0, 1, 0], :, 1), Bool[1, 0, 1, 0], [1, 1, 2, 2])))
        @test (@test_logs (:warn, r"my_fit: the maximum-likelihood estimate does not exist") NetworkCore.warn_separation("my_fit", v, ["a", "x"])) === true
        @test (@test_logs NetworkCore.warn_separation("my_fit", none)) === false
        c = NetworkCore.separation_caveat(v, ["(Intercept)", "x"])
        @test startswith(c, "separation: the coefficient(s) on `x` are infinite")
        @test NetworkCore.separation_caveat(none) === nothing
        @test startswith(NetworkCore.separation_caveat(["x"]),
                         "separation: the coefficient(s) on `x` are infinite")
        @test NetworkCore.separation_caveat(String[]) === nothing
        zp = z_pvalues([1.0, 20.0], [0.5, 900.0], v)
        @test all(isnan, zp.z) && all(isnan, zp.p)
        @test z_pvalues([1.0, 2.0], [0.5, 0.5], none) == z_pvalues([1.0, 2.0], [0.5, 0.5])
        @test_throws DimensionMismatch z_pvalues([1.0], [0.5, 0.5], v)
    end

    @testset "newton_fit says why it stopped" begin
        newton_fit = NetworkCore.newton_fit
        logistic_derivatives = NetworkCore.logistic_derivatives
        quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
        # converged
        poi = newton_fit(θ -> (7θ[1] - exp(θ[1]), [7 - exp(θ[1])], hcat(-exp(θ[1]))), [0.0])
        @test poi.converged && poi.stop === :converged
        # the iteration limit, and only then :maxiter
        lim = newton_fit(θ -> (7θ[1] - exp(θ[1]), [7 - exp(θ[1])], hcat(-exp(θ[1]))), [0.0];
                         maxiter=1)
        @test !lim.converged && lim.stop === :maxiter && lim.iterations == 1
        # a column of zeros: the information is singular at the start, there is
        # no Newton direction, and the fit stops at iteration 1 — which used to
        # be indistinguishable from an iteration limit
        rng = Xoshiro(3)
        Xz = hcat(ones(20), zeros(20)); yz = rand(rng, 20) .< 0.5
        z = quiet(() -> newton_fit(logistic_derivatives(Xz, yz), zeros(2)))
        @test !z.converged && z.stop === :singular_start && z.iterations == 1
        # exactly collinear columns, likewise
        Xc = hcat(ones(6), [1.0, 2, 3, 1, 2, 3], [2.0, 4, 6, 2, 4, 6])
        c = quiet(() -> newton_fit(logistic_derivatives(Xc, Bool[1, 0, 1, 0, 1, 1]), zeros(3)))
        @test !c.converged && c.stop === :singular_start
        # collinear up to rounding: it may move, but it never reports an
        # iteration limit or convergence
        for _ in 1:5
            x = randn(rng, 200)
            r = quiet(() -> newton_fit(logistic_derivatives(hcat(ones(200), x, 3 .* x .+ 1),
                                                            rand(rng, 200) .< 0.4), zeros(3)))
            @test !r.converged && r.stop in (:singular_start, :singular, :singular_information)
        end
        # settled with a singular information at the solution
        flat = quiet(() -> newton_fit(θ -> (-(θ[1] - 1)^2, [-2 * (θ[1] - 1), 0.0],
                                            [-2.0 0.0; 0.0 0.0]), [0.0, 0.0]))
        @test !flat.converged && flat.stop in (:singular_start, :singular_information)
        # no ascent: every step decreases the objective (a minimum, not a maximum)
        bowl = quiet(() -> newton_fit(θ -> (θ[1]^2 + 1, [2θ[1]], fill(2.0, 1, 1)), [1.0]))
        @test !bowl.converged && bowl.stop === :no_ascent
        # the positional destructuring older callers use still works
        θ̂, se, V, ll, conv, it = poi
        @test conv && it == poi.iterations
    end

    @testset "Shared task helpers: spawn_all and unwrap_task_exception" begin
        spawn_all = NetworkCore.spawn_all
        unwrap = NetworkCore.unwrap_task_exception
        out = zeros(Int, 8)
        @test spawn_all(c -> (out[c] = c^2), 8) === nothing
        @test out == (1:8) .^ 2
        @test spawn_all(c -> error("never"), 0) === nothing
        # the first failing task's OWN exception, after every task finished
        done = Threads.Atomic{Int}(0)
        err = try
            spawn_all(8) do c
                c in (3, 6) && throw(ArgumentError("task $c"))
                sleep(0.01); Threads.atomic_add!(done, 1)
            end
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && err.msg == "task 3"
        @test done[] == 6
        # The failing task's backtrace survives: the TaskFailedException is the
        # cause on the exception stack, and its task holds the frames inside
        # the task (here `_spawn_all_inner_fail`), which an uncaught error prints.
        stack = try
            spawn_all(c -> c == 2 ? _spawn_all_inner_fail(c) : nothing, 3)
            nothing
        catch
            current_exceptions()
        end
        @test stack !== nothing && length(stack) == 2
        @test stack[end].exception isa ArgumentError &&
              stack[end].exception.msg == "inner failure 2"
        cause = stack[1].exception
        @test cause isa TaskFailedException
        @test cause.task.exception === stack[end].exception
        @test any(fr -> fr.func === :_spawn_all_inner_fail, stacktrace(cause.task.backtrace))
        shown = sprint(Base.display_error, stack)
        @test occursin("_spawn_all_inner_fail", shown) && occursin("caused by", shown)
        # nested wrappers are peeled to the root cause
        t = Threads.@spawn throw(DomainError(-1, "negative"))
        e = try wait(t) catch x; x end
        @test e isa TaskFailedException && unwrap(e) isa DomainError
        @test unwrap(CompositeException([e])) isa DomainError
        @test unwrap(CompositeException([])) isa CompositeException
        ex = ErrorException("x")
        @test unwrap(ex) === ex
        # bootstrap_cov uses the same helper
        @test_throws DomainError bootstrap_cov(_ -> throw(DomainError(0)), (rng, B) -> fill(1, B),
                                               [0.0]; n_boot=4, rng=Xoshiro(1))
    end

    @testset "DataFrames-extension functions explain themselves without DataFrames" begin
        # A bare MethodError used to be all the user got. The test
        # environment loads DataFrames (the extension is active here), so the
        # hint is checked in a fresh process that does not.
        code = """
            using NetworkCore
            for f in (() -> as_dataframe(network(3)), () -> network_from_dataframe(1))
                try f() catch e; print(sprint(showerror, e)); println("|") end
            end
            """
        out = read(`$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $code`, String)
        @test count("NetworkCore.jl's DataFrames extension, which is not loaded", out) == 2
        @test occursin("using DataFrames", out)
        # with DataFrames loaded (as here) the hint stays silent
        msg = sprint(showerror, MethodError(as_dataframe, (1,)))
        @test !occursin("not loaded", msg)
    end

    @testset "Dataset directory is resolved at run time" begin
        # A precompile-time `const` path broke relocated depots and
        # sysimages; the directory is now asked of the loaded module.
        @test !isdefined(NetworkCore, :_DATA_DIR)
        @test NetworkCore._data_dir() == pkgdir(NetworkCore, "data")
        @test isfile(joinpath(NetworkCore._data_dir(), "florentine_marriage_edges.tsv"))
    end

    @testset "Graphs.jl algorithms on the Network, never on net.graph" begin
        # The undirected backing store holds both arcs, so Graphs.jl
        # on `net.graph` double-counts. The documented routes — the Network
        # itself, or `SimpleGraph(net)` — agree with a hand-built SimpleGraph.
        flo = load_dataset(:florentine_marriage)
        g = SimpleGraph(nv(flo))
        for e in edges(flo); add_edge!(g, src(e), dst(e)); end
        @test SimpleGraph(flo) == g
        @test Graphs.core_number(flo) == Graphs.core_number(g)
        @test Graphs.degree(flo) == Graphs.degree(g)
        @test Graphs.betweenness_centrality(flo) ≈ Graphs.betweenness_centrality(g)
        @test Graphs.global_clustering_coefficient(flo) ≈ Graphs.global_clustering_coefficient(g)
        # ... and the trap the docs warn about is real
        @test Graphs.core_number(flo.graph) == 2 .* Graphs.core_number(g)
        @test Graphs.ne(flo.graph) == 2 * ne(flo)
    end

    @testset "read_pajek: list sections, self-loops, two-mode header" begin
        dir = mktempdir()
        rd(lines) = (p = joinpath(dir, "f$(rand(UInt32)).net"); write(p, join(lines, "\n")); read_pajek(p))
        el(net) = [(src(e), dst(e)) for e in edges(net)]
        # *Arcslist: `1 2 3 4` is three arcs (it used to be the arc 1→2 of weight 3)
        a = rd(["*Vertices 4", "*Arcslist", "1 2 3 4", "2 3"])
        @test is_directed(a) && el(a) == [(1, 2), (1, 3), (1, 4), (2, 3)]
        @test isempty(list_edge_attributes(a))
        e = rd(["*Vertices 4", "*Edgeslist", "1 2 3", "4 2"])
        @test !is_directed(e) && el(e) == [(1, 2), (1, 3), (2, 4)]
        # self-loops survive (they were dropped: the network had loops=false)
        l = rd(["*Vertices 3", "*Arcs", "1 1", "1 2"])
        @test l.loops && el(l) == [(1, 1), (1, 2)]
        lw = rd(["*Vertices 3", "*Arcs", "1 1 2.5", "1 2 1"])        # used to throw
        @test get_edge_attribute(lw, :weight, 1, 1) == 2.5
        # the two-mode header is honoured, and write_pajek writes it back
        b = rd(["*Vertices 5 2", "*Edges", "1 3", "2 5"])
        @test b.bipartite == 2 && el(b) == [(1, 3), (2, 5)]
        p = joinpath(dir, "b.net"); write_pajek(b, p)
        @test first(readlines(p)) == "*Vertices 5 2"
        rb = read_pajek(p)
        @test rb.bipartite == 2 && el(rb) == el(b)
        # nothing is dropped silently
        @test_throws ArgumentError rd(["*Vertices 5 2", "*Edges", "1 2"])        # within mode 1
        @test_throws ArgumentError rd(["*Vertices 3", "*Arcs", "1 7"])           # outside 1:n
        @test_throws ArgumentError rd(["*Vertices 2", "*Matrix", "0 1", "1 0"])  # unsupported section
        d = @test_logs (:warn, r"listed more than once") rd(["*Vertices 3", "*Edges", "1 2 1", "2 1 4"])
        @test ne(d) == 1 && get_edge_attribute(d, :weight, 1, 2) == 4.0
        # a *Network title line is skipped
        @test ne(rd(["*Network test", "*Vertices 2", "*Arcs", "1 2"])) == 1
    end

    @testset "Aqua.jl quality assurance" begin
        # ambiguities, unbound type parameters, undefined exports, stale deps,
        # [compat] for every dependency AND test-only extra, type piracy,
        # persistent tasks, undocumented names
        Aqua.test_all(NetworkCore)
        Aqua.test_undocumented_names(NetworkCore)
    end


    @testset "Every docstring example runs" begin
        # Each ```julia block of an exported or public docstring runs in a
        # fresh module that has done `using NetworkCore` (as Documenter's
        # DocTestSetup does); anything else it needs, it loads itself. Last in
        # the file, because some examples extend the shared generics.
        @test length(DOC_BLOCKS) >= 120
        failed = String[]
        for (nm, code) in DOC_BLOCKS
            mod = Module(Symbol("DocExample_", nm))
            try
                redirect_stdout(devnull) do
                    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                        Core.eval(mod, :(using NetworkCore))
                        Base.include_string(mod, code, "docstring:$nm")
                    end
                end
            catch err
                push!(failed, "$nm: " * sprint(showerror, err))
            end
        end
        foreach(msg -> println(stderr, "docstring example failed: ", msg), failed)
        @test isempty(failed)
    end
end
