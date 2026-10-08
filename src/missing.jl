# Missing-dyad (unobserved tie) mask.
#
# A masked dyad means the tie status is UNOBSERVED — the analyst does not
# know whether the tie exists. This is distinct from "no tie": the backing
# graph still stores a face value (edge present or absent) for every dyad,
# and all structural queries (`has_edge`, `ne`, `degree`, `as_matrix`, ...)
# report that face value. The mask is metadata that statistical packages
# (ERGM.jl, ...) consult to exclude unobserved dyads from estimation, in the
# spirit of R statnet's NA edges.
#
# Storage: a `Set{Tuple{T,T}}` on the `Network` struct, keyed by
# `_canonical_edge` (ordered pair for directed networks, `minmax` pair for
# undirected) so lookups are O(1) in samplers' hot paths and undirected
# handling is symmetric.

"""
    set_missing_dyad!(net::Network, i, j) -> Network

Mark the dyad `(i, j)` as **missing (unobserved)**: the analyst does not know
whether the tie exists. This does *not* add or remove an edge — the network
keeps a "face value" for the dyad (edge present or absent), which structural
queries continue to report. Statistical packages built on NetworkCore.jl use the
mask to exclude unobserved dyads from estimation (e.g. ERGM.jl drops masked
dyads from the MPLE design matrix and never toggles them during MCMC).

For undirected networks the mask is symmetric: marking `(i, j)` also masks
`(j, i)`. For directed networks each ordered pair (arc) is masked
independently. Marking an already-masked dyad is a no-op.

Throws `ArgumentError` if either vertex does not exist, or for a self-dyad
`(i, i)` when the network does not allow loops.

See also [`is_missing_dyad`](@ref), [`delete_missing_dyad!`](@ref),
[`clear_missing_dyads!`](@ref), [`missing_dyads`](@ref),
[`n_missing_dyads`](@ref).

# Example
```julia
net = network(4)
add_edge!(net, 1, 2)
set_missing_dyad!(net, 3, 4)   # tie 3→4 is unobserved
is_missing_dyad(net, 3, 4)     # true
has_edge(net, 3, 4)            # false — the face value, NOT an observation
```
"""
function set_missing_dyad!(net::Network{T}, i, j) where T
    has_vertex(net, i) && has_vertex(net, j) ||
        throw(ArgumentError("both vertices of a missing dyad must exist " *
                            "(got ($i, $j) with $(nv(net)) vertices)"))
    if i == j && !net.loops
        throw(ArgumentError("cannot mark the self-dyad ($i, $i) as missing: " *
                            "the network does not allow self-loops"))
    end
    if !isnothing(net.bipartite)
        k = net.bipartite
        ((i <= k) ⊻ (j <= k)) ||
            throw(ArgumentError("cannot mark the within-mode dyad ($i, $j) of " *
                                "a two-mode network as missing: within-mode " *
                                "ties are structurally impossible, not unobserved"))
    end
    push!(net.missing_dyads, _canonical_edge(net, T(i), T(j)))
    return net
end

"""
    is_missing_dyad(net::Network, i, j) -> Bool

Return whether the tie status of dyad `(i, j)` is missing (unobserved); see
[`set_missing_dyad!`](@ref). O(1). For undirected networks the check is
symmetric in `i` and `j`.

# Example
```julia
net = network(3)
set_missing_dyad!(net, 1, 2)
is_missing_dyad(net, 1, 2)    # true
is_missing_dyad(net, 2, 1)    # false — directed: each arc is masked separately
```
"""
is_missing_dyad(net::Network{T}, i, j) where T =
    _canonical_edge(net, T(i), T(j)) in net.missing_dyads

"""
    delete_missing_dyad!(net::Network, i, j) -> Bool

Remove the missing-dyad mark from dyad `(i, j)`, declaring its tie status
observed again (at its current face value). Returns `true` if the dyad was
masked, `false` if it was not.

# Example
```julia
net = network(3)
set_missing_dyad!(net, 1, 2)
delete_missing_dyad!(net, 1, 2)    # true
delete_missing_dyad!(net, 1, 2)    # false — no longer masked
n_missing_dyads(net)               # 0
```
"""
function delete_missing_dyad!(net::Network{T}, i, j) where T
    key = _canonical_edge(net, T(i), T(j))
    present = key in net.missing_dyads
    delete!(net.missing_dyads, key)
    return present
end

"""
    clear_missing_dyads!(net::Network) -> Network

Remove every missing-dyad mark, declaring all dyads observed. This is the
route a routine that refuses masked data (`missing=:error`) points at when
the analyst has decided the face values are the data.

# Example
```julia
net = network(3)
set_missing_dyad!(net, 1, 2)
set_missing_dyad!(net, 2, 3)
clear_missing_dyads!(net)
n_missing_dyads(net)    # 0
```
"""
function clear_missing_dyads!(net::Network)
    empty!(net.missing_dyads)
    return net
end

"""
    missing_dyads(net::Network{T}) -> iterator of Tuple{T,T}

Iterator over the masked (unobserved) dyads of the network, in canonical
form: `(i, j)` as marked for directed networks, `minmax(i, j)` for
undirected. Iteration order is unspecified. The iterator reflects the live
mask; do not add or delete marks while iterating. The Julia counterpart of
R's `is.na(net)` / `as.matrix(net) |> is.na`.

# Example
```julia
net = network(3; directed=false)
set_missing_dyad!(net, 3, 1)
collect(missing_dyads(net))    # [(1, 3)] — canonical minmax form
```
"""
missing_dyads(net::Network) = (d for d in net.missing_dyads)

"""
    n_missing_dyads(net::Network) -> Int

Number of masked (unobserved) dyads: ordered pairs for directed networks,
unordered pairs for undirected (R's `network.naedgecount`).

# Example
```julia
net = network(4)
set_missing_dyad!(net, 1, 2)
set_missing_dyad!(net, 3, 4)
n_missing_dyads(net)    # 2
```
"""
n_missing_dyads(net::Network) = length(net.missing_dyads)

# ============================================================================
# The ecosystem missing-data contract
# ============================================================================
#
# A masked dyad is *unobserved*, not *absent*. Every routine in the ecosystem
# therefore has to make a choice about it, and the dangerous choice is the
# implicit one: reading the face value and returning a plausible number.
#
# The contract has two halves:
#
#   `supports_missing(f)` — a trait a routine opts into by defining a method.
#       The default is `false`, so a routine that has never thought about
#       missingness is, correctly, one that refuses to see it.
#
#   `require_observed(net, policy; context)` — the guard the routine calls.
#       With the default `policy = :error` a masked network raises rather than
#       being silently read at face value; `policy = :face` is the explicit,
#       auditable opt-in to the old behaviour.
#
# The point is that no exported routine may return a number computed from the
# face value of an unobserved tie without the caller having said so in writing.

"""
    supports_missing(x) -> Bool

Whether the routine, model type, or estimator `x` implements a principled
treatment of missing (unobserved) dyads.

This is the ecosystem's missing-data trait. It defaults to `false`: a routine
that has not declared a missing-data method is assumed *not* to have one, and
callers guard with [`require_observed`](@ref) so that a partially observed
network cannot be silently analysed at face value.

A routine opts in by defining a method:

```julia
my_density(net; missing=:error) = (require_observed(net, missing; context="my_density");
                                   ne(net) / (nv(net) * (nv(net) - 1)))
NetworkCore.supports_missing(::typeof(my_density)) = true
supports_missing(my_density)    # true
supports_missing(sum)           # false — the default
```

See also [`require_observed`](@ref), [`n_missing_dyads`](@ref).
"""
supports_missing(::Any) = false

"""
    MISSING_POLICIES

The missing-dyad policies accepted by [`require_observed`](@ref):

- `:error` — refuse to analyse a network with masked dyads (the default).
- `:face`  — explicitly opt in to using the stored face value of each masked
  dyad, i.e. treat "unobserved" as "as recorded". Auditable, never implicit.

# Example
```julia
MISSING_POLICIES              # (:error, :face)
:face in MISSING_POLICIES     # true
```
"""
const MISSING_POLICIES = (:error, :face)

"""
    missing_policies(f) -> Tuple{Vararg{Symbol}}

The missing-dyad policies **routine `f` actually accepts** through its
`missing=` keyword — the routine-level vocabulary of the missing-data contract.

[`MISSING_POLICIES`](@ref) is the vocabulary [`require_observed`](@ref)
understands; `missing_policies(f)` is the subset (or, for a package with its own
vocabulary, the *variant*) that a particular routine offers. Tooling such as the
ecosystem's capability-matrix generator prints `missing_policies(f)` for each
routine instead of assuming that `:face` is the opt-in everywhere — which it
deliberately is not.

- The default is `(:error,)`: a routine that exposes no `missing=` keyword
  refuses a masked network and offers no escape hatch. That is the correct
  answer for an estimator (TERGM, ERGMCount, ERGMMulti, ...) for which reading
  face values would enumerate unobserved dyads as observed rows.
- A routine that offers `missing=:face` declares
  `NetworkCore.missing_policies(::typeof(f)) = (:error, :face)` — as
  [`network_density`](@ref) does, and as every SNA descriptive measure with a
  `missing=` keyword is expected to.
- A package with its own vocabulary declares it here too: ERGM.jl's MCMC entry
  points accept `(:error, :condition_on_face)`, and `:face` is *not* among them,
  so the generator (and the error message) must not say that it is.

The declared tuple must start with `:error`, the policy every routine honours.

    missing_policies(f, ::Type{R}) -> Tuple{Vararg{Symbol}}

**Per-result-type form** for a *shared generic*. A verb such as [`gof`](@ref)
has one method per package (`gof(::ERGMResult)`, `gof(::STERGMResult)`, ...),
so a one-argument `missing_policies(gof)` cannot speak for any of them; a
package declares `NetworkCore.missing_policies(::typeof(gof), ::Type{<:MyResult})`
instead, and tooling asks the two-argument form for shared verbs. The fallback
`missing_policies(f, ::Type) = missing_policies(f)` makes the two-argument
query uniform: for a routine that owns its name it answers exactly what the
one-argument form does.

# Example
```julia
using NetworkCore
missing_policies(network_density)   # (:error, :face)
missing_policies(sum)               # (:error,) — no `missing=` keyword at all
missing_policies(network_density, Network)   # (:error, :face) — the fallback

# a routine declaring the policies it takes
function my_measure(net; missing::Symbol=:error)
    require_observed(net, missing; context="my_measure")
    # ...
end
NetworkCore.missing_policies(::typeof(my_measure)) = (:error, :face)

# a package's method of a shared generic declares per result type
struct MyResult end
NetworkCore.missing_policies(::typeof(gof), ::Type{MyResult}) = (:error, :condition_on_face)
missing_policies(gof, MyResult)     # (:error, :condition_on_face)
missing_policies(gof)               # (:error,) — the generic itself promises nothing
```

See also [`supports_missing`](@ref), [`require_observed`](@ref).
"""
missing_policies(::Any) = (:error,)
missing_policies(f, ::Type) = missing_policies(f)

"""
    require_observed(net, policy::Symbol=:error; context::AbstractString="this routine") -> net

Guard a routine against silently reading unobserved ties at face value.

Returns `net` unchanged when the network has no masked dyads, or when
`policy === :face` (the caller has explicitly accepted face values). Throws an
`ArgumentError` naming `context` when the network has masked dyads and
`policy === :error`.

`face_ok` says whether the CALLING routine actually exposes a `missing=` keyword
the user could set to `:face`. Most estimators do not, and must not: for
TERGM/ERGMCount/ERGMMulti a `:face` fit would enumerate unobserved dyads as
observed, which is exactly the silent data-invention this contract exists to
stop — so they call `require_observed(net; face_ok=false)` and the error omits an
instruction the caller could not follow. Suggesting a keyword that does not exist
is its own small lie.

`hint` is appended verbatim to the refusal, after the shared bullets: a routine
with its OWN opt-in vocabulary (ERGM's `missing=:condition_on_face` /
`missing=:mle`) adds the bullets naming the policies it does take, in the one
message, instead of catching the shared error and rethrowing it with text
glued on. Format it as the shared bullets are (`"  • ...\\n"`) so the list
reads as one.

This is the enforcement half of the ecosystem missing-data contract; see
[`supports_missing`](@ref) for the declaration half.

# Example
```julia
# a routine that CAN read face values, on request
function my_measure(net; missing::Symbol=:error)
    require_observed(net, missing; context="my_measure")
    # ... safe to read face values from here
end

# an estimator that cannot represent an unobserved dyad at all
function my_estimator(net)
    require_observed(net; context="my_estimator", face_ok=false)
end

# an estimator with its own opt-in, named in the same message
function my_mcmc_fit(net)
    require_observed(net; context="my_mcmc_fit", face_ok=false,
                     hint="  • pass `missing=:condition_on_face` to hold the " *
                          "masked dyads at their stored value during MCMC.\\n")
end
```
"""
function require_observed(net, policy::Symbol=:error;
                          context::AbstractString="this routine",
                          face_ok::Bool=true, hint::AbstractString="")
    policy in MISSING_POLICIES ||
        throw(ArgumentError("invalid missing-dyad policy $(repr(policy)); " *
                            "expected one of $(join(map(repr, MISSING_POLICIES), ", "))"))
    policy === :face && return net

    n = n_missing_dyads(net)
    n == 0 && return net

    face_hint = face_ok ?
        "  • pass `missing=:face` to analyse the stored face values anyway " *
        "(the tie is treated as recorded); or\n" : ""

    throw(ArgumentError(
        "$context does not support missing (unobserved) dyads, but the network " *
        "has $n masked $(n == 1 ? "dyad" : "dyads"). A masked dyad is unobserved, " *
        "not absent, so reading its face value would silently invent data.\n" *
        face_hint *
        "  • call `clear_missing_dyads!(net)` to declare every dyad observed; or\n" *
        "  • use a routine that implements missing-data handling " *
        "(`NetworkCore.supports_missing(f) == true`)." *
        (isempty(hint) ? "" : "\n" * hint)))
end

# BipartiteNetwork forwards the missing-dyad API to the wrapped Network.
set_missing_dyad!(bnet::BipartiteNetwork, i, j) =
    (set_missing_dyad!(bnet.network, i, j); bnet)
is_missing_dyad(bnet::BipartiteNetwork, i, j) = is_missing_dyad(bnet.network, i, j)
delete_missing_dyad!(bnet::BipartiteNetwork, i, j) =
    delete_missing_dyad!(bnet.network, i, j)
clear_missing_dyads!(bnet::BipartiteNetwork) =
    (clear_missing_dyads!(bnet.network); bnet)
missing_dyads(bnet::BipartiteNetwork) = missing_dyads(bnet.network)
n_missing_dyads(bnet::BipartiteNetwork) = n_missing_dyads(bnet.network)
