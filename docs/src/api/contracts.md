# Ecosystem Contracts

NetworkCore.jl is the foundation every model package in the ecosystem depends on,
so it is where the *shared contracts* live: the small vocabularies that let
`ERGM`, `Siena`, `REM`, `Revel`, `TERGM` and the rest mean the same thing by
the same keyword. A contract is defined once, here, and imported by name
wherever it is honoured — never redefined locally.

There are four, each documented below or linked from here:

| Contract | Keyword | Vocabulary |
|----------|---------|------------|
| [Missing data](#Missing-Data-Contract) | `missing=` | [`MISSING_POLICIES`](@ref) |
| [Tied events](#Tied-Event-Contract) | `ties=` | [`TIE_POLICIES`](@ref) |
| [Conversion](../guide/conversion_invariants.md) | `report=` | [`ConversionReport`](@ref) |
| [Separation](#Separation-Contract) | none: a fit property | [`SeparationVerdict`](@ref NetworkCore.SeparationVerdict) |

## Missing-Data Contract

A masked dyad is **unobserved**, not absent. The distinction is the whole point:
reading the stored face value of an unobserved tie silently invents data, and
turns a partially observed network into plausible, wrong output.

The contract has two halves. A routine *declares* whether it has a principled
missing-data treatment with [`supports_missing`](@ref), and *enforces* the
caller's choice with [`require_observed`](@ref).

The rule it enforces: **no exported routine may return a number computed from
the face value of an unobserved tie unless the caller asked for it in writing.**

```julia
function my_measure(net; missing::Symbol=:error)
    require_observed(net, missing; context="my_measure")
    # ... safe to read face values from here
end
```

Downstream packages expose this as the `missing=` keyword. `missing=:error` (the
default) refuses masked data; `missing=:face` is the explicit, auditable opt-in
to face values **where a routine offers it** — and not every routine does. An
estimator for which reading face values would enumerate unobserved dyads as
observed rows (TERGM, ERGMCount, ERGMMulti, ...) exposes no `missing=` keyword
at all, and ERGM.jl's MCMC entry points use their own vocabulary
(`:condition_on_face`).

So each routine *declares* the policies its keyword actually takes with
[`missing_policies`](@ref). The default `(:error,)` is the honest answer for a
routine without the keyword. The ecosystem's capability-matrix generator prints
`missing_policies(f)` for every routine it lists, instead of assuming `:face`
everywhere, and the shared error message raised by `require_observed` mentions
`missing=:face` only when the caller passed `face_ok=true`. A routine with its
own opt-in vocabulary appends its bullets to that one message through
`require_observed(...; hint=)` (ERGM's `missing=:condition_on_face` /
`missing=:mle`) rather than rethrowing it with text glued on.

For a *shared generic* — `gof`, which has one method per package — the
one-argument form cannot speak for any package's method, so packages declare
`NetworkCore.missing_policies(::typeof(gof), ::Type{<:TheirResult})` and tooling
asks the two-argument form `missing_policies(gof, ResultType)`; the fallback
`missing_policies(f, ::Type) = missing_policies(f)` makes that query uniform.

```@docs
supports_missing
missing_policies
require_observed
MISSING_POLICIES
```

See also the [Missing Data guide](../guide/missing_data.md) and the
missing-dyad mask API on the [Types page](types.md).

## Tied-Event Contract

For the models with an event clock (`REM.fit_rem`,
`Revel.fit_revel`), two events sharing a timestamp is not a nuisance to be
sorted away — an exact-time or exact-order likelihood evaluated on tied data is
evaluated on information the data does not contain.

So the vocabulary starts from refusal (`ties=:error` is the default) and every
alternative is an explicit, recorded claim about what a simultaneous pair
*means*. A policy a model cannot honour is refused loudly by
[`check_tie_policy`](@ref) rather than silently no-oping.

```@docs
TIE_POLICIES
check_tie_policy
```

[`tie_method`](@ref) reports the policy that actually **ran** — which is why
`:error` can never appear there: a fit that would have had to break a tie under
`ties=:error` threw instead of returning.

## Separation Contract

A logistic, Poisson or conditional-logit likelihood whose data are
**separated** has no finite maximiser. Some combination of the coefficients
can run to ±∞ and keep improving the fit, because along that direction a set
of observations is predicted perfectly. Newton's method does not notice. Its
increments shrink along the asymptote, and it stops at a large, finite point
with enormous standard errors and `converged == true`, exactly as R's `glm`
does. Every z value and p-value computed from that point is meaningless.

**One verdict.** Whether a finite maximiser exists is a property of the data,
not of where an optimizer stopped. Each likelihood family has one entry
point, and all three reduce to one criterion:

| Family | Entry point | Used by |
|--------|-------------|---------|
| logistic (binary or binomial rows) | [`logistic_separation`](@ref NetworkCore.logistic_separation) | ERGM's MPLE and its variants, TERGM's CMPLE, SNA's `netlogit` |
| Poisson log-linear | [`poisson_separation`](@ref NetworkCore.poisson_separation) | Poisson GLMs |
| conditional logit (strata of alternatives) | [`clogit_separation`](@ref NetworkCore.clogit_separation) | REM's case-control partial likelihood, relational-event models, ERGMCount's count pseudo-likelihood (a conditional logit over each dyad's support) |
| any of the above, as margin rows | [`separation_from_margins`](@ref NetworkCore.separation_from_margins) | the core criterion |

The criterion: the log-likelihood has a **direction of recession** `d`, one
along which it never decreases and somewhere increases. For these families
that is a set of linear inequalities on `d`, one per observation. A success
row needs `x'd ≥ 0` and a failure row `x'd ≤ 0`. A Poisson row needs
`x'd ≤ 0`, and `x'd = 0` if its count is positive. A chosen alternative needs
`(x_c − x_j)'d ≥ 0` against every alternative in its stratum. Separation means
all of them hold with at least one strict. The verdict decides this with a
linear programme, as R ergm's `mple.existence` does. It then certifies the
direction it found in exact arithmetic. A certified verdict is a proof,
independent of tolerances. The certificate is built on the minimal set of
coefficients, where most margins are zero, so the verdict costs about as much
as a few floating-point linear programmes: well under a second on a
relational-event design with 640,000 alternatives, whether or not it is
separated. A direction that leaves every margin at
exactly zero is a collinear design, not separation; the optimizer's rank
check reports that case.

The verdict names the coefficients to flag (`terms`). They form an
inclusion-minimal set that carries a direction of recession on its own, so a
design separated by one covariate names that covariate rather than every
coefficient the optimum happens to move. It also lists the observations
predicted perfectly (`units`).

**One policy.** It is decided once, here, and every fitting function follows
it, as R's `glm`, `coxph` and ergm do:

1. **warn**, with the one sentence of [`separation_message`](@ref NetworkCore.separation_message)
   (emitted by [`warn_separation`](@ref NetworkCore.warn_separation));
2. return the fit with **`converged == false`**;
3. **flag** the separated terms on the result;
4. **withhold inference**. z values, p-values and confidence intervals are
   `NaN`. [`z_pvalues`](@ref) with the verdict as third argument does this.
   [`separation_caveat`](@ref NetworkCore.separation_caveat) is the entry for
   `approximations(fit)`.

The estimates stay in the result, because they are what the user has to
diagnose with, but nothing that pretends they are a maximum is reported.

```julia
using NetworkCore
# A nodematch-like design: no tie in the second dyad class
X = [1.0 0.0; 1.0 1.0]
v = NetworkCore.logistic_separation(X, [40.0, 10.0], [8.0, 0.0])
v.separated, v.terms, v.certified           # (true, [2], true)
names = ["edges", "nodematch"]
NetworkCore.warn_separation("my_mple", v, names; estimate="MPLE")   # warns
NetworkCore.separation_caveat(v, names)     # the approximations() entry
z_pvalues([-1.39, -20.1], [0.4, 2900.0], v) # (z = [NaN, NaN], p = [NaN, NaN])
```

```@docs
NetworkCore.SeparationVerdict
NetworkCore.logistic_separation
NetworkCore.poisson_separation
NetworkCore.clogit_separation
NetworkCore.separation_from_margins
NetworkCore.separation_message
NetworkCore.warn_separation
NetworkCore.separation_caveat
```
