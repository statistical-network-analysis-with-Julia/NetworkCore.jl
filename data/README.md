# Bundled teaching datasets

## Licences at a glance

NetworkCore.jl's code is MIT-licensed; this directory holds data extracted from R
packages distributed under the GNU GPL, redistributed with the provenance
recorded below. Whether the extracted values carry the upstream licence is not
settled here; anyone redistributing them should consult the upstream terms.

| Files | Dataset | Taken from | Upstream licence |
|---|---|---|---|
| `florentine_*.tsv` | Padgett's Florentine families | statnet `ergm` (`florentine.RData`) | GPL-3 |
| `sampson_*.tsv` | Sampson's monastery | statnet `ergm` (`sampson.RData`) | GPL-3 |
| `faux_mesa_high_*.tsv` | faux.mesa.high (simulated) | statnet `ergm` (`faux.mesa.high.RData`) | GPL-3 |
| `wtc_police_calls_*.tsv` | WTC police radio calls | statnet `relevent` workshop (`relevent_workshop.Rdata`) | GPL-2+ |
| `s50/*.tsv` | Teenage Friends and Lifestyle Study, s50 excerpt | RSiena | GPL-2+ |

Tab-separated files loaded by `load_dataset(name::Symbol)` (see
`src/datasets.jl`). Every value was extracted programmatically from the
upstream R objects and cross-checked against this ecosystem's own golden
fixtures (SNA.jl, ERGM.jl, ERGMEgo.jl, Siena.jl). The Florentine and Sampson
files were parsed from the statnet `ergm` package's `data/florentine.RData`
and `data/sampson.RData` (https://github.com/statnet/ergm); the three
datasets added in 2026-09 are produced by the checked-in
script `r/extract_datasets.R`, which needs only base R (>= 4.5) and the
CRAN `network` package:

```bash
Rscript data/r/extract_datasets.R      # from the package root; downloads the two R objects
```

The script prints the SHA-256 of every upstream file it consumed and of every
TSV it wrote; the hashes below are from the run that produced the committed
files (R 4.6.1, network 1.20.0, 2026-09-09), so a regeneration can be checked
against them.

## Florentine families (Padgett)

- `florentine_vertices.tsv` — 16 families with `name`, `wealth` (net wealth
  in thousands of lira, 1427), `priorates` (seats on the civic council,
  1282–1344), `totalties` (total number of business or marriage ties in the
  full Padgett dataset of 116 families).
- `florentine_marriage_edges.tsv` — 20 undirected marriage ties
  (statnet `flomarriage`).
- `florentine_business_edges.tsv` — 15 undirected business ties
  (statnet `flobusiness`).

Source: Padgett, J. F. (1994); Padgett, J. F. and Ansell, C. K. (1993),
"Robust action and the rise of the Medici, 1400–1434", *American Journal of
Sociology* 98(6), 1259–1319. Distributed with statnet's `ergm` package
(GPL-3); the underlying data are historical facts collected by Padgett and
are in wide public use.

## Sampson monastery (liking)

- `sampson_vertices.tsv` — 18 monks with `name`, `group`
  (Loyal / Outcasts / Turks, the three-way faction assignment used by
  statnet's `samplike`), `cloisterville` (attended the minor seminary of
  "Cloisterville" before entering the monastery).
- `sampson_edges.tsv` — 88 directed "liking" ties: i → j if monk i named
  monk j among his top-three liked monks at *any* of the three time points
  (statnet's aggregated `samplike` network). `nominations` counts the number
  of time points (1–3) at which the tie was named.

Source: Sampson, S. F. (1968), *A novitiate in a period of change: An
experimental and case study of relationships*, PhD thesis, Cornell
University. Distributed with statnet's `ergm` package (GPL-3).

## faux.mesa.high (statnet `ergm`)

- `faux_mesa_high_vertices.tsv` — 205 students: `id`, `Grade` (integer,
  7–12), `Race` (`Black`, `Hisp`, `NatAm`, `Other`, `White`), `Sex` (`F`,
  `M`). Attribute names are statnet's, capitalised, so `nodematch("Grade")`
  reads the same in both languages.
- `faux_mesa_high_edges.tsv` — 203 undirected friendship ties (`source`,
  `target`, smaller id first, sorted). 163 join students of the same grade —
  the `nodematch.Grade` target that the provenanced
  `ERGMEgo.jl/test/fixtures/fauxmesa_ego_census.toml` pins to `[203, 163]`.

Source: `https://raw.githubusercontent.com/statnet/ergm/master/data/faux.mesa.high.RData`
(statnet `ergm` 4.x data file; the object is unchanged since ergm 2.x),
SHA-256 `7057b11cd8f1731ddef3e98c133e29e97c811037142b3543867e3ed6b7474c48`.
Written files: vertices
`0bf4369c17da449e0ffc689094656656081d4c6267642290050358f39dd4b8f6`, edges
`ccf73d580b5d6ffbbda9b3e565e3dbb9547313583b6fecc0a2bed2d286bef839`.
Licence: GPL-3 (ergm). The network is *simulated*: an ERGM fitted to an
Add Health school ("Mesa High", a pseudonym) and then drawn from, so it
carries no real student's data.

Citations: Resnick, M. D. et al. (1997), "Protecting adolescents from harm:
Findings from the National Longitudinal Study on Adolescent Health", *JAMA*
278(10), 823–832 (the Add Health source); Hunter, D. R., Goodreau, S. M. and
Handcock, M. S. (2008), "Goodness of fit of social network models", *JASA*
103(481), 248–258 (the fitted model the network was simulated from); Hunter,
D. R. et al. (2008), "ergm: A package to fit, simulate and diagnose
exponential-family models for networks", *JSS* 24(3) (the distributing
package).

## WTC police radio calls (Butts, Petrescu-Prahova & Cross 2007)

- `wtc_police_calls_events.tsv` — 481 radio calls among the 37 named
  communicants of a police unit at the World Trade Center on the morning of
  2001-09-11: `number` (1…481, the observed order — the data are *ordinal*,
  there is no clock), `source`, `recipient` (actor ids 1–37; no self-calls;
  35 of the 37 actors ever call or are called).
- `wtc_police_calls_actors.tsv` — the 37 actors: `id`, `is_icr` (1 for the
  three actors, 13, 14 and 32, who occupied *institutionalised coordinator
  roles*; R's `WTCPoliceIsICR`).

Source: the CRAN `relevent` package (1.2-1) ships **no** data directory; the
objects `WTCPoliceCalls` and `WTCPoliceIsICR` are distributed with statnet's
relevent workshop, downloaded from
`https://raw.githubusercontent.com/statnet/Workshops/master/relevent/relevent_workshop.Rdata`,
SHA-256 `92055b15d49584327046e7de6f30d5128d930c6703053570f1867fa38e5f0e15`
(the same file is served as `docs/relevent_tutorial_files/relevent_workshop.Rdata`
and used by the tutorial that accompanies `relevent`). Written files: events
`271a2753f0f071da5610890f2995f5c746f20b42e366be04c7d5c05c3dd95553`, actors
`6a903d560167577922178d265fc73984f563c4b4af75e3e9a4a8be9f62f32189`.
Licence: GPL-2+ (relevent and its workshop materials; author C. T. Butts).

Citation: Butts, C. T., Petrescu-Prahova, M. and Cross, B. R. (2007),
"Responder communication networks in the World Trade Center disaster:
Implications for modeling of communication within emergency settings",
*Journal of Mathematical Sociology* 31(2), 121–147. Model reference: Butts,
C. T. (2008), "A relational event framework for social action",
*Sociological Methodology* 38(1), 155–200.

## s50 (RSiena; Teenage Friends and Lifestyle Study)

Files under `s50/`, all header-less integer matrices with one row per actor:

- `friendship_wave1.tsv`, `friendship_wave2.tsv`, `friendship_wave3.tsv` —
  50 × 50 directed adjacency matrices (RSiena `s501`, `s502`, `s503`),
  entry (i, j) = 1 if girl i named girl j as a friend; 113, 116 and 122
  ties respectively, no self-ties.
- `alcohol.tsv` — 50 × 3 (RSiena `s50a`), alcohol use per wave coded
  1 = none, 2 = once or twice a year, 3 = once a month, 4 = once a week,
  5 = more than once a week.
- `smoke.tsv` — 50 × 3 (RSiena `s50s`), smoking per wave coded 1 = none,
  2 = occasional, 3 = regular (more than once a week).

The 50 girls are an excerpt of the 160-pupil Glasgow cohort of the Teenage
Friends and Lifestyle Study, observed at three yearly waves from 1995 (age
13). Source: the `s50` objects bundled with RSiena (GPL-2+). The committed
files were converted from the comma-separated copies Siena.jl carried under
`test/data/` — the same matrices from which the provenanced
`Siena.jl/test/fixtures/s50_siena07.toml` was generated with **RSiena
1.6.6** (its script `test/fixtures/r/s50_siena07.R` reads `s501`… from
RSiena directly). Input hashes (`s501.csv`, `s502.csv`, `s503.csv`,
`s50a.csv`, `s50s.csv`):
`aaf6cea3ace631bedeb750497bbfe85b65630a7ae1368480feaf0dd6cf713a08`,
`61dc36d2248f8bcf0a69eb01b4760db9079498b81fec40524608fb34cf879686`,
`c8abff615a2d3f0fa329be765a3fb80967d4266ee986f91b2de739aea9436c22`,
`f258e8e26bbda2999440c91925f63956c834ac7fec3cba8bfc8405add0c95c41`,
`0c586e028c9a3da6a4c143002fb780633989d61c01dd2ca36ceb159797f5527a`.
Written files: waves
`e56c16c1ae6f7d4df9e77c58a137aae1340582c1fd465b304d32df2ff7e40725`,
`416c4ac0e0802876ff886b04b4775383cc56902e4283486ac785fd588d7ddc29`,
`052828181b9d1034a5c78644b4d91930f0419932a32e57a80ccdb75583f777d7`; alcohol
`c587ae719afadc1d15744f5d90b8b64abaa16076a4421cf206f26f70a6d00928`; smoke
`1063cc8adecc088912ccddef768270a83aaf02a0653eb1de4a1915269480a53f`.
When RSiena is installed the extraction script reads the objects from it
instead (and the hashes of the written files must come out identical).

Citations: West, P. and Sweeting, H. (1996), *Background rationale and design
of the West of Scotland 11–16 Study*, MRC Medical Sociology Unit Working
Paper 52, Glasgow; Michell, L. and Amos, A. (1997), "Girls, pecking order and
smoking", *Social Science & Medicine* 44(12), 1861–1869; Ripley, R. M.,
Snijders, T. A. B., Boda, Z., Vörös, A. and Preciado, P. (2026), *Manual for
RSiena*, University of Oxford / University of Groningen (the distributing
package, section "s50 data").
