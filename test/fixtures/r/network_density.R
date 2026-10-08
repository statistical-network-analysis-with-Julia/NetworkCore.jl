# Golden fixture: R `network::network.density` on small hand-built networks, so
# that NetworkCore.jl's `network_density` (and every package that inherits it —
# SNA.gden, TSNA.t_density, ERGM's MCMC starting density) is pinned to the
# reference implementation rather than to a formula somebody remembered.
#
# Covers the whole keyword surface of the Julia routine:
#   * directed / undirected  ×  loops allowed or not      (R keys on has.loops)
#   * two-mode directed / undirected  ×  discount.bipartite FALSE / TRUE
#     (Julia `discount_bipartite=`; R's default FALSE counts within-mode dyads)
#   * one NA edge with na.omit = FALSE  (Julia `missing=:face`: the masked tie
#     is counted at its stored face value over the full dyad count). R's
#     default na.omit = TRUE is deliberately NOT reproduced — it reads an
#     unobserved tie as absent, which the ecosystem's missing-data contract
#     forbids — so the Julia side offers no such policy.
#   * a one-vertex network (NaN)
#   * `diag=FALSE` on the loops-allowed networks, from `sna::gden(diag=FALSE)`
#     (`network::network.density` has no `diag` argument): self-loops leave
#     the numerator AND the denominator. Julia `diag=false`. The loops networks
#     are chosen so that the default and `diag=false` values DIFFER (a 4-arc
#     directed example gives 4/16 == 3/12, which could not tell a silent no-op
#     from the real thing). `sna::gden(diag=TRUE)` is NOT recorded: on an
#     undirected graph it divides the symmetric matrix sum by n², which is not
#     `network.density`'s n(n+1)/2, and Julia's `diag=true` follows `network`.
#
# Needs `network` and `sna`. Regenerate from the package root:
#
#   Rscript test/fixtures/r/network_density.R > test/fixtures/network_density.toml

suppressMessages(library(network))
suppressMessages(library(sna))

seed <- 20260909
set.seed(seed)   # nothing here is stochastic; recorded so the provenance block is complete

fmt <- function(x) sprintf("%.17g", x)

# --- one-mode, no loops -----------------------------------------------------
dir_noloops <- network.initialize(4, directed = TRUE, loops = FALSE)
add.edges(dir_noloops, c(1, 2, 3), c(2, 3, 4))            # 3 arcs / 12

undir_noloops <- network.initialize(4, directed = FALSE, loops = FALSE)
add.edges(undir_noloops, c(1, 2, 3), c(2, 3, 4))          # 3 edges / 6

# --- one-mode, loops allowed (one self-loop present) -------------------------
dir_loops <- network.initialize(4, directed = TRUE, loops = TRUE)
add.edges(dir_loops, c(1, 1, 2, 3), c(2, 1, 3, 1))        # 4 arcs incl. 1->1 / 16

undir_loops <- network.initialize(4, directed = FALSE, loops = TRUE)
add.edges(undir_loops, c(1, 2, 3), c(2, 3, 3))            # 3 edges incl. 3-3 / 10

# --- loops allowed, chosen so diag=FALSE differs from the default ---------------
# directed: 5 arcs incl. 1->1: 5/16 by default, 4/12 with diag=FALSE
dir_loops_diag <- network.initialize(4, directed = TRUE, loops = TRUE)
add.edges(dir_loops_diag, c(1, 1, 2, 3, 4), c(2, 1, 3, 1, 2))
# undirected: 2 edges incl. 2-2: 2/10 by default, 1/6 with diag=FALSE
undir_loops_diag <- network.initialize(4, directed = FALSE, loops = TRUE)
add.edges(undir_loops_diag, c(2, 1), c(2, 2))

# --- two-mode (bipartite = 2 mode-1 vertices out of 5) -----------------------
bip_undir <- network.initialize(5, directed = FALSE, bipartite = 2)
add.edges(bip_undir, c(1, 1, 2), c(3, 4, 5))              # 3 edges

bip_dir <- network.initialize(5, directed = TRUE, bipartite = 2)
add.edges(bip_dir, c(1, 1, 4), c(3, 4, 2))                # 3 arcs, one mode2 -> mode1

# --- one NA edge -------------------------------------------------------------
na_net <- network.initialize(4, directed = TRUE, loops = FALSE)
add.edges(na_net, c(1, 2, 3), c(2, 3, 4))
set.edge.attribute(na_net, "na", TRUE, e = 2)             # 2->3 is unobserved (face: present)

one <- network.initialize(1, directed = TRUE)

cat('name = "network_density"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf('sna_version = "%s"\n', as.character(packageVersion("sna"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/network_density.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('dataset = "hand-built 4- and 5-vertex networks (see script)"\n\n')

cat("[tolerance]\n")
cat("# Deterministic ratio of two integers: agreement to floating-point rounding.\n")
cat("default = 1e-12\n\n")

cat("[values]\n")
cat("# one-mode: edge count / possible dyads, has.loops decides the diagonal\n")
cat(sprintf("directed_noloops = %s\n",   fmt(network.density(dir_noloops))))
cat(sprintf("undirected_noloops = %s\n", fmt(network.density(undir_noloops))))
cat(sprintf("directed_loops = %s\n",     fmt(network.density(dir_loops))))
cat(sprintf("undirected_loops = %s\n",   fmt(network.density(undir_loops))))
cat("# two-mode: discount.bipartite=FALSE (R default) counts within-mode dyads\n")
cat(sprintf("bipartite_undirected = %s\n",          fmt(network.density(bip_undir))))
cat(sprintf("bipartite_undirected_discount = %s\n", fmt(network.density(bip_undir, discount.bipartite = TRUE))))
cat(sprintf("bipartite_directed = %s\n",            fmt(network.density(bip_dir))))
cat(sprintf("bipartite_directed_discount = %s\n",   fmt(network.density(bip_dir, discount.bipartite = TRUE))))
cat("# one NA edge, na.omit=FALSE == Julia missing=:face (face value, full denominator)\n")
cat(sprintf("na_face = %s\n", fmt(network.density(na_net, na.omit = FALSE))))
cat("# diag: network.density keys on has.loops (default == Julia diag=true);\n")
cat("# sna::gden(diag=FALSE) drops self-loops from numerator and denominator (Julia diag=false)\n")
cat(sprintf("directed_loops_diag = %s\n",   fmt(network.density(dir_loops_diag))))
cat(sprintf("directed_loops_nodiag = %s\n", fmt(gden(dir_loops_diag, mode = "digraph", diag = FALSE))))
cat(sprintf("undirected_loops_diag = %s\n",   fmt(network.density(undir_loops_diag))))
cat(sprintf("undirected_loops_nodiag = %s\n", fmt(gden(undir_loops_diag, mode = "graph", diag = FALSE))))
cat(sprintf("undirected_loops_fixture_nodiag = %s\n", fmt(gden(undir_loops, mode = "graph", diag = FALSE))))
cat("# a single vertex has no dyads\n")
cat(sprintf("single_vertex_is_nan = %s\n", if (is.nan(network.density(one))) "true" else "false"))
