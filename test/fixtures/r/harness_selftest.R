# Golden-fixture harness self-test: freeze R `sna` output on the Padgett
# Florentine marriage network in the canonical TOML layout every fixture in this
# ecosystem uses.
#
# This is the reference EXAMPLE for writing a fixture generator. Emit three
# blocks:
#   [provenance] — enough to regenerate this file from scratch. Without it a
#                  frozen number is folklore: when it disagrees with the Julia
#                  code, nobody can tell whether the code broke or R moved.
#   [tolerance]  — what agreement is deemed acceptable, and why. For a Monte
#                  Carlo estimator this is a scientific claim, so it belongs
#                  somewhere a reviewer can see it — not buried at a call site.
#   [values]     — the frozen output.
#
# Regenerate from the package root:
#
#   Rscript test/fixtures/r/harness_selftest.R > test/fixtures/florentine_sna.toml

suppressMessages({
  library(network)
  library(sna)
})

seed <- 20260713
set.seed(seed)

data(flo, package = "network")
g <- network(flo, directed = FALSE)

deg  <- sna::degree(g, gmode = "graph")
dens <- sna::gden(g, mode = "graph")
tri  <- sna::gtrans(g, mode = "graph", measure = "weak")

cat('name = "florentine_sna"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('sna_version = "%s"\n', as.character(packageVersion("sna"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/harness_selftest.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('dataset = "network::flo (Padgett Florentine marriage)"\n\n')

cat("[tolerance]\n")
cat("# Deterministic descriptive statistics: exact agreement is required.\n")
cat("default = 1e-10\n\n")

cat("[values]\n")
cat(sprintf("degree = [%s]\n", paste(deg, collapse = ", ")))
cat(sprintf("density = %.17g\n", dens))
cat(sprintf("transitivity_weak = %.17g\n", tri))
