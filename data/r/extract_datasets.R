# Extract the bundled teaching datasets from their upstream R objects into the
# tab-separated files that `NetworkCore.load_dataset` reads. Needs only base R
# (>= 4.5, for tools::sha256sum) and the CRAN `network` package; no other
# statnet package is required, so the extraction runs on a machine that cannot
# build `ergm`, `relevent` or `RSiena`.
#
# Regenerate from the package root (network access required for the two
# downloads):
#
#   Rscript data/r/extract_datasets.R
#
# The script prints the SHA-256 of every upstream file it consumed; those
# hashes are recorded in data/README.md so a regeneration can be checked
# against the copy that produced the committed TSVs.
#
# Three sources, three datasets:
#
#   1. faux.mesa.high  -- statnet/ergm's data/faux.mesa.high.RData (GPL-3):
#      205 students, 203 undirected friendship ties, vertex attributes
#      Grade / Race / Sex. Written as faux_mesa_high_vertices.tsv and
#      faux_mesa_high_edges.tsv.
#
#   2. WTCPoliceCalls  -- the World Trade Center police radio-call event
#      stream (Butts, Petrescu-Prahova & Cross 2007), 481 ordinal events among
#      37 responders, plus the WTCPoliceIsICR indicator of the three
#      institutionalised coordinators. The CRAN `relevent` tarball (1.2-1)
#      ships NO data directory; the objects are distributed with statnet's
#      relevent workshop, https://github.com/statnet/Workshops
#      (relevent/relevent_workshop.Rdata), which is what is downloaded here.
#      Written as wtc_police_calls_events.tsv and wtc_police_calls_actors.tsv.
#
#   3. s50            -- the RSiena "s50" excerpt of the Teenage Friends and
#      Lifestyle Study (West & Sweeting 1996; Michell & Amos 1997): three
#      50 x 50 friendship waves (s501, s502, s503), alcohol (s50a, 50 x 3) and
#      smoking (s50s, 50 x 3). Taken from `RSiena` when it is installed;
#      otherwise from the comma-separated copies Siena.jl carried under
#      test/data/ (the same matrices the provenanced
#      Siena.jl/test/fixtures/s50_siena07.toml was generated from with RSiena
#      1.6.6). Written under s50/ as friendship_wave{1,2,3}.tsv, alcohol.tsv
#      and smoke.tsv.

suppressPackageStartupMessages(library(network))

script_dir <- local({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
    if (length(file_arg) == 0) getwd() else dirname(normalizePath(file_arg[1]))
})
data_dir <- normalizePath(file.path(script_dir, ".."))
snwj_root <- normalizePath(file.path(data_dir, "..", ".."), mustWork = FALSE)

write_tsv <- function(x, file) {
    path <- file.path(data_dir, file)
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !is.null(colnames(x)))
    cat(sprintf("  wrote %-40s %s\n", file, tools::sha256sum(path)))
}

fetch <- function(url, dest) {
    cat(sprintf("downloading %s\n", url))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    cat(sprintf("  sha256 %s  %s\n", tools::sha256sum(dest), basename(dest)))
    dest
}

tmp <- tempdir()

## ---------------------------------------------------------------- faux.mesa.high
cat("== faux.mesa.high (statnet/ergm)\n")
fmh_file <- fetch(
    "https://raw.githubusercontent.com/statnet/ergm/master/data/faux.mesa.high.RData",
    file.path(tmp, "faux.mesa.high.RData"))
fmh_env <- new.env()
load(fmh_file, envir = fmh_env)
fmh <- get("faux.mesa.high", envir = fmh_env)
stopifnot(inherits(fmh, "network"),
          network.size(fmh) == 205L,
          network.edgecount(fmh) == 203L,
          !is.directed(fmh))

grade <- as.integer(get.vertex.attribute(fmh, "Grade"))
stopifnot(all(grade %in% 7:12))
vertices <- data.frame(id = seq_len(network.size(fmh)),
                       Grade = grade,
                       Race = get.vertex.attribute(fmh, "Race"),
                       Sex = get.vertex.attribute(fmh, "Sex"),
                       stringsAsFactors = FALSE)
# R's edgelist for an undirected network lists each edge once with the
# smaller endpoint first; sort for a stable file.
el <- as.matrix(fmh, matrix.type = "edgelist")
el <- t(apply(el, 1, sort))
el <- el[order(el[, 1], el[, 2]), , drop = FALSE]
stopifnot(nrow(el) == 203L, !any(el[, 1] == el[, 2]))
edges <- data.frame(source = el[, 1], target = el[, 2])
write_tsv(vertices, "faux_mesa_high_vertices.tsv")
write_tsv(edges, "faux_mesa_high_edges.tsv")
cat(sprintf("  same-Grade edges: %d (ERGMEgo fixture nodematch.Grade target)\n",
            sum(grade[el[, 1]] == grade[el[, 2]])))

## ---------------------------------------------------------------- WTCPoliceCalls
cat("== WTCPoliceCalls (statnet relevent workshop)\n")
wtc_file <- fetch(
    "https://raw.githubusercontent.com/statnet/Workshops/master/relevent/relevent_workshop.Rdata",
    file.path(tmp, "relevent_workshop.Rdata"))
wtc_env <- new.env()
load(wtc_file, envir = wtc_env)
calls <- get("WTCPoliceCalls", envir = wtc_env)
is_icr <- get("WTCPoliceIsICR", envir = wtc_env)
stopifnot(is.data.frame(calls),
          identical(colnames(calls), c("number", "source", "recipient")),
          nrow(calls) == 481L,
          all(calls$number == seq_len(nrow(calls))),
          all(calls$source %in% 1:37), all(calls$recipient %in% 1:37),
          is.logical(is_icr), length(is_icr) == 37L)
events <- data.frame(number = as.integer(calls$number),
                     source = as.integer(calls$source),
                     recipient = as.integer(calls$recipient))
actors <- data.frame(id = seq_along(is_icr), is_icr = as.integer(is_icr))
write_tsv(events, "wtc_police_calls_events.tsv")
write_tsv(actors, "wtc_police_calls_actors.tsv")
cat(sprintf("  ICR actors: %s\n", paste(which(is_icr), collapse = ", ")))

## ---------------------------------------------------------------- s50 (RSiena)
cat("== s50 (RSiena)\n")
dir.create(file.path(data_dir, "s50"), showWarnings = FALSE)
s50 <- NULL
if (requireNamespace("RSiena", quietly = TRUE)) {
    cat(sprintf("  from RSiena %s\n", as.character(utils::packageVersion("RSiena"))))
    s50 <- list(friendship = list(RSiena::s501, RSiena::s502, RSiena::s503),
                alcohol = RSiena::s50a, smoke = RSiena::s50s)
} else {
    csv_dir <- file.path(snwj_root, "Siena.jl", "test", "data")
    if (!dir.exists(csv_dir))
        stop("RSiena is not installed and ", csv_dir, " does not exist: ",
             "install RSiena (the s50 objects are RSiena::s501/s502/s503/s50a/s50s) ",
             "or check out Siena.jl beside NetworkCore.jl")
    cat(sprintf("  RSiena not installed; reading Siena.jl's CSV copies from %s\n", csv_dir))
    read_csv_matrix <- function(f) {
        path <- file.path(csv_dir, f)
        cat(sprintf("  sha256 %s  %s\n", tools::sha256sum(path), f))
        as.matrix(utils::read.csv(path, header = FALSE))
    }
    s50 <- list(friendship = list(read_csv_matrix("s501.csv"),
                                  read_csv_matrix("s502.csv"),
                                  read_csv_matrix("s503.csv")),
                alcohol = read_csv_matrix("s50a.csv"),
                smoke = read_csv_matrix("s50s.csv"))
}
for (w in 1:3) {
    m <- s50$friendship[[w]]
    m <- matrix(as.integer(m), nrow(m), ncol(m))
    stopifnot(dim(m) == c(50L, 50L), all(m %in% 0:1), all(diag(m) == 0L))
    write_tsv(m, sprintf("s50/friendship_wave%d.tsv", w))
    cat(sprintf("  wave %d ties: %d\n", w, sum(m)))
}
alc <- matrix(as.integer(s50$alcohol), 50L, 3L)
smo <- matrix(as.integer(s50$smoke), 50L, 3L)
stopifnot(all(alc %in% 1:5), all(smo %in% 1:3))
write_tsv(alc, "s50/alcohol.tsv")
write_tsv(smo, "s50/smoke.tsv")
cat("done\n")
