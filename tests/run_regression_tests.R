# Run from any directory: Rscript --vanilla /path/to/DBKey/tests/run_regression_tests.R
# All spectra below are synthetic; no experimental libraries are needed.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)))
Sys.setenv(DBKEY_HOME = root, DBKEY_FIXED_CARBAMIDOMETHYL = "FALSE",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressPackageStartupMessages(source(file.path(root, "lib", "msp2db.R")))
suppressPackageStartupMessages(source(file.path(root, "lib", "SpXLibraryParser.R")))
data.table::setDTthreads(1L)

record <- function(sequence = "ACDEFGHIKLMNPQRSTVWY", mods = "Mods=0",
                   n = 100L, annotation = rep("y1/0.01", n),
                   intensity = (seq_len(n) * 37) %% 101,
                   masses = rev(200 + seq_len(n)), id = 0L) {
  name <- paste0(sequence, "/2")
  c(paste("Name:", name), paste("LibID:", id), "MW: 1000.0000",
    "PrecursorMZ: 501.0073", "Status: Normal", paste("FullName:", name, "(HCD)"),
    paste("Comment: CollisionEnergy=35.0", mods, "RetentionTime=123.4,123.4,123.4"),
    paste("NumPeaks:", n), paste(masses, intensity, annotation, sep = "\t"), "")
}
parse_library <- function(lines, filter = FALSE, top = 20L, cutoff = 0, ions = NULL) {
  SpXLibraryParser(lines, "Read From file", "IT", "Read from file", NULL,
                  Filter = filter, TMTPro = FALSE, Source = "SpectraST", topX = top,
                  cutoff = cutoff, massOffset = NULL, IonTypes = ions,
                  deltaFragment = FALSE, oldMod = 0, newMod = 0)
}
unpack <- function(x) readBin(x, numeric(), n = length(x) / 8L, size = 8L, endian = "little")
expect_error <- function(code, pattern) {
  error <- tryCatch({ force(code); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}

for (n in c(1L, 3L)) {
  input <- unlist(lapply(seq_len(n), function(i) record(id = i)))
  all <- parse_library(input)
  top100 <- parse_library(input, TRUE, 100L)
  top20 <- parse_library(input, TRUE, 20L)
  stopifnot(nrow(all) == n, nrow(top100) == n, nrow(top20) == n,
            sum(names(all) == "CompoundClass") == 1L)
  for (i in seq_len(n)) {
    stopifnot(length(unpack(all$blobInt[[i]])) == 100L,
              identical(all$blobMass[[i]], top100$blobMass[[i]]),
              identical(all$blobInt[[i]], top100$blobInt[[i]]),
              identical(sort(unpack(top20$blobInt[[i]])), as.double(81:100)),
              !is.unsorted(unpack(all$blobMass[[i]])))
  }
}
cat("PASS: one/multiple spectra with equal peak counts, peak blobs, strongest Top N\n")

baseline <- parse_library(record())
bypassed <- parse_library(record(), top = 0, cutoff = 100, ions = "b")
stopifnot(identical(baseline$blobMass, bypassed$blobMass),
          identical(baseline$blobInt, bypassed$blobInt))
for (top in c(0, -1, 1.5, NA_real_))
  expect_error(parse_library(record(), TRUE, top), "Top N peaks must be a positive integer")
ties <- parse_library(record(n = 4L, intensity = c(50, 50, 50, 10)), TRUE, 2L)
stopifnot(identical(unpack(ties$blobMass[[1L]]), c(203, 204)))
invalid <- parse_library(record(n = 5L, masses = c(Inf, 202, 203, 204, 205),
                                intensity = c(1, 0, -1, Inf, 50)))
stopifnot(identical(unpack(invalid$blobMass[[1L]]), 205),
          identical(unpack(invalid$blobInt[[1L]]), 50))
cat("PASS: Filter OFF bypass, Top N validation, deterministic ties, invalid peaks\n")

annotations <- parse_library(record(n = 5L, masses = 201:205,
                                    annotation = c("b18", "y18^2", "b2+18", "?", "y3/0.02")))
# DBKey represents unknown/excluded ions as empty fields after sanitization.
stopifnot(grepl("ions:b18;y18^2;;;y3", annotations$Tags, fixed = TRUE))
missing_errors <- parse_library(record(annotation = rep("?", 100L)))
stopifnot(nrow(missing_errors) == 1L, length(unpack(missing_errors$blobMass[[1L]])) == 100L)
long_header <- parse_library(c(rep("### synthetic header", 120L), record()))
stopifnot(identical(baseline$blobMass, long_header$blobMass))
expect_error(parse_library("### no spectra"), "No Name: entry found")
cat("PASS: b18/y18, literal +18, missing mass errors, long preambles\n")

unmodified <- parse_library(record(sequence = "AC"))
# Independent elemental composition of neutral AC: C6 H12 N2 O3 S1.
expected_ac <- (6 * 12 + 12 * 1.0078250321 + 2 * 14.0030740052 +
                  3 * 15.9949146221 + 31.97207069 + 2 * 1.007276466) / 2
stopifnot(abs(unmodified$PrecursorMasses - expected_ac) < 1e-8)
explicit_cam <- parse_library(record(sequence = "AC", mods = "Mods=1/1,C,Carbamidomethyl"))
stopifnot(abs(explicit_cam$PrecursorMasses - expected_ac - 57.021464 / 2) < 1e-8)
Sys.setenv(DBKEY_FIXED_CARBAMIDOMETHYL = "TRUE")
legacy <- parse_library(record(sequence = "AC"))
Sys.setenv(DBKEY_FIXED_CARBAMIDOMETHYL = "FALSE")
stopifnot(abs(legacy$PrecursorMasses - expected_ac - 57.021464 / 2) < 1e-5)
no_mod <- parse_library(record(sequence = "QC"))
named <- parse_library(record(sequence = "QC", mods = "Mods=1/0,Q,Gln->pyro-Glu"))
numeric <- parse_library(record(sequence = "QC", mods = "Mods=1/0,Q,-17.026549"))
positive <- parse_library(record(sequence = "QC", mods = "Mods=1/0,Q,+15.5"))
stopifnot(identical(named$PrecursorMasses, numeric$PrecursorMasses),
          abs(no_mod$PrecursorMasses - named$PrecursorMasses - 17.026549 / 2) < 1e-8,
          abs(positive$PrecursorMasses - no_mod$PrecursorMasses - 15.5 / 2) < 1e-8,
          identical(unname(named$Formula), "-17.026549@0"))
label <- parse_library(record(sequence = "AC", mods = "Mods=1/-1,A,iTRAQ8plex"))
stopifnot(abs(label$PrecursorMasses - expected_ac - 304.207146 / 2) < 1e-8)
expect_error(parse_library(record(mods = "Mods=1/0,A,UnrecognizedTestMod")),
             "Unknown modification(s): UnrecognizedTestMod")
cat("PASS: explicit/legacy cysteine, exact named and signed numeric modifications\n")

build_and_check <- function(input, expected_count) {
  source_path <- tempfile(fileext = ".sptxt")
  database_path <- tempfile(fileext = ".db")
  on.exit(unlink(c(source_path, database_path)), add = TRUE)
  writeLines(input, source_path)
  before_threads <- data.table::getDTthreads()
  result <- DBbuilder(
    Library = data.frame(name = basename(source_path), datapath = source_path),
    FragmentationMode = "Read From file", MassAnalyzer = "IT", CollisionEnergy = "Read from file",
    CompoundClass = NULL, Filter = FALSE, DBoutput = database_path,
    topX = 20, cutoff = 0, TMTPro = FALSE, massOffset = NULL, IonTypes = NULL,
    deltaFragment = FALSE, oldMod = 0, newMod = 0)
  stopifnot(identical(result, database_path), identical(data.table::getDTthreads(), before_threads))
  db <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  on.exit(DBI::dbDisconnect(db), add = TRUE, after = FALSE)
  stopifnot(identical(DBI::dbGetQuery(db, "PRAGMA integrity_check")[[1L]], "ok"),
            nrow(DBI::dbGetQuery(db, "PRAGMA foreign_key_check")) == 0L,
            DBI::dbGetQuery(db, "SELECT COUNT(*) AS n FROM CompoundTable")$n == expected_count,
            DBI::dbGetQuery(db, "SELECT COUNT(*) AS n FROM SpectrumTable")$n == expected_count,
            DBI::dbGetQuery(db, "SELECT COUNT(*) AS n FROM SpectrumTable WHERE length(blobMass)=800 AND length(blobIntensity)=800 AND MassAnalyzer='IT' AND FragmentationMode='HCD' AND CAST(CollisionEnergy AS REAL)=35 AND RetentionTime=123.4")$n == expected_count,
            DBI::dbGetQuery(db, "SELECT version FROM HeaderTable")$version == 5L)
  spectra <- DBI::dbGetQuery(db, "SELECT blobMass, blobIntensity FROM SpectrumTable ORDER BY SpectrumId")
  stopifnot(all(vapply(spectra$blobMass, identical, logical(1), baseline$blobMass[[1L]])),
            all(vapply(spectra$blobIntensity, identical, logical(1), baseline$blobInt[[1L]])))
}
build_and_check(record(), 1L)
build_and_check(rep(record(), 3L), 3L)
# 999 full records plus one final record span the actual 100,000-line read boundary.
build_and_check(rep(record(), 1000L), 1000L)
cat("PASS: SQLite schema, metadata, peak blobs, small inputs and 100,000-line boundary\n")

# Run the public entry point and verify its refusal to replace an existing output.
input <- tempfile(fileext = ".sptxt")
output <- tempfile(fileext = ".db")
writeLines(record(sequence = "AC"), input)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
runner <- file.path(root, "scripts", "run_dbkey_local.R")
invoke <- function(options = c("--mass-analyzer", "IT")) {
  suppressWarnings(system2(rscript,
    c("--vanilla", shQuote(runner), shQuote(input), shQuote(output), options),
    stdout = TRUE, stderr = TRUE))
}
for (options in list(character(), c("--mass-analyzer", "invalid"))) {
  rejected <- invoke(options)
  stopifnot(attr(rejected, "status") != 0L,
            any(grepl("Specify --mass-analyzer IT or FT", rejected, fixed = TRUE)),
            !file.exists(output))
}
for (analyzer in c("IT", "FT")) {
  run <- invoke(c("--mass-analyzer", analyzer))
  if (!is.null(attr(run, "status"))) stop(paste(run, collapse = "\n"))
  stopifnot(file.exists(output))
  db <- DBI::dbConnect(RSQLite::SQLite(), output)
  stored_analyzer <- DBI::dbGetQuery(db, "SELECT DISTINCT MassAnalyzer FROM SpectrumTable")$MassAnalyzer
  DBI::dbDisconnect(db)
  stopifnot(identical(stored_analyzer, analyzer))
  before <- unname(tools::md5sum(output))
  repeat_run <- invoke(c("--mass-analyzer", analyzer))
  stopifnot(attr(repeat_run, "status") != 0L,
            any(grepl("Output already exists", repeat_run, fixed = TRUE)),
            identical(before, unname(tools::md5sum(output))))
  unlink(output)
}
unlink(input)
cat("PASS: portable CLI, required IT/FT analyzer selection, and existing-output protection\nAll regression checks passed.\n")
