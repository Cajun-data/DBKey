# Usage: Rscript --vanilla run_dbkey_local.R input.sptxt output.db --mass-analyzer IT|FT [--fixed-carbamidomethyl]
# By default, only modifications explicitly listed in Mods= alter peptide masses.
# The optional flag adds implicit carbamidomethyl to every C for legacy files;
# do not combine that flag with explicit carbamidomethyl Mods= entries.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript --vanilla run_dbkey_local.R input.sptxt output.db --mass-analyzer IT|FT [--fixed-carbamidomethyl]"
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(status = 0L)
}
if (length(args) < 2L) stop(usage)
options <- args[-c(1L, 2L)]
mass_analyzer <- NULL
fixed_carbamidomethyl <- FALSE
while (length(options)) {
  if (options[1L] == "--mass-analyzer") {
    if (length(options) < 2L || !is.null(mass_analyzer))
      stop("Supply --mass-analyzer exactly once with IT or FT")
    mass_analyzer <- toupper(options[2L])
    options <- options[-c(1L, 2L)]
  } else if (options[1L] == "--fixed-carbamidomethyl") {
    fixed_carbamidomethyl <- TRUE
    options <- options[-1L]
  } else {
    stop("Unknown option: ", options[1L], "\n", usage)
  }
}
if (is.null(mass_analyzer) || !mass_analyzer %in% c("IT", "FT"))
  stop("Specify --mass-analyzer IT or FT to match the fragment mass analyzer")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
root <- dirname(dirname(script_path))
Sys.setenv(DBKEY_HOME = root, OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
Sys.setenv(DBKEY_FIXED_CARBAMIDOMETHYL = if (fixed_carbamidomethyl) "TRUE" else "FALSE")
input <- normalizePath(args[1L], mustWork = TRUE)
output <- normalizePath(args[2L], mustWork = FALSE)
if (file.exists(output)) stop("Output already exists; choose a new output path")
if (!dir.exists(dirname(output))) stop("Output directory does not exist")
if (!grepl("\\.sptxt$", input, ignore.case = TRUE)) stop("This runner accepts SpectraST .sptxt files")
source(file.path(root, "lib", "msp2db.R"))
temporary_output <- tempfile(pattern = paste0(basename(output), ".building-"), tmpdir = dirname(output), fileext = ".db")
message("Input: ", input)
message("Output in progress: ", temporary_output)
DBbuilder(
  Library = data.frame(name = basename(input), datapath = input),
  FragmentationMode = "Read From file", MassAnalyzer = mass_analyzer,
  CollisionEnergy = "Read from file", CompoundClass = NULL,
  Filter = FALSE, DBoutput = temporary_output, topX = 0,
  cutoff = 0, TMTPro = FALSE, massOffset = NULL, IonTypes = NULL,
  deltaFragment = FALSE, oldMod = 0, newMod = 0)
db <- DBI::dbConnect(RSQLite::SQLite(), temporary_output)
check <- DBI::dbGetQuery(db, "PRAGMA integrity_check")[[1L]]
counts <- DBI::dbGetQuery(db, "SELECT (SELECT COUNT(*) FROM CompoundTable) AS compounds, (SELECT COUNT(*) FROM SpectrumTable) AS spectra")
DBI::dbDisconnect(db)
if (!identical(check, "ok")) stop("SQLite integrity check failed")
if (counts$compounds != counts$spectra || counts$spectra == 0) stop("Invalid output spectrum count")
if (!file.rename(temporary_output, output)) stop("Could not publish completed database: ", temporary_output)
message("Created ", output, " with ", counts$spectra, " spectra; SQLite integrity_check=ok")
