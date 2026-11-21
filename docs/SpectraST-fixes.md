# SpectraST conversion fixes

These changes build on SchweppeLab/DBKey commit
`658dbec26b861ce6b0940f93c2a25209fa0ed078`.
The changes focus on SpectraST conversion. The Prosit and Skyline parsers remain
unchanged and have not been comprehensively validated with this update.

## Corrected behavior

- Equal-length spectrum peak lists stay as lists, including single-spectrum input.
  Missing fragment mass-error annotations no longer cause column-index failures.
- Filter OFF retains finite m/z values with finite, positive intensities, sorted
  by m/z. It bypasses Top N, intensity cutoff and ion-type filtering. Filter ON
  selects the strongest Top N peaks, with input order breaking intensity ties.
- Modification names use exact lookup in the existing custom UniMod table.
  Numeric shifts retain their signs; unknown names fail with a clear message.
  The existing `iTRAQ8plex` mapping to 304.207146 Da is unchanged.
- Cysteine is unmodified unless the source `Mods=` explicitly specifies a
  modification. This avoids silently adding carbamidomethyl to unmodified C or
  counting an explicit carbamidomethyl modification twice.
- The `+18` annotation exclusion is literal, preserving valid b18/y18 labels.
- Complete records are read in 100,000-line blocks without a PSOCK worker pool.
  A record crossing a block boundary is carried into the next block. Parsed
  output is still accumulated in memory before SQLite writing; this is not a
  fully constant-memory database writer.
- `DBKEY_HOME` selects the source/data directory; its default remains
  `~/Repos/MSPtoDB` for the existing Docker layout. SQLite is disconnected once.
- The version-5 SQLite schema and double-precision peak blobs are unchanged.

## Native conversion

Install the R dependencies listed in `dependencies/Dockerfile`, including
BiocParallel. The native runner determines the repository root from its own
location, so it can be launched from any working directory:

```sh
Rscript --vanilla scripts/run_dbkey_local.R input.sptxt output.db --mass-analyzer FT
```

Choose `IT` or `FT` to match the fragment mass analyzer; there is no implicit
analyzer default. The runner reads fragmentation mode and collision energy
from each input file, retains all valid peaks, builds a temporary database, checks
SQLite integrity and spectrum counts, then renames the completed output.
It refuses an existing output path. A failed run can leave a `.building-*.db`
file for diagnosis.

PowerShell users can run:

```powershell
./scripts/Run-DBKey-Fixed.ps1 -InputLibrary input.sptxt -OutputDatabase output.db -MassAnalyzer FT
```

Use `-RscriptPath 'C:\path\to\Rscript.exe'` if Rscript is not on PATH.

For legacy files that deliberately omit a known fixed cysteine modification:

```sh
Rscript --vanilla scripts/run_dbkey_local.R input.sptxt output.db --mass-analyzer FT --fixed-carbamidomethyl
```

The equivalent setting for direct parser or Docker use is
`DBKEY_FIXED_CARBAMIDOMETHYL=TRUE`. Do not enable it when carbamidomethyl is
already present in `Mods=`, because it adds a fixed modification to every C
in addition to explicit modifications.

## Docker

Build the fork locally to include these source changes; an existing upstream
image does not contain them:

```sh
docker build -t dbkey-fixed .
docker run --rm -p 3838:3838 dbkey-fixed
```

The optional `scripts/Install-DBKey-Fix.ps1 -Container dbkey` copies the two
patched R sources into an existing container after backing up its sources and
checking R syntax. It restarts that container, so run it only when no conversion
is active. The installer was syntax-checked but not executed against Docker
during validation.

## Regression checks

```sh
Rscript --vanilla tests/run_regression_tests.R
```

The tests generate synthetic spectra and temporary databases. They check
single/equal-size peak lists, strongest-peak selection and deterministic ties,
Filter OFF, missing annotations, b18/y18 retention, named/signed modifications,
explicit and legacy cysteine handling, a 100,000-line read boundary, SQLite
schema/integrity/metadata, exact peak blobs, explicit IT/FT analyzer selection,
invalid or missing analyzer arguments, and refusal to overwrite an output.

Validation covers native R conversion with synthetic inputs. Docker execution,
mzVault opening and instrument acquisition compatibility have not been tested.
