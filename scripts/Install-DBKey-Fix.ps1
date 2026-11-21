# Apply the reviewed R source fixes to an existing DBKey container.
# This restarts only the named container. Run when no DBKey conversion is active.
[CmdletBinding()]
param([string]$Container = 'dbkey')
$ErrorActionPreference = 'Stop'
$appRoot = '/root/Repos/MSPtoDB'
$files = @('msp2db.R', 'SpXLibraryParser.R')
$patchedRoot = Join-Path (Split-Path -Parent $PSScriptRoot) 'lib'
foreach ($name in $files) {
    if (-not (Test-Path -LiteralPath (Join-Path $patchedRoot $name) -PathType Leaf)) {
        throw "Missing patched source: $name"
    }
}
& docker inspect --format '{{.State.Running}}' $Container
if ($LASTEXITCODE -ne 0) { throw 'Cannot access the requested container.' }
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$backup = Join-Path $PSScriptRoot "container_backup_$stamp"
New-Item -ItemType Directory -Path $backup | Out-Null
foreach ($name in $files) {
    $remoteOriginal = '{0}:{1}/lib/{2}' -f $Container, $appRoot, $name
    & docker cp $remoteOriginal (Join-Path $backup $name)
    if ($LASTEXITCODE -ne 0) { throw "Backup failed for $name; no source files have been replaced." }
}
$stage = '/tmp/dbkey_fix_' + [Guid]::NewGuid().ToString('N')
& docker exec $Container mkdir -p $stage
if ($LASTEXITCODE -ne 0) { throw 'Could not create the staging directory.' }
foreach ($name in $files) {
    & docker cp (Join-Path $patchedRoot $name) ('{0}:{1}/{2}' -f $Container, $stage, $name)
    if ($LASTEXITCODE -ne 0) { throw "Staging failed for $name." }
}
& docker exec $Container Rscript --vanilla -e 'invisible(lapply(commandArgs(TRUE), function(x) parse(file=x)))' "$stage/msp2db.R" "$stage/SpXLibraryParser.R"
if ($LASTEXITCODE -ne 0) { throw 'R syntax check failed; original sources were retained.' }
foreach ($name in $files) {
    & docker cp (Join-Path $patchedRoot $name) ('{0}:{1}/lib/{2}' -f $Container, $appRoot, $name)
    if ($LASTEXITCODE -ne 0) { throw "Installing $name failed. Originals are saved at $backup." }
}
& docker restart $Container
if ($LASTEXITCODE -ne 0) { throw "Restart failed. Originals are saved at $backup." }
Write-Host "Installed DBKey fixes. Original R sources: $backup"
Write-Host 'Use Filter OFF to retain every positive finite peak, or Filter ON for strongest Top N peaks.'
Write-Host 'Cysteine masses now follow explicit Mods= values. Implicit carbamidomethyl is disabled by default.'
