[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$InputLibrary,
    [Parameter(Mandatory=$true)][ValidateSet('IT', 'FT')][string]$MassAnalyzer,
    [string]$OutputDatabase,
    [switch]$FixedCarbamidomethyl,
    [string]$RscriptPath = 'Rscript'
)
$ErrorActionPreference = 'Stop'
$rCommand = Get-Command $RscriptPath -ErrorAction Stop
$resolvedInput = (Resolve-Path -LiteralPath $InputLibrary).Path
if (-not $OutputDatabase) {
    $OutputDatabase = [System.IO.Path]::ChangeExtension($resolvedInput, '.db')
}
$rArguments = @('--vanilla', (Join-Path $PSScriptRoot 'run_dbkey_local.R'), $resolvedInput, $OutputDatabase, '--mass-analyzer', $MassAnalyzer)
if ($FixedCarbamidomethyl) { $rArguments += '--fixed-carbamidomethyl' }
& $rCommand.Source @rArguments
if ($LASTEXITCODE -ne 0) {
    throw "DBKey conversion failed (exit code $LASTEXITCODE). See the R error above."
}
