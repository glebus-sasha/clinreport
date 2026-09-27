[CmdletBinding()]
param(
    [string]$OutputDir = "results/dev-report",
    [int]$Port = 3838,
    [string]$BindHost = "127.0.0.1",
    [string]$PatientId = "PTA_27",
    [double]$PadjCutoff = 0.05,
    [double]$Log2fcCutoff = 1.0,
    [int]$MaxNeighbors = 40,
    [string]$RscriptPath,
    [switch]$PrepareOnly
)

$ErrorActionPreference = "Stop"
$projectDir = Split-Path -Parent $PSCommandPath
Set-Location $projectDir

if (-not $RscriptPath) {
    $rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rscriptCommand) {
        $RscriptPath = $rscriptCommand.Source
    } else {
        $rCandidates = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName "bin\Rscript.exe" } |
            Where-Object { Test-Path -LiteralPath $_ }
        if ($rCandidates) { $RscriptPath = @($rCandidates)[0] }
    }
}
if (-not $RscriptPath -or -not (Test-Path -LiteralPath $RscriptPath)) {
    throw "Rscript was not found. Supply -RscriptPath 'C:\Program Files\R\R-x.y.z\bin\Rscript.exe'."
}
if ($Port -lt 1 -or $Port -gt 65535) { throw "Port must be from 1 to 65535." }
if ($MaxNeighbors -lt 1) { throw "MaxNeighbors must be positive." }

$resultsDir = [IO.Path]::GetFullPath((Join-Path $projectDir "results"))
$outputPath = [IO.Path]::GetFullPath((Join-Path $projectDir $OutputDir))
if (-not $outputPath.StartsWith($resultsDir + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
    throw "OutputDir must be inside $resultsDir"
}

$rawDir = Join-Path $projectDir "raw"
$currentInputDir = Join-Path $rawDir "pta_27"
$drugDetailsDir = Join-Path $currentInputDir "drugs"
$gseaDir = Join-Path $currentInputDir "gsea"
$gmtListFile = Join-Path $gseaDir "gmt_list.csv"
if (-not (Test-Path -LiteralPath $gmtListFile)) { throw "GMT collection list is missing: $gmtListFile" }
$collections = @(Import-Csv -LiteralPath $gmtListFile)
if ($collections.Count -eq 0) { throw "GMT collection list is empty: $gmtListFile" }
$gmtFiles = @($collections | ForEach-Object { Join-Path $gseaDir ([IO.Path]::GetFileName(([uri]$_.url).AbsolutePath)) })
$gseaReports = @(Get-ChildItem -LiteralPath $gseaDir -File -Filter '*gsea_report_for_*.tsv' | Sort-Object Name)
if ($gmtFiles.Count -ne @(Get-ChildItem -LiteralPath $gseaDir -File -Filter '*.gmt').Count) {
    throw "GMT files do not match the collections in $gmtListFile"
}
foreach ($gmtFile in $gmtFiles) {
    if (-not (Test-Path -LiteralPath $gmtFile)) { throw "GMT input is missing: $gmtFile" }
    $stem = [IO.Path]::GetFileNameWithoutExtension($gmtFile)
    $matches = @($gseaReports | Where-Object { $_.Name.Contains("$stem.gsea_report_for_") })
    if ($matches.Count -eq 0) { throw "No GSEA reports match GMT collection $stem" }
}
foreach ($drugSource in @("chembl", "opentargets", "pubchem")) {
    if (-not (Get-ChildItem -LiteralPath $drugDetailsDir -Directory -Filter "*_$drugSource" -ErrorAction SilentlyContinue | Select-Object -First 1)) {
        throw "Drug annotations are missing for $drugSource in $drugDetailsDir"
    }
}
$requiredInputs = @(
    "--drug-network-file", (Join-Path $currentInputDir "drug_network.csv"),
    "--gene-file", (Join-Path $currentInputDir "gene_data.tsv"),
    "--interaction-network-file", (Join-Path $currentInputDir "network\string.human_links_v12_0_min900.Ensembl_edges.tsv"),
    "--tx2gene-file", (Join-Path $currentInputDir "tx2gene.tsv"),
    "--gsea-files") + @($gseaReports | ForEach-Object FullName) + @(
    "--gmt-files") + $gmtFiles + @(
    "--gmt-list-file", $gmtListFile,
    "--patient-id", $PatientId,
    "--padj-cutoff", $PadjCutoff.ToString([Globalization.CultureInfo]::InvariantCulture),
    "--log2fc-cutoff", $Log2fcCutoff.ToString([Globalization.CultureInfo]::InvariantCulture),
    "--interaction-network-max-neighbors", $MaxNeighbors
)
$requiredFiles = @(
    (Join-Path $currentInputDir "drug_network.csv"),
    (Join-Path $currentInputDir "gene_data.tsv"),
    (Join-Path $currentInputDir "network\string.human_links_v12_0_min900.Ensembl_edges.tsv"),
    (Join-Path $currentInputDir "tx2gene.tsv"),
    $gmtListFile) + $gmtFiles + @($gseaReports | ForEach-Object FullName)
foreach ($requiredFile in $requiredFiles) {
    if (-not (Test-Path -LiteralPath $requiredFile)) { throw "Required development input is missing: $requiredFile" }
}

$optionalInputs = @()
foreach ($optionalInput in @(
    @("--clinreport-dir", $drugDetailsDir),
    @("--wes-vcf-file", (Join-Path $currentInputDir "wes.pcgr.grch38.pass.vcf.gz")),
    @("--tf-file", (Join-Path $currentInputDir "tf.tsv"))
)) {
    if (Test-Path -LiteralPath $optionalInput[1]) {
        $optionalInputs += $optionalInput
    } else {
        Write-Host "Optional input not found; disabling $($optionalInput[0])."
    }
}

if (Test-Path -LiteralPath $outputPath) {
    Remove-Item -LiteralPath $outputPath -Recurse -Force
}

Write-Host "Preparing $outputPath"
& $RscriptPath (Join-Path $projectDir "run.R") "--prepare" $outputPath @requiredInputs @optionalInputs
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

if ($PrepareOnly) {
    Write-Host "Prepared. Start it with: & '$RscriptPath' '$outputPath\run.R'"
    exit 0
}

Write-Host "Starting ClinReport at http://${BindHost}:$Port (Ctrl+C stops it)"
& $RscriptPath (Join-Path $outputPath "run.R") "--host" $BindHost "--port" $Port
exit $LASTEXITCODE
