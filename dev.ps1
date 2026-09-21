[CmdletBinding()]
param(
    [string]$OutputDir = "results/dev-report",
    [int]$Port = 3838,
    [string]$BindHost = "127.0.0.1",
    [string]$PatientId = "PATIENT-001",
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
$inputs = @(
    "--drug-network-file", (Join-Path $rawDir "all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv"),
    "--clinreport-dir", (Join-Path $rawDir "clinreport"),
    "--gene-file", (Join-Path $rawDir "carcinoma_vs_normal_gene_names_added.tsv"),
    "--wes-vcf-file", (Join-Path $rawDir "R_PTA_22.pcgr.grch38.pass.vcf.gz"),
    "--interaction-network-file", (Join-Path $rawDir "network\string.human_links_v12_0_min900.Ensembl.edges.tsv"),
    "--tx2gene-file", (Join-Path $rawDir "tx2gene.tsv"),
    "--gsea-carcinoma-report-file", (Join-Path $rawDir "gsea\carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_carcinoma.tsv"),
    "--gsea-normal-report-file", (Join-Path $rawDir "gsea\carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_normal.tsv"),
    "--pathway-gmt-file", (Join-Path $rawDir "gsea\carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_h_all_v2026_1_Hs_symbols.gmt"),
    "--tf-file", (Join-Path $rawDir "carcinoma_vs_normal_significant_tfs.tsv"),
    "--patient-id", $PatientId,
    "--padj-cutoff", $PadjCutoff.ToString([Globalization.CultureInfo]::InvariantCulture),
    "--log2fc-cutoff", $Log2fcCutoff.ToString([Globalization.CultureInfo]::InvariantCulture),
    "--interaction-network-max-neighbors", $MaxNeighbors
)
foreach ($index in 1, 3, 5, 7, 9, 11, 13, 15, 17, 19) {
    if (-not (Test-Path -LiteralPath $inputs[$index])) { throw "Required development input is missing: $($inputs[$index])" }
}

if (Test-Path -LiteralPath $outputPath) {
    Remove-Item -LiteralPath $outputPath -Recurse -Force
}

Write-Host "Preparing $outputPath"
& $RscriptPath (Join-Path $projectDir "run.R") "--prepare" $outputPath @inputs
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

if ($PrepareOnly) {
    Write-Host "Prepared. Start it with: & '$RscriptPath' '$outputPath\run.R'"
    exit 0
}

Write-Host "Starting ClinReport at http://${BindHost}:$Port (Ctrl+C stops it)"
& $RscriptPath (Join-Path $outputPath "run.R") "--host" $BindHost "--port" $Port
exit $LASTEXITCODE
