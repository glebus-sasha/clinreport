# ClinReport — modular Shiny app

This is a structural refactor of the original single-file Shiny app. The application logic is intentionally kept behaviorally close to the original code; the main change is separation into configuration, data loading, domain helpers, UI components/pages, and server logic.

## Run

From the project root:

```r
shiny::runApp()
```

`app.R` is the composition root and sources the modules in dependency order.

## Prepare a runnable application

`clinreport --prepare` writes a standalone application directory and exits. The
output contains `app.R`, `run.R`, `R/`, `www/`, and copied `data/`. The exact
arguments supplied for the report are embedded in the application itself; there
is no JSON or external configuration file. Run it later in the **same R environment
or container**. No image build or running Shiny server is needed during preparation.

From this repository (using the included local data):

```sh
Rscript run.R --prepare results/clinreport \
  --drug-network-file raw/all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv \
  --clinreport-dir raw/clinreport \
  --gene-file raw/carcinoma_vs_normal_gene_names_added.tsv \
  --wes-vcf-file raw/R_PTA_22.pcgr.grch38.pass.vcf.gz \
  --interaction-network-file raw/network/string.human_links_v12_0_min900.Ensembl.edges.tsv \
  --tx2gene-file raw/tx2gene.tsv \
  --gsea-carcinoma-report-file raw/gsea/carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_carcinoma.tsv \
  --gsea-normal-report-file raw/gsea/carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_normal.tsv \
  --pathway-gmt-file raw/gsea/carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_h_all_v2026_1_Hs_symbols.gmt \
  --tf-file raw/carcinoma_vs_normal_significant_tfs.tsv \
  --patient-id PATIENT-001
Rscript results/clinreport/run.R
```

Open http://127.0.0.1:3838. The launcher works from any working directory.
An existing output directory is rejected to avoid overwriting a previous report.
If copying fails, the partial output is retained for inspection; use a new output
directory after fixing the error.

### Windows development shortcut

For day-to-day development with the files in `raw/`, run this from PowerShell in
the repository root:

```powershell
.\dev.ps1
```

If PowerShell blocks local scripts, double-click `dev.bat` in Explorer or run:

```powershell
.\dev.bat
```

The batch launcher applies `ExecutionPolicy Bypass` only to the PowerShell process
it starts; it does not change the system or user execution policy.

It recreates `results/dev-report` and starts the result immediately at
http://127.0.0.1:3838. It deliberately replaces only that directory. Stop Shiny
with `Ctrl+C`. Use `-PrepareOnly` to rebuild without starting the server, or set
`-Port`, `-PatientId`, `-PadjCutoff`, `-Log2fcCutoff`, `-MaxNeighbors`, and
`-OutputDir` as needed; use `-BindHost` to change the listening address.
`-OutputDir` is restricted to `results/`.

The launcher finds the newest R installation under `C:\Program Files\R`.
If R is installed elsewhere, provide its executable explicitly:

```powershell
.\dev.bat -RscriptPath "C:\path\to\Rscript.exe"
```

All input paths, display metadata, and thresholds are command-line options. Paths
are relative to the calling working directory. Seven core input paths are required
for `--prepare` and `--check`. `--clinreport-dir`, `--wes-vcf-file`, `--tf-file`,
and `--dorothea-file` are optional. Omitting the first three disables the
corresponding analysis; they also accept an explicit empty value, such as
`--wes-vcf-file=`. Metadata and thresholds use the documented defaults unless
explicitly supplied.

| Parameter | Required for `--prepare` / `--check` | Default |
| --- | --- | --- |
| `--drug-network-file` | Yes | — |
| `--clinreport-dir` | No | Annotation tabs are unavailable |
| `--gene-file` | Yes | — |
| `--wes-vcf-file` | No | WES views and variant evidence are unavailable |
| `--interaction-network-file` | Yes | — |
| `--tx2gene-file` | Yes | — |
| `--gsea-carcinoma-report-file` | Yes | — |
| `--gsea-normal-report-file` | Yes | — |
| `--pathway-gmt-file` | Yes | — |
| `--tf-file` | No | TF analysis is unavailable |
| `--dorothea-file` | No | Downloaded at application startup |
| `--patient-id` | No | `PATIENT-001` |
| `--wes-display-name` | No | `WES · PCGR` |
| `--network-display-name` | No | `Drug–gene interaction network` |
| `--interaction-network-name` | No | `STRING` |
| `--pathway-collection-name` | No | `Hallmark` |
| `--pathway-id-prefix` | No | `HALLMARK_` |
| `--padj-cutoff` | No | `0.05` |
| `--log2fc-cutoff` | No | `1.0` |
| `--interaction-network-max-neighbors` | No | `40` |

```sh
Rscript run.R --help
Rscript run.R --check [all seven required input arguments]
Rscript run.R --prepare results/patient-002 [all seven required input arguments] --patient-id PATIENT-002 --padj-cutoff 0.01
Rscript results/patient-002/run.R --port 3839
```

`--check` checks file existence and parameter values, without loading the tables,
downloading resources, or starting Shiny. Table schemas are checked by the app's
existing loaders at startup. Required data: drug network CSV, gene expression TSV,
interaction network TSV, transcript mapping TSV, two GSEA reports and GMT. Drug
annotations, PCGR VCF.gz, and TF results are optional. Supply
`--dorothea-file /path/to/dorothea_hs.rda` to bundle a local resource and avoid its
startup download. Otherwise the existing internet download remains in place.

## Nextflow and Apptainer

The updated Dockerfile installs the Shiny dependencies and registers
`clinreport` as the image entrypoint. It does not include patient data. Build and
publish the image separately to Quay, then pin its version or digest in the pipeline.

`examples/clinreport.nf` is an example DSL2 module with explicitly staged inputs.
Set `params.clinreport_container` (for example `quay.io/your-org/clinreport:1.0.0`)
and `params.outdir`, include the module, and pass the seven input channels in the
declared order. Add staged inputs and `--clinreport-dir`, `--wes-vcf-file`, or
`--tf-file` when those optional data sources are available. Add metadata and
threshold arguments directly to the command. For
multiple samples, use unique per-sample published directories.
The module copies its output so it can survive removal of Nextflow's `work/`.

Run a generated report using the same container, from the report directory:

```sh
apptainer exec --bind "$PWD:/report" --pwd /report clinreport.sif clinreport --host 127.0.0.1
# Or use the published OCI image directly:
apptainer exec --bind "$PWD:/report" --pwd /report docker://quay.io/your-org/clinreport:1.0.0 clinreport --host 127.0.0.1
```

The server defaults to loopback. On a remote host, forward port 3838 with SSH to
view it locally, or explicitly select `--host` for your deployment.

Packaging regression test (base R and jsonlite; no network):
`Rscript tests/cli_integration.R`. This checks explicit arguments, relocation,
nested file copies, invalid inputs and overwrite protection.

## Project structure

- `R/config/` — paths, thresholds, validation
- `R/data/` — CSV/TSV/JSON loading and discovery
- `R/domain/` — drug identifiers and gene logic
- `R/utils/` — generic helpers and URL builders
- `R/ui/` — reusable components, tabs, pages, layout
- `R/server/` — Shiny server orchestration and volcano plot
- `www/css/` — application stylesheet


## Drug / STRING network visualization

The main screen includes an interactive `visNetwork` graph containing the
selected drug, its target genes, and one-hop STRING neighbors. Drug targets
are blue, STRING neighbors are grey, and the selected drug is a dark diamond.
Hover over a node to view its Ensembl identifier; displayed labels use gene
symbols whenever they are available.

The graph is limited to one-hop neighbors and `network_max_neighbors` (default 40) to keep the browser responsive. The `tx2gene.tsv` mapping is used as a fallback source for gene display names; the DE gene table is preferred when a name is available there.

Install the additional package if needed:

```r
install.packages("visNetwork")
```

## TF overlay

Switch the network side panel between Pathways and TF. Select individual imported TFs or use All TFs / Clear. The TF tab contains the raw activity statistics; clicking a row toggles that TF. Activity differences are not RNA log2 fold changes. FDR determines significance, regardless of the input filename.

Selecting a drug automatically selects all pathways containing its targets and all imported TFs regulating its targets. Use All pathways / Clear to change the pathway selection. Target-linked genes, next to Direction, is enabled by default in both overlay modes: pathway genes are limited to drug targets and their displayed one-hop interaction neighbours; TF overlays additionally retain regulators of drug targets, with regulon links limited to that same drug subgraph. Turn it off to display all genes in the selected pathways or regulons. The filter preserves your pathway and TF selections.

The optional TF analysis input is configured with `--tf-file` (columns TF, logFC, AveExpr, t, P.Value, adj.P.Val). When it is omitted, the TF overlay and tab are unavailable. No preview-project files, saved R sessions, or sample activity matrices are used. TF targets are restricted to genes present in the RNA table.

Unless `dorothea_file` is supplied, at each app startup the official human DoRothEA resource is downloaded from https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda and decompressed in memory. A/B/C signed interactions are retained. No automatic disk cache is created; a running process holds the resource in memory. This current resource may differ from the version used for the original activity inference. If loading fails, the app and TF table remain available and the panel reports the failure.

TFs are triangles. Drug-target links remain blue; STRING links remain grey dashed. TF activation uses green arrows and repression uses pink dashed links with a terminal bar. WES centre dots are preserved. In Direction mode selected TFs use activity difference while other genes use RNA expression.

Integration check (requires internet and installed app dependencies): `Rscript --vanilla tests/tf_integration.R`.
