# ClinReport — modular Shiny app

This is a structural refactor of the original single-file Shiny app. The application logic is intentionally kept behaviorally close to the original code; the main change is separation into configuration, data loading, domain helpers, UI components/pages, and server logic.

## Run

From the project root:

```r
shiny::runApp()
```

`app.R` is the composition root and sources the modules in dependency order.

## Project structure

- `R/config/` — paths, thresholds, validation
- `R/data/` — CSV/TSV/JSON loading and discovery
- `R/domain/` — drug identifiers, names, and gene logic
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

The only TF analysis input is `raw/carcinoma_vs_normal_significant_tfs.tsv` (TF, logFC, AveExpr, t, P.Value, adj.P.Val). No preview-project files, saved R sessions, or sample activity matrices are used. TF targets are restricted to genes present in the existing raw RNA table.

At each app startup the official human DoRothEA resource is downloaded from https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda and decompressed in memory. A/B/C signed interactions are retained. No disk cache is created; a running process holds the resource in memory. This current resource may differ from the version used for the original activity inference. If the download fails, the app and TF table remain available and the panel reports the failure.

TFs are triangles. Drug-target links remain blue; STRING links remain grey dashed. TF activation uses green arrows and repression uses pink dashed links with a terminal bar. WES centre dots are preserved. In Direction mode selected TFs use activity difference while other genes use RNA expression.

Integration check (requires internet and installed app dependencies): `Rscript --vanilla tests/tf_integration.R`.
