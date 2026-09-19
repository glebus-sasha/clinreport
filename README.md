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
