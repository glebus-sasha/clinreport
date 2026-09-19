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
