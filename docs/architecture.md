# Architecture

A high-level tour of the codebase for new contributors. For deployment
instructions, see [deployment.md](deployment.md).

## Top-level entry points

- [`app.R`](../app.R) — Shiny entry point. Sources every module in `R/` in
  dependency order, configures `future::plan()` for async work, calls
  `prewarm_on_boot()`, then defines the `ui` and `server` functions.
- [`R/mainUI.R`](../R/mainUI.R) — Top-level UI composition (sidebar +
  tabset layout, dataset selector, trait search).

## Server modules (R/)

Each module is a paired `<thing>UI()` / `<thing>Server()` function. Modules
do not share state through globals — they communicate via reactive
arguments passed in from `app.R`.

| Module | What it does |
|---|---|
| [`scanPlotModule.R`](../R/scanPlotModule.R) | Genome-wide LOD scans, including additive vs. interactive overlays and on-the-fly difference plots. |
| [`peaksTableModule.R`](../R/peaksTableModule.R) | Peaks table with hover/click metadata; ties into the scan plot and allele-effects view. |
| [`manhattanPlotApp.R`](../R/manhattanPlotApp.R) | Manhattan plot view. |
| [`cisTransPlotApp.R`](../R/cisTransPlotApp.R) | Cis/trans (local/distal) plot. |
| [`correlationApp.R`](../R/correlationApp.R) | Trait-trait correlations against a chosen target dataset. |
| [`profilePlotApp.R`](../R/profilePlotApp.R) | Per-trait phenotype profile plots faceted by diet/sex. |
| [`mediationTab.R`](../R/mediationTab.R) | Mediation analysis around a selected QTL peak. |
| [`snpAssociationTab.R`](../R/snpAssociationTab.R) | Local SNP association inside a peak window. |
| [`importApp.R`](../R/importApp.R) | Wraps `import_data()` as a Shiny module. |
| [`traitApp.R`](../R/traitApp.R), [`traitProcessingModule.R`](../R/traitProcessingModule.R) | Trait search and pre-processing. |
| [`downloadApp.R`](../R/downloadApp.R) | Download handlers for plots and tables. |

## Computation and data access

| File | Role |
|---|---|
| [`R/import_data.R`](../R/import_data.R) | Loads `file_index.csv`, gene symbols, chromosome breaks, and lazy-binds annotation/marker RDS files. |
| [`R/trait_scan.R`](../R/trait_scan.R) | Per-trait LOD scan: iterates chromosome FSTs, slices by trait, returns the genome-wide scan table. |
| [`R/peak_finder.R`](../R/peak_finder.R), [`R/peak_info.R`](../R/peak_info.R), [`R/peak_info_ui.R`](../R/peak_info_ui.R) | Peak detection and metadata extraction. |
| [`R/fst_rows.R`](../R/fst_rows.R) | Builds and reads the row-index sidecar files that make per-trait FST slicing fast. |
| [`R/QTL_plot_visualizer.R`](../R/QTL_plot_visualizer.R) | Normalizes scan data for plotting (chr/position joins). |
| [`R/data_handling.R`](../R/data_handling.R) | Validation, cache-key construction, format helpers. |

## Plotting

| File | Role |
|---|---|
| [`R/ggplot_qtl_scan.R`](../R/ggplot_qtl_scan.R), [`R/ggplotly_qtl_scan.R`](../R/ggplotly_qtl_scan.R) | LOD scan plotting (ggplot2 / plotly variants). |
| [`R/ggplot_alleles.R`](../R/ggplot_alleles.R) | Founder allele effects plot. |
| [`R/plot_enhancements.R`](../R/plot_enhancements.R), [`R/ui_styles.R`](../R/ui_styles.R) | Themes, styles, custom CSS. |
| [`R/plot_null.R`](../R/plot_null.R) | Empty-state placeholder plots. |

## The data layer

The most architecturally interesting part of the codebase is the
S3-backed data layer. Every module that needs to read a file calls
`local_path("relative/key")` and gets back an absolute filesystem path —
the layer transparently handles either a local mount or an S3 bucket
fetch + cache.

```
caller (any R/*.R module)
     │
     ▼
local_path("foo.fst")            ─┐
data_exists("foo.csv")            │  R/data_source.R
as_rel_key("/path/to/foo.fst")   │
     │                            ┘
     ▼
  cache_get / cache_put           ─┐  R/s3_cache.R
  cache_pin / cache_touch          │  (local disk + LRU + ETag)
     │                            ┘
     ▼
  s3_head / s3_get_to_path         ─┐  R/s3_io.R
  s3_list                            │  (paws.storage wrapper)
     │                             ┘
     ▼
    Local mount  OR  S3 bucket
```

| File | Role |
|---|---|
| [`R/config.R`](../R/config.R) | Reads `QTLAPP_*` environment variables; exports `qtl_config()`. |
| [`R/s3_io.R`](../R/s3_io.R) | Thin paws wrapper: `s3_head()`, `s3_get_to_path()`, `s3_list()`. Per-process client cache. Retries with 404 fast-fail. |
| [`R/s3_cache.R`](../R/s3_cache.R) | Local LRU cache with size cap and pinning. RDS-backed index protected by a `filelock`. |
| [`R/data_source.R`](../R/data_source.R) | Public API: `local_path()`, `data_exists()`, `as_rel_key()`. Backend dispatch (local vs. s3) plus optional ETag re-verification. |
| [`R/prewarm.R`](../R/prewarm.R) | Boot-time fetch of small reference files so the first user-facing request doesn't pay cold S3 latency. |

If you're adding a new file read, the rule is: **never hardcode a path**.
Use `local_path("scans/foo.fst")` and the caching, S3 fetch, and
local-fallback all happen for free.

## Configuration

The app is configured exclusively through environment variables. See
[deployment.md](deployment.md#configuration) for the full list.

The two most important:

- `QTLAPP_DATA_BACKEND` — `local` (default) or `s3`
- `QTLAPP_DATA_ROOT` — for `local`; the directory that contains
  `file_index.csv` and the data files referenced by it

For `s3` backend, additionally: `QTLAPP_S3_ENDPOINT`, `QTLAPP_S3_BUCKET`,
`QTLAPP_S3_PREFIX`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`.

## Reactive patterns

- Heavy reactives (trait scan, correlation table) use `shiny::bindCache()`
  so re-selecting the same trait is instant.
- Slider inputs that drive plots are wrapped with `shiny::debounce()` to
  avoid mid-drag re-renders.
- SNP association runs in a background worker via `promises::future_promise()`
  so the UI stays responsive during long scans.

## Deployment artifacts

- [`deployment/latest_app/`](../deployment/latest_app/) — Dockerfile,
  entrypoint, launcher, and `.env.example`.
- [`deployment/posit_connect/`](../deployment/posit_connect/) — Connect
  publishing recipe.
- [`install_packages.R`](../install_packages.R) — Bootstrap installer used
  by the Dockerfile.

## Utilities

- [`utils/slim_fst_scans.R`](../utils/slim_fst_scans.R) — slim and round
  scan FSTs in place (column projection + numeric rounding).
- [`utils/convert_csv_to_fst.R`](../utils/convert_csv_to_fst.R) — batch
  convert CSV correlation files to FST.
- [`scripts/`](../scripts/), [`preprocess_pheno_data.R`](../preprocess_pheno_data.R),
  [`deployment/latest_app/update_annotation_list_for_*.R`](../deployment/latest_app/) —
  one-off data-preparation scripts. They expect a `<your-data-root>` you
  set via `QTLAPP_DATA_ROOT`. Not needed for end users.
