# AttieLab-QTL-Viewer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

An interactive R/Shiny web application for exploring quantitative trait loci
(QTL) in Diversity Outbred mouse data — genome-wide LOD scans, Manhattan and
cis/trans plots, peak exploration with founder allele effects, sex- and
diet-stratified comparisons, trait search, correlation, profile plots,
mediation analysis, and SNP association.

> **Live application**: <https://connect.doit.wisc.edu/content/08a89532-32d3-4b0b-bff3-a16f3a05c926>
>
> Hosted on UW–Madison's Posit Connect (Research Computing Initiative). The
> source for the deployed version is this repository.

---

## Citing this software

If you use AttieLab-QTL-Viewer in your work, please cite both the software
archive and the manuscript that describes it:

- **Software archive (Zenodo)**: <https://doi.org/10.5281/zenodo.19961001>
  (concept DOI — always resolves to the latest version)
- **Manuscript**: *citation to be added.*

A machine-readable citation is available in [`CITATION.cff`](CITATION.cff).

## Features

- **LOD scans** — fast genome-wide scans across multiple datasets, with
  additive and interactive (sex / diet / sex × diet) modes and on-the-fly
  difference plots.
- **Peak exploration** — Manhattan and cis/trans plots; click any peak to
  inspect its founder allele effects.
- **Trait search** — large-scale fuzzy search across genes, isoforms,
  clinical traits, lipids, plasma metabolites, and splice junctions.
- **Correlation & profile plots** — trait-trait correlations against any
  matched dataset; per-trait profile plots faceted by diet/sex.
- **Mediation analysis** — explore candidate mediators between QTL peaks
  and downstream traits.
- **SNP association** — local SNP scans within a peak window, with founder
  allele effects and gene context.

## Architecture

A modular Shiny app with focused server modules in [`R/`](R/), a single
top-level [`app.R`](app.R) entry point, and a pluggable data backend that can
read either from a local directory or any S3-compatible bucket (AWS, MinIO,
Ceph, on-prem). See [`docs/architecture.md`](docs/architecture.md) for a
module-by-module map.

## Quick start (local development)

Requires R ≥ 4.2.0 and a few system libraries (curl, ssl, xml2, blas,
lapack, sodium, gfortran). On Ubuntu/Debian those are
`libcurl4-gnutls-dev libssl-dev libxml2-dev libblas-dev liblapack-dev gfortran libsodium-dev cmake`.

```r
# Install R-package dependencies (one-time)
source("install_packages.R")

# Point the app at a local data root, then launch
Sys.setenv(QTLAPP_DATA_BACKEND = "local",
           QTLAPP_DATA_ROOT    = "<your-data-root>")
shiny::runApp("app.R", port = 3838)
```

If you don't have a local data root, the simpler path is to point the app
at an S3-compatible bucket — see the [Deployment](#deployment) section.

## Deployment

Two deployment paths are supported and documented:

- **Docker on a server** — see [`deployment/latest_app/README.md`](deployment/latest_app/README.md).
  One-command launcher that builds the image and starts the container with
  S3-backend env vars supplied via a `.env` file.
- **Posit Connect** — see [`deployment/posit_connect/README.md`](deployment/posit_connect/README.md).
  Git-backed deployment using `renv` + `manifest.json`; env vars supplied
  via the Connect Vars tab.

Both paths use the same data layer (`R/config.R`, `R/s3_io.R`,
`R/s3_cache.R`, `R/data_source.R`, `R/prewarm.R`) and the same `.env`-style
configuration. See [`docs/deployment.md`](docs/deployment.md) for an
overview comparing them.

## Project layout

```
.
├── app.R                       Shiny entry point
├── R/                          Shiny modules + data layer
├── inst/shinyApp/              Alternate entry point + server config
├── deployment/                 Docker & Posit Connect deployment
├── utils/                      Maintenance and data-prep utilities
├── scripts/                    One-off data conversion scripts
├── docs/                       Architecture and deployment guides
├── DESCRIPTION                 R package metadata
├── NAMESPACE                   Auto-generated exports
├── LICENSE                     MIT
└── install_packages.R          Bootstrap install script
```

## Contributing

We welcome bug reports, feature requests, and pull requests. Please read
[`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or PR. By
participating in this project you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Authors

- **Charles Opara** (maintainer) — University of Wisconsin–Madison
- **Brian S Yandell** — University of Wisconsin–Madison
- **Chris Emfinger** — University of Wisconsin–Madison
- **Kalynn Willis** — University of Wisconsin–Madison
- **Kelly Mitok** — University of Wisconsin–Madison

## Acknowledgments

This work is supported by the Attie Lab at the University of Wisconsin–Madison.

## License

Released under the [MIT License](LICENSE).
