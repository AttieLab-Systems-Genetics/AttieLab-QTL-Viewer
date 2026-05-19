# Contributing

Thank you for your interest in improving AttieLab-QTL-Viewer. We welcome
bug reports, feature requests, documentation fixes, and code contributions.

## Reporting bugs

Use the [issue tracker](https://github.com/AttieLab-Systems-Genetics/AttieLab-QTL-Viewer/issues)
and the **Bug report** template. A useful report includes:

1. The exact thing you did (which dataset, which trait, which tab).
2. What you expected to happen.
3. What actually happened, including any error message or relevant log
   snippet from `docker logs <container>` if you're running locally.
4. Your environment: R version, OS, deployment target (Docker / Posit Connect).

## Requesting features

Open an issue using the **Feature request** template. Describe the user
problem the feature would solve before describing a specific implementation.

## Pull requests

1. Fork the repository and create a branch from `main`.
2. Make your change. Keep diffs focused — one feature or fix per PR.
3. If you change behavior, update the relevant docs in `docs/`.
4. Make sure `R CMD check` passes (the GitHub Actions workflow runs this
   on every push).
5. Open a PR against `main`. Use the PR template; reference any issues.

### Coding conventions

- Follow the existing style in [`R/`](R/). Roughly: 1-space indent (yes,
  really — that's the existing convention), `<-` for assignment,
  `lintr` config in [`.lintr`](.lintr).
- New Shiny modules should follow the existing pattern: a paired
  `<thing>UI` / `<thing>Server` function, no global state, all data access
  routed through `local_path()` from [`R/data_source.R`](R/data_source.R).
- New file reads should use `local_path("relative/key")` rather than any
  absolute path. The data layer handles caching and S3 fetch transparently.

### Local dev environment

```sh
# Install package dependencies
Rscript install_packages.R

# Point at a local data root (or set up S3 — see docs/deployment.md)
export QTLAPP_DATA_BACKEND=local
export QTLAPP_DATA_ROOT=<your-data-root>

# Run the app
R -e 'shiny::runApp("app.R", port = 3838)'
```

## Code of Conduct

By participating, you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).

## License

Contributions are accepted under the [MIT License](LICENSE) of this project.
