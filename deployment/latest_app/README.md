# Docker deployment

This directory contains the artifacts to run AttieLab-QTL-Viewer as a Shiny
Server inside Docker, with data fetched from any S3-compatible bucket
(AWS S3, MinIO, Ceph, on-prem, etc.).

## TL;DR

```sh
# 1. From the repo root, copy the env template and fill in your values
cp deployment/latest_app/.env.example .env
$EDITOR .env

# 2. Launch
./deployment/latest_app/launch_Rshiny.sh
```

The app comes up on `http://<your-host>:3838/` by default.

## Configuration

All deployment configuration lives in a single `.env` file at the repository
root. See [.env.example](.env.example) for the template; required keys:

- `QTLAPP_S3_ENDPOINT` — full HTTPS URL of your S3 endpoint
- `QTLAPP_S3_BUCKET` — bucket name
- `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` — credentials (read-only is
  sufficient)

Optional knobs (defaults in the launch script):

- `HOST_PORT` (default `3838`)
- `CONTAINER_NAME` (default `qtl-viewer`)
- `IMAGE_NAME` (default `qtl-viewer:latest`)
- `HOST_CACHE_DIR` (default `~/qtl_viewer_cache`) — persistent on-disk cache
  for files fetched from S3
- `QTLAPP_S3_PREFIX`, `QTLAPP_S3_REGION`, `QTLAPP_S3_FORCE_PATH_STYLE`,
  `QTLAPP_CACHE_MAX_GB`, `QTLAPP_ETAG_CHECK_MODE`, `QTLAPP_FUTURE_WORKERS`

## Common operations

```sh
# Logs
docker logs -f qtl-viewer

# Shell into the container
docker exec -it qtl-viewer /bin/bash

# Restart (reuses cache)
docker restart qtl-viewer

# Stop + remove
docker rm -f qtl-viewer

# Full rebuild (forgets image cache, reinstalls R packages)
docker rm -f qtl-viewer
docker rmi qtl-viewer:latest
./deployment/latest_app/launch_Rshiny.sh

# Clear the file cache without rebuilding
docker rm -f qtl-viewer
rm -rf "$HOME/qtl_viewer_cache"
./deployment/latest_app/launch_Rshiny.sh
```

## Running multiple instances

Use separate env files + container names. Each instance also needs its own
`HOST_CACHE_DIR`:

```sh
# .env.prod with HOST_PORT=3838 CONTAINER_NAME=qtl-viewer-prod IMAGE_NAME=qtl-viewer:prod HOST_CACHE_DIR=$HOME/qtl_viewer_cache_prod
# .env.dev  with HOST_PORT=3839 CONTAINER_NAME=qtl-viewer-dev  IMAGE_NAME=qtl-viewer:dev  HOST_CACHE_DIR=$HOME/qtl_viewer_cache_dev

./deployment/latest_app/launch_Rshiny.sh .env.prod
./deployment/latest_app/launch_Rshiny.sh .env.dev
```

## Posit Connect

For deployment on Posit Connect (managed or self-hosted), see
[../posit_connect/README.md](../posit_connect/README.md). The same env-var
set applies; values are entered in Connect's content `Vars` tab instead of
a `.env` file.
