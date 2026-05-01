# Deployment

AttieLab-QTL-Viewer supports two deployment paths. Both share the same
underlying data layer and the same `.env`-style configuration; they differ
only in how the Shiny process is hosted.

| Path | When to use |
|---|---|
| **Docker** ([`deployment/latest_app/`](../deployment/latest_app/)) | You control a server with Docker installed. One-command `launch_Rshiny.sh` starts a Shiny Server container with all dependencies pre-installed. Best for self-hosted institutional servers. |
| **Posit Connect** ([`deployment/posit_connect/`](../deployment/posit_connect/)) | Your institution provides Posit Connect (managed or self-hosted). Git-backed deployment, no server admin required. Best for managed academic platforms. |

The detailed step-by-step instructions live in those subdirectories; this
document covers the cross-cutting concepts.

## Configuration

All deployment configuration is via environment variables. The app reads
them through [`R/config.R`](../R/config.R).

### Backend selection

| Variable | Default | Purpose |
|---|---|---|
| `QTLAPP_DATA_BACKEND` | `local` | `local` reads from a directory; `s3` reads from an S3-compatible bucket and caches to disk |
| `QTLAPP_DATA_ROOT` | `<your-data-root>` placeholder | Used only when backend is `local` |

### S3 backend

| Variable | Default | Purpose |
|---|---|---|
| `QTLAPP_S3_BUCKET` | *(empty — required)* | Bucket name (no `s3://`) |
| `QTLAPP_S3_PREFIX` | *(empty)* | Sub-path within the bucket; trailing slash optional |
| `QTLAPP_S3_ENDPOINT` | *(empty — AWS default)* | HTTPS URL for non-AWS endpoints (MinIO, Ceph, on-prem) |
| `QTLAPP_S3_REGION` | `us-east-1` | Required by paws even when ignored by the endpoint |
| `QTLAPP_S3_FORCE_PATH_STYLE` | `true` | Required for MinIO/Ceph; harmless on AWS |
| `AWS_ACCESS_KEY_ID` | *(empty — required)* | Credentials |
| `AWS_SECRET_ACCESS_KEY` | *(empty — required)* | Credentials |

### Caching

| Variable | Default | Purpose |
|---|---|---|
| `QTLAPP_CACHE_DIR` | `<tempdir>/qtlapp_cache` | On-disk cache for files fetched from S3 |
| `QTLAPP_CACHE_MAX_GB` | `15` | LRU eviction triggers above this |
| `QTLAPP_ETAG_CHECK_MODE` | `on_stale` | `never` / `on_stale` / `always` — how often to re-verify cached files against S3 |
| `QTLAPP_ETAG_CHECK_TTL_SECONDS` | `600` | Used only with `on_stale` |
| `QTLAPP_PREWARM` | `true` | Set `false` to skip boot-time small-file fetch |
| `QTLAPP_FUTURE_WORKERS` | `8` | Parallel workers for I/O-bound prefetch |

### Launch-script options (Docker only)

These are read by [`launch_Rshiny.sh`](../deployment/latest_app/launch_Rshiny.sh).

| Variable | Default | Purpose |
|---|---|---|
| `HOST_PORT` | `3838` | Host port to publish |
| `CONTAINER_NAME` | `qtl-viewer` | Docker container name |
| `IMAGE_NAME` | `qtl-viewer:latest` | Docker image tag |
| `HOST_CACHE_DIR` | `~/qtl_viewer_cache` | Host-side bind-mount target for the cache |

## Bucket layout expectations

The app expects a directory of files matching what `file_index.csv`
references. At minimum:

```
<bucket>/<prefix>/
├── file_index.csv                          required
├── gene_symbols.csv                        recommended (trait search)
├── chromosomal_sep_mm11.csv                required (axis breaks)
├── annotation_list.rds                     required
├── CHTC_dietDO_markers_RDSgrcm39.rds       required (marker positions)
├── chromosome<N>_<dataset>_..._with_symbols.fst       scan FSTs
├── chromosome<N>_..._rows.fst                         row-index sidecars
├── DO1200_<dataset>_..._peaks.csv                     peak tables
└── *_corr.fst, *_pval.fst, *_num_mice.fst             correlation tables
```

Row-index sidecars (`*_rows.fst`) make per-trait slicing fast and **must**
exist alongside their parent FSTs. They can be regenerated locally with
[`R/fst_rows.R`](../R/fst_rows.R) before uploading to the bucket.

## Memory considerations

The app pins several large files in cache (the cross RDS, pheno FSTs,
founder variants SQLite) so they survive LRU eviction once loaded. Plan
for ~8–16 GB of RAM per Shiny process and several GB of disk for the
cache directory.

## Updating data

Three options when bucket contents change:

1. **`QTLAPP_ETAG_CHECK_MODE=on_stale`** (default) — the app re-verifies
   each cached file against S3 every `QTLAPP_ETAG_CHECK_TTL_SECONDS` and
   refetches when ETags differ.
2. **`QTLAPP_ETAG_CHECK_MODE=always`** — re-verify on every access.
   Strongest consistency, highest latency.
3. **`QTLAPP_ETAG_CHECK_MODE=never`** — fastest; restart the container
   (or clear the cache directory) to pick up bucket changes.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `qtl_config: backend=local` when you expected s3 | Env vars not reaching R | On Docker via Shiny Server, pass them through `entrypoint.sh` (already done in `deployment/latest_app/`); on Posit Connect, set them in the Vars tab |
| `s3_head(...) failed: 403` | IAM credentials lack access | Re-check policy and key; verify `QTLAPP_S3_BUCKET` and `QTLAPP_S3_PREFIX` |
| `Row index not found in backend` | Missing `*_rows.fst` sidecar | Generate locally with `fst_rows()` and upload |
| Memory grows unbounded | Cache cap too high or too many workers | Lower `QTLAPP_CACHE_MAX_GB` and/or `QTLAPP_FUTURE_WORKERS`; restart |
