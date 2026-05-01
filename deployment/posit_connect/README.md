# Deploying AttieLab-QTL-Viewer to Posit Connect

Posit Connect (https://posit.co/products/enterprise/connect/) is the
recommended deployment target for institutional Shiny hosting. This document
walks through publishing the app from a GitHub repository to a Connect
instance.

## Prerequisites

- A Posit Connect account with publisher permissions on your target server.
- The repository pushed to GitHub (or another Git host that Connect can reach).
- An S3-compatible bucket containing the data files referenced by the app.
- An IAM user (or service account) with read-only access to that bucket.

## 1. Prepare the repo

From the repo root, create an `renv` lockfile so Connect can reproduce the
package environment:

```r
install.packages("renv")
renv::init(bare = TRUE)
renv::install(c(
  "shiny", "bslib", "dplyr", "purrr", "stringr",
  "qtl2", "qtl2fst",
  "ggplot2", "plotly", "ggiraph",
  "data.table", "DT", "fst",
  "shinyjs", "shinycssloaders", "shinybusy", "waiter",
  "promises", "future", "future.apply",
  "paws.storage", "filelock",
  "writexl", "fontawesome", "reshape2", "htmltools", "RSQLite"
))
renv::snapshot()
```

Generate a Connect manifest (regenerate any time `renv.lock` or `app.R` changes):

```r
rsconnect::writeManifest(appDir = ".")
```

Commit `renv.lock`, `manifest.json`, and `.Rprofile`.

## 2. Publish to Connect

1. In the Connect UI → **Content → Publish → From Git**.
2. Repository URL: the GitHub HTTPS URL.
3. Branch: `main` (or whichever you want Connect to track).
4. Subdirectory: leave blank if the repo root is the app, otherwise the
   subdirectory containing `app.R`.
5. Entrypoint: `app.R`.
6. Click **Deploy**.

## 3. Bucket access

Create a dedicated IAM user (or equivalent) with **read-only** access scoped
to just the bucket/prefix this app uses. A minimal AWS IAM policy:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "ReadObjects",
      "Effect": "Allow",
      "Action": ["s3:GetObject", "s3:HeadObject"],
      "Resource": "arn:aws:s3:::<your-bucket>/<your-prefix>/*"
    },
    {
      "Sid": "ListPrefix",
      "Effect": "Allow",
      "Action": ["s3:ListBucket"],
      "Resource": "arn:aws:s3:::<your-bucket>",
      "Condition": {
        "StringLike": { "s3:prefix": ["<your-prefix>/*"] }
      }
    }
  ]
}
```

Generate an access key pair for this principal. Treat the secret as a
credential — never paste it into committed source.

## 4. Per-content environment variables

In the Connect UI → **Content → Settings → Vars** add the following.
Mark the AWS keys as encrypted.

| Variable | Value | Notes |
|---|---|---|
| `QTLAPP_DATA_BACKEND` | `s3` | required |
| `QTLAPP_S3_BUCKET` | `<your-bucket>` | no `s3://`, no trailing slash |
| `QTLAPP_S3_PREFIX` | `<your-prefix>/` | trailing slash OK; empty if data is at bucket root |
| `QTLAPP_S3_ENDPOINT` | `https://<your-s3-endpoint>` | leave empty for AWS regional default |
| `QTLAPP_S3_REGION` | `us-east-2` | whichever region applies |
| `QTLAPP_S3_FORCE_PATH_STYLE` | `true` | required for MinIO/Ceph; harmless on AWS |
| `AWS_ACCESS_KEY_ID` | *(encrypted)* | from step 3 |
| `AWS_SECRET_ACCESS_KEY` | *(encrypted)* | from step 3 |
| `QTLAPP_CACHE_DIR` | `/tmp/qtl_viewer_cache` | override to a persistent path if your Connect tier supports it |
| `QTLAPP_CACHE_MAX_GB` | `15` | LRU eviction kicks in above this |
| `QTLAPP_ETAG_CHECK_MODE` | `on_stale` | `never` / `on_stale` / `always` |
| `QTLAPP_PREWARM` | `true` | set to `false` to skip boot-time downloads |

## 5. Runtime settings

Recommended starting values; tune with telemetry:

- **Min/Max processes**: 1 / 2
- **Max connections per process**: 5
- **Idle timeout per process**: 3600 s (keeps cache warm between users)
- **Initial timeout**: 180 s (prewarm runs here)
- **Process memory**: 8–16 GB if your tier allows. The cross object and
  pheno FSTs hold multi-hundred-MB working sets.

## 6. Verification

After deploy, open the content URL and tail the Connect logs while it boots:

```
qtl_config: backend=s3 ...
prewarm_on_boot: starting.
prewarm_on_boot: done in N s (M files cached).
```

In the app, pick a trait and trigger a scan. Watch for `s3_get` lines on
first use; subsequent scans of the same trait should be cache hits.

## 7. Updating the dataset

If you change a file in the bucket:

- With `QTLAPP_ETAG_CHECK_MODE=on_stale` or `always`, the app refetches on
  the next access once the TTL elapses (default 10 minutes).
- To force immediate refresh, restart the content from the Connect UI; the
  cache will re-verify ETags on next access.
- To force a complete rebuild, clear `QTLAPP_CACHE_DIR` before restart.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `QTLAPP_DATA_BACKEND=s3 but QTLAPP_S3_BUCKET is empty` at boot | Missing env var | Set it in the Vars tab; restart |
| `s3_head(...) failed: 403` | IAM policy too narrow or wrong keys | Re-check policy; verify access keys |
| First request is slow (90+ s) | Prewarm running on cold cache | Subsequent requests are fast |
| Memory grows unbounded | Cache cap too high or too many `future` workers | Lower `QTLAPP_CACHE_MAX_GB` and/or `QTLAPP_FUTURE_WORKERS`; restart |
