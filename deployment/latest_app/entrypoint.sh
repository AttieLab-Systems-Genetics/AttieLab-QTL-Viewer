#!/bin/bash
# Container entrypoint — propagates Docker `-e` env vars into a .Renviron file
# so the R process started by shiny-server (which uses `su --login shiny`,
# stripping the environment) can still see them.
#
# Without this, every QTLAPP_* and AWS_* variable would be empty inside R,
# and qtl_config() would silently fall back to backend=local.

set -e

RENV=/srv/shiny-server/.Renviron
: > "$RENV"

for var in \
    QTLAPP_DATA_BACKEND QTLAPP_DATA_ROOT \
    QTLAPP_S3_BUCKET QTLAPP_S3_PREFIX QTLAPP_S3_REGION \
    QTLAPP_S3_ENDPOINT QTLAPP_S3_FORCE_PATH_STYLE \
    QTLAPP_CACHE_DIR QTLAPP_CACHE_MAX_GB \
    QTLAPP_ETAG_CHECK_MODE QTLAPP_ETAG_CHECK_TTL_SECONDS \
    QTLAPP_PREWARM QTLAPP_FUTURE_WORKERS \
    AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION; do
  if [ -n "${!var:-}" ]; then
    printf '%s=%s\n' "$var" "${!var}" >> "$RENV"
  fi
done

chmod 640 "$RENV"
chown shiny:shiny "$RENV" 2>/dev/null || true

# v3: bind-mounted cache dir comes in with the host's UID/GID; the R process
# runs as `shiny` inside the container and won't be able to write to it
# without explicit ownership. The entrypoint runs as root, so we can fix it.
CACHE_DIR="${QTLAPP_CACHE_DIR:-/cache}"
mkdir -p "$CACHE_DIR"
if ! chown -R shiny:shiny "$CACHE_DIR" 2>/dev/null; then
  echo "qtlapp entrypoint: chown of $CACHE_DIR failed; falling back to chmod 777"
  chmod -R 777 "$CACHE_DIR"
fi

# Log how many vars made it (don't echo values — secrets)
echo "qtlapp entrypoint: wrote $(wc -l < "$RENV") env vars to $RENV"
echo "qtlapp entrypoint: backend=${QTLAPP_DATA_BACKEND:-<unset>}, bucket=${QTLAPP_S3_BUCKET:-<unset>}, prefix=${QTLAPP_S3_PREFIX:-<unset>}"
echo "qtlapp entrypoint: cache dir $CACHE_DIR ready ($(stat -c '%U:%G %a' "$CACHE_DIR" 2>/dev/null))"

exec /usr/bin/shiny-server
