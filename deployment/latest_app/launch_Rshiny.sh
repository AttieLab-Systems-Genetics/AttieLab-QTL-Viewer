#!/bin/bash
# AttieLab-QTL-Viewer launch script — builds the Docker image and runs the
# container with S3-backend environment variables.
#
# Usage:
#   ./deployment/latest_app/launch_Rshiny.sh                       # reads ./.env
#   ./deployment/latest_app/launch_Rshiny.sh /path/to/my.env       # explicit env file
#
# The env file must define (see deployment/latest_app/.env.example):
#   QTLAPP_S3_ENDPOINT, QTLAPP_S3_BUCKET,
#   AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY
#
# Optional overrides:
#   QTLAPP_S3_PREFIX, QTLAPP_S3_REGION, QTLAPP_S3_FORCE_PATH_STYLE,
#   QTLAPP_CACHE_MAX_GB, QTLAPP_ETAG_CHECK_MODE,
#   HOST_PORT, CONTAINER_NAME, IMAGE_NAME, HOST_CACHE_DIR

set -euo pipefail

# --- locate the repo root (this script lives at deployment/latest_app/) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

# --- load env file ---
ENV_FILE="${1:-$REPO_ROOT/.env}"
if [[ ! -f "$ENV_FILE" ]]; then
  echo "Error: env file not found at $ENV_FILE"
  echo ""
  echo "Copy the template and fill in your S3 credentials:"
  echo "  cp deployment/latest_app/.env.example .env"
  echo "  \$EDITOR .env"
  exit 1
fi
echo "Loading env from: $ENV_FILE"
set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a

# --- required vars ---
missing=()
for v in QTLAPP_S3_ENDPOINT QTLAPP_S3_BUCKET AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY; do
  if [[ -z "${!v:-}" ]]; then missing+=("$v"); fi
done
if [[ ${#missing[@]} -gt 0 ]]; then
  echo "Error: missing required env vars in $ENV_FILE:"
  printf '  %s\n' "${missing[@]}"
  exit 2
fi

# --- defaults ---
HOST_PORT="${HOST_PORT:-3838}"
CONTAINER_NAME="${CONTAINER_NAME:-qtl-viewer}"
IMAGE_NAME="${IMAGE_NAME:-qtl-viewer:latest}"
HOST_CACHE_DIR="${HOST_CACHE_DIR:-$HOME/qtl_viewer_cache}"
QTLAPP_CACHE_MAX_GB="${QTLAPP_CACHE_MAX_GB:-20}"
QTLAPP_ETAG_CHECK_MODE="${QTLAPP_ETAG_CHECK_MODE:-on_stale}"
QTLAPP_S3_PREFIX="${QTLAPP_S3_PREFIX:-}"
QTLAPP_S3_REGION="${QTLAPP_S3_REGION:-us-east-1}"
QTLAPP_S3_FORCE_PATH_STYLE="${QTLAPP_S3_FORCE_PATH_STYLE:-true}"

mkdir -p "$HOST_CACHE_DIR"

# --- stop existing container ---
echo "Stopping any existing container '$CONTAINER_NAME'..."
docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true

# --- build ---
echo "Building image '$IMAGE_NAME'..."
docker build -t "$IMAGE_NAME" -f deployment/latest_app/Dockerfile .

# --- run ---
echo "Starting container '$CONTAINER_NAME' on host port $HOST_PORT..."
echo "  S3 endpoint:  $QTLAPP_S3_ENDPOINT"
echo "  S3 bucket:    $QTLAPP_S3_BUCKET"
echo "  S3 prefix:    ${QTLAPP_S3_PREFIX:-<root>}"
echo "  Cache dir:    $HOST_CACHE_DIR (-> /cache in container)"
echo "  Cache cap:    ${QTLAPP_CACHE_MAX_GB} GB"

docker run --memory=35g --memory-swap=35g --restart always -d \
  -p "${HOST_PORT}:3838" \
  -e QTLAPP_DATA_BACKEND=s3 \
  -e QTLAPP_S3_BUCKET \
  -e QTLAPP_S3_PREFIX \
  -e QTLAPP_S3_REGION \
  -e QTLAPP_S3_ENDPOINT \
  -e QTLAPP_S3_FORCE_PATH_STYLE \
  -e AWS_ACCESS_KEY_ID \
  -e AWS_SECRET_ACCESS_KEY \
  -e QTLAPP_CACHE_DIR=/cache \
  -e QTLAPP_CACHE_MAX_GB \
  -e QTLAPP_ETAG_CHECK_MODE \
  -e QTLAPP_FUTURE_WORKERS \
  -v "$HOST_CACHE_DIR:/cache" \
  --name "$CONTAINER_NAME" "$IMAGE_NAME" >/dev/null

echo ""
echo "Container: $CONTAINER_NAME"
echo "App URL:   http://$(hostname -f 2>/dev/null || hostname):${HOST_PORT}/"
echo ""
echo "Tail logs:     docker logs -f $CONTAINER_NAME"
echo "Container RAM: docker stats --no-stream $CONTAINER_NAME"
echo "Stop:          docker rm -f $CONTAINER_NAME"
