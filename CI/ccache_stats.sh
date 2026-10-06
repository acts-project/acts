#!/usr/bin/env bash
# Emit a record only after both collection and JSON validation succeed.
# Reporting must never change the build outcome.
set -uo pipefail

if ! command -v jq >/dev/null 2>&1; then
  echo '::warning::ccache JSON statistics unavailable: jq is not installed' >&2
  exit 0
fi
if ! stats=$(ccache --print-stats --format=json); then
  echo '::warning::ccache JSON statistics unavailable: collection failed (JSON support required)' >&2
  exit 0
fi
if ! stats=$(printf '%s\n' "$stats" | jq -ce 'select(type == "object" and has("cache_miss") and has("direct_cache_hit"))'); then
  echo '::warning::ccache JSON statistics unavailable: invalid statistics' >&2
  exit 0
fi
printf 'CCACHE_STATS %s\n' "$stats"
