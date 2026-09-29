#!/usr/bin/env bash
set -euo pipefail

# Publish completion only after startup and statistics initialization finish.
status=0
if sccache --start-server && sccache --zero-stats; then
  status=0
else
  status=$?
fi
echo "$status" > "$RUNNER_TEMP/sccache-startup.status.tmp"
mv "$RUNNER_TEMP/sccache-startup.status.tmp" "$RUNNER_TEMP/sccache-startup.status"
