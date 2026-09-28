#!/usr/bin/env bash
set -euo pipefail

# Do not invoke sccache here until startup finishes: it could launch a second
# daemon without the publishing credentials supplied to the original process.
deadline=$((SECONDS + 65))
while [[ ! -f "$RUNNER_TEMP/sccache-startup.status" && $SECONDS -lt $deadline ]]; do
  sleep 1
done

if [[ -f "$RUNNER_TEMP/sccache-startup.status" ]] \
  && [[ "$(cat "$RUNNER_TEMP/sccache-startup.status")" == 0 ]]; then
  echo 'launcher=sccache' >> "$GITHUB_OUTPUT"
  cat "$RUNNER_TEMP/sccache-startup-client.log"
else
  echo 'launcher=' >> "$GITHUB_OUTPUT"
  echo '::warning::sccache startup failed; building without compiler caching.'
  echo 'sccache unavailable: building without compiler caching.' >> "$GITHUB_STEP_SUMMARY"
  echo '::group::sccache startup diagnostics'
  cat "$RUNNER_TEMP/sccache-startup-client.log" "$RUNNER_TEMP/sccache-startup.log"
  echo '::endgroup::'
fi
