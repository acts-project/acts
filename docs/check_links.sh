#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LINKCHECK_IGNORE_REPO="${LINKCHECK_IGNORE_REPO:-acts-project/linkcheck-ignore}"
LINKCHECK_IGNORE_REF="${LINKCHECK_IGNORE_REF:-main}"
LINKCHECK_IGNORE_PATH="${LINKCHECK_IGNORE_PATH:-data.json}"
LINKCHECK_IGNORE_URL="${LINKCHECK_IGNORE_URL:-https://raw.githubusercontent.com/${LINKCHECK_IGNORE_REPO}/${LINKCHECK_IGNORE_REF}/${LINKCHECK_IGNORE_PATH}}"
LINKCHECK_THREADS="${LINKCHECK_THREADS:-4}"
LINKCHECK_TIMEOUT="${LINKCHECK_TIMEOUT:-20}"

if [[ -n "${LINKCHECK_START_URL:-}" ]]; then
  START_URL="${LINKCHECK_START_URL}"
else
  LINKCHECK_START_PATH="${LINKCHECK_START_PATH:-${REPO_ROOT}/build/docs/html/index.html}"
  START_URL="${LINKCHECK_START_PATH}"
fi

if command -v gh >/dev/null 2>&1; then
  if IGNORE_JSON="$(
    gh api \
      -H "Accept: application/vnd.github.raw" \
      "/repos/${LINKCHECK_IGNORE_REPO}/contents/${LINKCHECK_IGNORE_PATH}?ref=${LINKCHECK_IGNORE_REF}" \
      2>/dev/null
  )"; then
    echo "Loaded ignore list via gh api from ${LINKCHECK_IGNORE_REPO}@${LINKCHECK_IGNORE_REF}:${LINKCHECK_IGNORE_PATH}"
  else
    echo "gh api failed, falling back to curl from ${LINKCHECK_IGNORE_URL}"
    IGNORE_JSON="$(curl -fsSL "${LINKCHECK_IGNORE_URL}")"
  fi
else
  IGNORE_JSON="$(curl -fsSL "${LINKCHECK_IGNORE_URL}")"
fi
jq -e 'type == "array" and all(.[]; type == "string")' <<<"${IGNORE_JSON}" >/dev/null
mapfile -t IGNORE_PATTERNS < <(jq -r '.[]' <<<"${IGNORE_JSON}")

IGNORE_ARGS=(--ignore-url='^mailto:')
for pattern in "${IGNORE_PATTERNS[@]}"; do
  IGNORE_ARGS+=(--ignore-url="${pattern}")
done

echo "Loaded ${#IGNORE_PATTERNS[@]} external link ignore patterns"

uvx --from LinkChecker linkchecker \
  --check-extern \
  --threads="${LINKCHECK_THREADS}" \
  --timeout="${LINKCHECK_TIMEOUT}" \
  "${IGNORE_ARGS[@]}" \
  "${START_URL}"
