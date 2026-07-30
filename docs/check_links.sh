#!/usr/bin/env bash
set -euo pipefail

# Documentation link checker, backed by lychee (https://lychee.cli.rs).
#
# Two modes, selected via LINKCHECK_CHECK_EXTERN:
#   0 -> internal only: runs offline, checks that links between the built HTML
#        files resolve (including cross-file #anchors). Fast and deterministic,
#        used as the per-PR gate.
#   1 -> external too: also checks that external URLs are reachable, with
#        retries, on-disk caching and GitHub API handling to absorb transient
#        upstream failures. Used by the scheduled sweep.
#
# The set of external URLs to skip is fetched from a shared ignore list
# (acts-project/linkcheck-ignore) and passed to lychee as --exclude patterns.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LYCHEE_BIN="${LYCHEE_BIN:-lychee}"

LINKCHECK_IGNORE_REPO="${LINKCHECK_IGNORE_REPO:-acts-project/linkcheck-ignore}"
LINKCHECK_IGNORE_REF="${LINKCHECK_IGNORE_REF:-main}"
LINKCHECK_IGNORE_PATH="${LINKCHECK_IGNORE_PATH:-data.json}"
LINKCHECK_IGNORE_URL="${LINKCHECK_IGNORE_URL:-https://raw.githubusercontent.com/${LINKCHECK_IGNORE_REPO}/${LINKCHECK_IGNORE_REF}/${LINKCHECK_IGNORE_PATH}}"

# Check external links too (1) or internal links only (0).
LINKCHECK_CHECK_EXTERN="${LINKCHECK_CHECK_EXTERN:-1}"

# Check #fragments/anchors (1) in addition to target files (0). Off by default:
# the docs tree is Doxygen-generated, and Doxygen's own member-list pages link
# to detail anchors that only exist for members with a detailed description, so
# fragment checking produces false failures we cannot fix.
LINKCHECK_FRAGMENTS="${LINKCHECK_FRAGMENTS:-0}"

# Directory holding the built HTML documentation.
LINKCHECK_ROOT="${LINKCHECK_ROOT:-${REPO_ROOT}/build/docs/html}"

# External-check tuning (ignored in internal-only mode).
LINKCHECK_TIMEOUT="${LINKCHECK_TIMEOUT:-20}"
LINKCHECK_MAX_RETRIES="${LINKCHECK_MAX_RETRIES:-3}"
LINKCHECK_RETRY_WAIT="${LINKCHECK_RETRY_WAIT:-2}"
LINKCHECK_MAX_CONCURRENCY="${LINKCHECK_MAX_CONCURRENCY:-16}"
# Treat timed-out requests as OK, so a single slow-but-alive site is not an
# error (the main cause of transient false alarms).
LINKCHECK_ACCEPT_TIMEOUTS="${LINKCHECK_ACCEPT_TIMEOUTS:-1}"
# Extra status codes to accept as valid (comma separated), e.g. "429".
LINKCHECK_ACCEPT="${LINKCHECK_ACCEPT:-}"
# On-disk result cache (.lycheecache) to avoid re-hitting healthy links.
LINKCHECK_CACHE="${LINKCHECK_CACHE:-0}"
LINKCHECK_MAX_CACHE_AGE="${LINKCHECK_MAX_CACHE_AGE:-1d}"

# Optional markdown report path (used as an issue body by the scheduled job).
LINKCHECK_FAILURES_OUT="${LINKCHECK_FAILURES_OUT:-}"

if ! command -v "${LYCHEE_BIN}" >/dev/null 2>&1; then
  echo "error: lychee binary '${LYCHEE_BIN}' not found on PATH" >&2
  exit 127
fi

if [[ ! -d "${LINKCHECK_ROOT}" ]]; then
  echo "error: documentation root '${LINKCHECK_ROOT}' does not exist" >&2
  exit 1
fi
ROOT_ABS="$(cd "${LINKCHECK_ROOT}" && pwd)"

# --- load the shared external-link ignore list -----------------------------
# The ignore list is a public raw file, so a plain curl is enough (and, unlike
# an unauthenticated `gh api` call, is not subject to the API rate limit).
echo "Loading ignore list from ${LINKCHECK_IGNORE_URL}"
IGNORE_JSON="$(curl -fsSL "${LINKCHECK_IGNORE_URL}")"
jq -e 'type == "array" and all(.[]; type == "string")' <<<"${IGNORE_JSON}" >/dev/null
mapfile -t IGNORE_PATTERNS < <(jq -r '.[]' <<<"${IGNORE_JSON}")
echo "Loaded ${#IGNORE_PATTERNS[@]} external link ignore patterns"

# lychee treats --exclude values as regular expressions matched against the URL.
# mailto: links are excluded by lychee by default (opt in with --include-mail).
EXCLUDE_ARGS=(--exclude '^https?://localhost')
for pattern in "${IGNORE_PATTERNS[@]}"; do
  EXCLUDE_ARGS+=(--exclude "${pattern}")
done

# --- assemble the lychee invocation ----------------------------------------
LYCHEE_ARGS=(--no-progress --root-dir "${ROOT_ABS}")

if [[ "${LINKCHECK_CHECK_EXTERN}" == "1" ]]; then
  echo "External link checking: ENABLED"
  LYCHEE_ARGS+=(
    --timeout "${LINKCHECK_TIMEOUT}"
    --max-retries "${LINKCHECK_MAX_RETRIES}"
    --retry-wait-time "${LINKCHECK_RETRY_WAIT}"
    --max-concurrency "${LINKCHECK_MAX_CONCURRENCY}"
  )
  if [[ "${LINKCHECK_ACCEPT_TIMEOUTS}" == "1" ]]; then
    LYCHEE_ARGS+=(--accept-timeouts)
  fi
  if [[ -n "${LINKCHECK_ACCEPT}" ]]; then
    LYCHEE_ARGS+=(--accept "${LINKCHECK_ACCEPT}")
  fi
  if [[ "${LINKCHECK_CACHE}" == "1" ]]; then
    LYCHEE_ARGS+=(--cache --max-cache-age "${LINKCHECK_MAX_CACHE_AGE}")
  fi
  # Use the GitHub API for github.com links to avoid anonymous rate limiting.
  GH_API_TOKEN="${GITHUB_TOKEN:-${GH_TOKEN:-}}"
  if [[ -n "${GH_API_TOKEN}" ]]; then
    LYCHEE_ARGS+=(--github-token "${GH_API_TOKEN}")
  fi
else
  echo "External link checking: DISABLED (internal links only)"
  # Offline: only local files are checked; external URLs are skipped.
  LYCHEE_ARGS+=(--offline)
fi

# Anchor/fragment checking is opt-in (see LINKCHECK_FRAGMENTS above).
if [[ "${LINKCHECK_FRAGMENTS}" == "1" ]]; then
  LYCHEE_ARGS+=(--include-fragments)
fi

# Write a markdown report for the console (on failure) and, if requested, for
# use as a GitHub issue body.
REPORT="${LINKCHECK_FAILURES_OUT:-$(mktemp)}"
LYCHEE_ARGS+=(--format markdown --output "${REPORT}")

# Run from a stable cwd (not the docs root) so any on-disk cache (.lycheecache)
# lands somewhere cacheable rather than in the ephemeral build output.
echo "Running lychee over ${ROOT_ABS}"
set +e
"${LYCHEE_BIN}" "${LYCHEE_ARGS[@]}" "${EXCLUDE_ARGS[@]}" "${ROOT_ABS}/**/*.html"
rc=$?
set -e

if [[ ${rc} -ne 0 ]]; then
  echo "=== lychee report ==="
  cat "${REPORT}" || true
  echo "====================="
  echo "Link check failed (lychee exit ${rc})."
else
  echo "All checked links OK."
fi

exit ${rc}
