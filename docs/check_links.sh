#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LINKCHECK_IGNORE_REPO="${LINKCHECK_IGNORE_REPO:-acts-project/linkcheck-ignore}"
LINKCHECK_IGNORE_REF="${LINKCHECK_IGNORE_REF:-main}"
LINKCHECK_IGNORE_PATH="${LINKCHECK_IGNORE_PATH:-data.json}"
LINKCHECK_IGNORE_URL="${LINKCHECK_IGNORE_URL:-https://raw.githubusercontent.com/${LINKCHECK_IGNORE_REPO}/${LINKCHECK_IGNORE_REF}/${LINKCHECK_IGNORE_PATH}}"

# First-pass crawl settings.
LINKCHECK_THREADS="${LINKCHECK_THREADS:-4}"
LINKCHECK_TIMEOUT="${LINKCHECK_TIMEOUT:-20}"

# External link checking. Disable for fast, deterministic PR runs (internal
# links only); enable for the scheduled sweep that validates external URLs.
LINKCHECK_CHECK_EXTERN="${LINKCHECK_CHECK_EXTERN:-1}"

# Retry settings. Only the URLs that failed are retried, each re-checked in
# isolation (no crawl) with an escalated timeout and reduced concurrency, so a
# single slow or briefly-flaky external link does not fail the whole run. A
# link is only reported as broken if it fails the initial pass AND every retry.
LINKCHECK_RETRIES="${LINKCHECK_RETRIES:-0}"
LINKCHECK_RETRY_TIMEOUT="${LINKCHECK_RETRY_TIMEOUT:-60}"
LINKCHECK_RETRY_THREADS="${LINKCHECK_RETRY_THREADS:-1}"
LINKCHECK_RETRY_DELAY="${LINKCHECK_RETRY_DELAY:-15}"

# Optional: write the list of links that failed every attempt to this file (for
# reporting, e.g. an issue body in the scheduled workflow).
LINKCHECK_FAILURES_OUT="${LINKCHECK_FAILURES_OUT:-}"

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

EXTERN_ARGS=()
if [[ "${LINKCHECK_CHECK_EXTERN}" == "1" ]]; then
  EXTERN_ARGS+=(--check-extern)
  echo "External link checking: ENABLED"
else
  echo "External link checking: DISABLED (internal links only)"
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

# Extract the checked URLs from a linkchecker 'failures' output file. Each line
# has the form:  <count> "('<parent>', '<url>')"
extract_failed_urls() {
  local file="$1"
  [[ -f "${file}" ]] || return 0
  python3 - "${file}" <<'PY'
import ast, sys
seen = []
for line in open(sys.argv[1], encoding="utf-8"):
    line = line.strip()
    if not line:
        continue
    rest = line[line.find(" ") + 1:].strip()
    if rest.startswith('"') and rest.endswith('"'):
        rest = rest[1:-1]
    try:
        tup = ast.literal_eval(rest)
        url = tup[1] if isinstance(tup, tuple) and len(tup) > 1 else str(tup)
    except Exception:
        continue
    if url not in seen:
        seen.append(url)
for u in seen:
    print(u)
PY
}

# Run linkchecker. Args: <failures_file> <timeout> <threads> <recursion|""> <urls...>
run_linkchecker() {
  local failures_file="$1"; shift
  local timeout="$1"; shift
  local threads="$1"; shift
  local recursion="$1"; shift
  rm -f "${failures_file}"
  local rec_args=()
  if [[ -n "${recursion}" ]]; then
    rec_args=(--recursion-level="${recursion}")
  fi
  uvx --from LinkChecker linkchecker \
    --config="${SCRIPT_DIR}/linkcheckerrc" \
    "${EXTERN_ARGS[@]}" \
    --threads="${threads}" \
    --timeout="${timeout}" \
    "${rec_args[@]}" \
    --file-output="failures/utf-8/${failures_file}" \
    "${IGNORE_ARGS[@]}" \
    "$@"
}

FAILURES_FILE="${WORKDIR}/failures.txt"

echo "==> Initial link check (timeout=${LINKCHECK_TIMEOUT}s, threads=${LINKCHECK_THREADS})"
set +e
run_linkchecker "${FAILURES_FILE}" "${LINKCHECK_TIMEOUT}" "${LINKCHECK_THREADS}" "" "${START_URL}"
last_rc=$?
set -e

if [[ ${last_rc} -eq 0 ]]; then
  echo "All links OK."
  exit 0
fi

mapfile -t failed_urls < <(extract_failed_urls "${FAILURES_FILE}")
echo "Initial check reported ${#failed_urls[@]} failing link(s)."

if (( ${#failed_urls[@]} == 0 )); then
  echo "Link check failed but no specific URLs could be parsed; not retrying."
  exit 1
fi

round=0
while (( round < LINKCHECK_RETRIES )) && (( ${#failed_urls[@]} > 0 )); do
  round=$((round + 1))
  echo "==> Retry ${round}/${LINKCHECK_RETRIES} for ${#failed_urls[@]} link(s) after ${LINKCHECK_RETRY_DELAY}s (timeout=${LINKCHECK_RETRY_TIMEOUT}s, threads=${LINKCHECK_RETRY_THREADS})"
  printf '    %s\n' "${failed_urls[@]}"
  sleep "${LINKCHECK_RETRY_DELAY}"

  RETRY_FAILURES="${WORKDIR}/failures_retry_${round}.txt"
  set +e
  run_linkchecker "${RETRY_FAILURES}" "${LINKCHECK_RETRY_TIMEOUT}" "${LINKCHECK_RETRY_THREADS}" 0 "${failed_urls[@]}"
  last_rc=$?
  set -e

  if [[ ${last_rc} -eq 0 ]]; then
    echo "All previously-failing links recovered on retry ${round}."
    failed_urls=()
    break
  fi

  mapfile -t failed_urls < <(extract_failed_urls "${RETRY_FAILURES}")
  echo "Retry ${round} still has ${#failed_urls[@]} failing link(s)."
  if (( ${#failed_urls[@]} == 0 )); then
    echo "Retry ${round} failed but no specific URLs could be parsed; stopping."
    break
  fi
done

if [[ ${last_rc} -eq 0 ]]; then
  echo "All links OK after retries."
  exit 0
fi

echo "The following link(s) failed the initial check and every retry:"
printf '  %s\n' "${failed_urls[@]}"

if [[ -n "${LINKCHECK_FAILURES_OUT}" ]]; then
  printf '%s\n' "${failed_urls[@]}" > "${LINKCHECK_FAILURES_OUT}"
fi

exit 1
