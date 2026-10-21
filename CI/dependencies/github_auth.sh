#!/bin/bash
# Make git authenticate its github.com reads with GITHUB_TOKEN, so they count
# against the token's quota instead of the runners' shared anonymous one.
set -e
set -u

if [ -z "${GITHUB_TOKEN:-}" ]; then
  echo "GITHUB_TOKEN not set: github.com fetches stay anonymous"
  exit 0
fi

# An Authorization header, not credentials in the URL: git does not send
# credentials preemptively. It issues the request anonymously and only re-sends
# it authenticated after a 401 challenge -- and github.com never challenges for
# a *public* repo, it just answers. So the obvious-looking
# `url.https://x-access-token:${GITHUB_TOKEN}@github.com/.insteadOf` leaves
# every read of a public repo anonymous (verifiable: a deliberately bogus token
# in the URL still clones acts-project/ci-dependencies just fine), which is
# precisely the traffic GitHub sheds with "GitHub is temporarily limiting some
# unauthenticated downloads". A header goes out on the first request, so it
# actually authenticates. This is what actions/checkout does for the workspace
# checkout; here we do it globally, for the clones spack makes on its own.
#
# The config key is URL-matched, so it reaches https://github.com/ only -- no
# other host, and not the ssh form. --replace-all keeps it idempotent: extra
# headers accumulate, and a cached spack install can run this more than once.
auth=$(printf '%s' "x-access-token:${GITHUB_TOKEN}" | base64 | tr -d '\n')
if [ -n "${GITHUB_ACTIONS:-}" ]; then
  echo "::add-mask::${auth}"
fi
git config --global --replace-all \
  "http.https://github.com/.extraheader" \
  "AUTHORIZATION: basic ${auth}"
echo "Authenticating https://github.com/ fetches with GITHUB_TOKEN"

# Silently degrading to anonymous is the failure mode this whole file exists to
# prevent, so check rather than assume. The anonymous REST quota is 60/h; any
# accepted token is far above it. Warn only: a hiccup on this probe is no
# reason to fail the job, and the fetches themselves are retried anyway.
rate_limit=$(curl -sS -o /dev/null -D - -m 30 \
  -H "Authorization: basic ${auth}" \
  https://api.github.com/rate_limit 2>/dev/null |
  tr -d '\r' | awk 'tolower($1) == "x-ratelimit-limit:" { print $2 }') || true
if [ -z "${rate_limit}" ]; then
  echo "Warning: could not probe the REST rate limit; assuming the token works"
elif [ "${rate_limit}" -le 60 ]; then
  echo "Warning: GITHUB_TOKEN was not accepted (rate limit ${rate_limit}/h);" \
    "github.com fetches will be treated as anonymous"
else
  echo "GITHUB_TOKEN accepted (rate limit ${rate_limit}/h)"
fi
