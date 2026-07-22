#!/bin/bash
# Point apt at a mirror near the runner, when the runner fleet names one.
#
# archive.ubuntu.com is reachable but slow from some networks: measured from CERN
# it served ~48 kB/s against ~1.9 MB/s for the Swiss mirror, which made an
# apt-get the longest step in jobs that were not otherwise network-bound.
#
# The mirror URL is NOT hardcoded here, because the right answer depends on where
# the job is running and this repository cannot know that. Self-hosted fleets
# export HUSK_APT_MIRROR into the job container; GitHub-hosted runners do not set
# it, and then this script does nothing at all. That keeps the knowledge in the
# one place that has it, and keeps the hosted path provably unchanged.
#
# Safe to run more than once, on any image: it edits only lines that still point
# at the upstream archive.

set -u

[ -n "${HUSK_APT_MIRROR:-}" ] || exit 0
command -v apt-get > /dev/null 2>&1 || exit 0

if command -v sudo > /dev/null 2>&1; then SUDO="sudo"; else SUDO=""; fi

# 24.04 uses deb822 (sources.list.d/*.sources); older images use sources.list.
# Rewrite whichever exist, and both hosts: the country mirrors carry -security
# too, so leaving security.ubuntu.com alone would leave half the fetch slow.
for f in /etc/apt/sources.list /etc/apt/sources.list.d/*.sources \
  /etc/apt/sources.list.d/*.list; do
  [ -f "$f" ] || continue
  # One expression per scheme rather than `https\?://`: that is a GNU extension,
  # and silently matching nothing is the worst failure mode here — the job would
  # simply stay slow with no error to notice.
  $SUDO sed -i \
    -e "s|http://archive\.ubuntu\.com/ubuntu|${HUSK_APT_MIRROR}|g" \
    -e "s|https://archive\.ubuntu\.com/ubuntu|${HUSK_APT_MIRROR}|g" \
    -e "s|http://security\.ubuntu\.com/ubuntu|${HUSK_APT_MIRROR}|g" \
    -e "s|https://security\.ubuntu\.com/ubuntu|${HUSK_APT_MIRROR}|g" \
    "$f"
done

echo "apt mirror set to ${HUSK_APT_MIRROR}"
