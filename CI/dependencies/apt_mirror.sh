#!/bin/bash
# Point apt at a mirror near the runner, when the runner fleet names one.
#
# archive.ubuntu.com is reachable but slow from some networks: measured from CERN
# it served ~48 kB/s against ~1.9 MB/s for the Swiss mirror, which made an
# apt-get the longest step in jobs that were not otherwise network-bound.
#
# The mirror URL is NOT hardcoded here, because the right answer depends on where
# the job is running and this repository cannot know that. Self-hosted fleets
# export HUSK_UBUNTU_APT_MIRROR into the job container; GitHub-hosted runners do
# not set it, and then this script does nothing at all. That keeps the knowledge
# in the one place that has it, and keeps the hosted path provably unchanged.
#
# The UBUNTU in the name is load-bearing: a distribution mirror is only valid for
# the distribution whose suites it carries. Pointing a Debian image at an Ubuntu
# mirror does not degrade to slow, it 404s, because the suite names do not overlap
# (bookworm vs noble). The gate below enforces that rather than leaving it to the
# rewrite patterns to happen to miss.
#
# Safe to run more than once, on any image: it edits only lines that still point
# at the upstream archive.

set -u

if [ -z "${HUSK_UBUNTU_APT_MIRROR:-}" ]; then
  # Exiting mute is right for a hosted runner, which is never expected to name a
  # mirror. On a self-hosted slot it is a misconfiguration whose only other
  # symptom is an unexplained slow apt-get, so leave a trace to grep for.
  if [ "${RUNNER_ENVIRONMENT:-}" = "self-hosted" ]; then
    echo "apt mirror: HUSK_UBUNTU_APT_MIRROR unset on a self-hosted runner, using upstream archive"
  fi
  exit 0
fi

command -v apt-get > /dev/null 2>&1 || exit 0

# Debian and its derivatives all ship apt, so apt-get is not evidence of Ubuntu.
# Absent os-release, skip: an unidentifiable image is not one to rewrite blind.
distro=""
if [ -r /etc/os-release ]; then
  # shellcheck disable=SC1091
  distro=$(. /etc/os-release && echo "${ID:-}")
fi
if [ "$distro" != "ubuntu" ]; then
  echo "apt mirror: skipping, distro is '${distro:-unknown}', not ubuntu"
  exit 0
fi

if command -v sudo > /dev/null 2>&1; then SUDO="sudo"; else SUDO=""; fi

# 24.04 uses deb822 (sources.list.d/*.sources); older images use sources.list.
# Rewrite whichever exist, and both hosts: the country mirrors carry -security
# too, so leaving security.ubuntu.com alone would leave half the fetch slow.
#
# ports.ubuntu.com (arm64 and friends) is deliberately NOT rewritten: the country
# mirrors that carry the x86 archive frequently do not carry ports, so redirecting
# it would break rather than speed up a non-x86 job.
changed=0
for f in /etc/apt/sources.list /etc/apt/sources.list.d/*.sources \
  /etc/apt/sources.list.d/*.list; do
  [ -f "$f" ] || continue
  # One expression per scheme rather than `https\?://`: that is a GNU extension,
  # and silently matching nothing is the worst failure mode here — the job would
  # simply stay slow with no error to notice.
  before=$($SUDO cat "$f")
  $SUDO sed -i \
    -e "s|http://archive\.ubuntu\.com/ubuntu|${HUSK_UBUNTU_APT_MIRROR}|g" \
    -e "s|https://archive\.ubuntu\.com/ubuntu|${HUSK_UBUNTU_APT_MIRROR}|g" \
    -e "s|http://security\.ubuntu\.com/ubuntu|${HUSK_UBUNTU_APT_MIRROR}|g" \
    -e "s|https://security\.ubuntu\.com/ubuntu|${HUSK_UBUNTU_APT_MIRROR}|g" \
    "$f"
  [ "$before" = "$($SUDO cat "$f")" ] || changed=$((changed + 1))
done

# Report what actually happened. An unconditional success line would hide the
# case this script most needs to surface: the mirror was named, the image is
# ubuntu, and yet nothing matched — meaning the sources are in a form we do not
# recognise and the job is silently still on the slow path.
if [ "$changed" -gt 0 ]; then
  echo "apt mirror set to ${HUSK_UBUNTU_APT_MIRROR} (${changed} source file(s) rewritten)"
else
  echo "apt mirror: ${HUSK_UBUNTU_APT_MIRROR} named, but no upstream archive URLs found to rewrite"
fi
