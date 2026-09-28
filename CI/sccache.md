# CERN S3 compiler cache pilot

`Builds / linux_ubuntu` uses sccache 0.18.0 with the S3 backend. The composite
action `.github/actions/sccache` downloads the pinned Linux x64 binary and
verifies its SHA-256 digest. Other jobs
continue to use ccache and GitHub Actions cache storage.

## Setup action

The action installs, configures, and starts sccache. It currently supports Linux
x64, matching the pilot job. Its inputs are:

| Input | Default | Purpose |
| --- | --- | --- |
| `read-only` | `true` | Anonymous reads; set `false` to enable authenticated writes |
| `key-prefix` | Required | Prefix shared by compatible builds |
| `endpoint` | `https://s3.cern.ch` | S3 endpoint |
| `region` | `cern` | Signing region |
| `bucket` | `cache` | S3 bucket |
| `access-key-id` | Empty | Required for writes |
| `secret-access-key` | Empty | Required for writes |

The caller selects the GitHub environment and passes the `main`-only write
condition and credentials. The action exports non-secret `SCCACHE_*` settings
for subsequent steps, but does not select a CMake compiler launcher. Callers
must also include the final `always()` statistics/shutdown step, as shown in
`builds.yml`; the composite action does not stop the daemon automatically.

## GitHub setup

Configure the write environment before enabling the pilot:

- `s3-sccache`: allow only the **branch** `main` under deployment branches
  and tags. Store `S3_ACCESS_KEY_ID` and `S3_SECRET_ACCESS_KEY` here.
- `s3-sccache-read`: the job uses this unrestricted environment for branch,
  pull request, and merge queue builds. GitHub can create it automatically;
  it needs no S3 secrets or other setup.

Use a dedicated CI credential where possible, so it can be rotated or revoked
independently. A second OpenStack EC2 credential for the same user and project
does not by itself narrow permissions. Bucket or prefix restrictions require
storage-side authorization supported by the CERN service.

Only a `push` to `refs/heads/main` selects the write environment and enables
authenticated writes. All other runs use anonymous `READ_ONLY` access, including
pushes to feature branches, pull requests from forks, and merge queue builds.
Environment branch restrictions provide an additional server-side boundary;
do not put the write credentials in repository-wide secrets.

## Storage and retention

| Setting | Value |
| --- | --- |
| Endpoint | `https://s3.cern.ch` |
| Region used for signing | `cern` |
| Bucket | `cache` |
| Pilot object prefix | `acts-sccache/linux_ubuntu/v1/` |

The bucket must allow anonymous object reads for the pilot prefix and deny
anonymous writes. sccache writes individual cache entries during compilation;
there is no archive restore or upload step. Its cache format is separate from
ccache's, so the pilot starts cold.

The bucket lifecycle rule `acts-sccache-expire-after-30-days` expires objects
under `acts-sccache/` after 30 days. This is configured in CERN S3, not by the
workflow. Reads do not refresh object age. Expired entries become misses, and
only a subsequent `main` build can repopulate them. This is not a total-size cap.

Write credentials are passed only to the daemon startup step, not written to
configuration files. The daemon stays alive for the job and is stopped in an
`always()` cleanup step, which also publishes statistics in the job summary.

## Evaluate the pilot

1. Run a `main` build after merging the workflow change. Confirm the job summary
   reports `READ_WRITE`, an S3 cache location, and no cache write errors.
2. Run a PR build with overlapping source and the same compiler configuration.
   Confirm `READ_ONLY`, an S3 cache location, and cache hits. A PR before the
   first `main` writer will compile normally with a cold cache.
3. Compare build time and cache errors across a few runs before extending this
   to other jobs. Check S3 usage separately; sccache's local cache size setting
   does not limit remote storage.

In sccache 0.18.0, a read-only miss increments `Cache write errors`: its
read-only storage wrapper rejects the write locally, without sending an S3
upload. This is expected for readers; write errors on `main` need investigation.
The smoke test confirmed that read-only hits and misses leave S3 unchanged.

The shared CMake presets remain unchanged. The pilot overrides the C++ compiler
launcher on its configure command. To roll back, restore the job's ccache
restore/save steps and remove that override and its sccache setup/cleanup steps.
