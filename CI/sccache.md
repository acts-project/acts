# CERN S3 compiler cache

`Builds / linux_ubuntu` uses the [sccache action](../.github/actions/sccache/action.yml)
(Linux x64). Only `main` pushes publish; branches, PRs, and merge queue builds
read anonymously.

GitHub environments:

- `s3-sccache`: restrict to branch `main`; set `S3_ACCESS_KEY_ID` and
  `S3_SECRET_ACCESS_KEY`.
- `s3-sccache-read`: unrestricted, no secrets required.

Storage: `https://s3.cern.ch`, region `cern`, bucket `cache`, prefix
`acts-sccache/linux_ubuntu/v1/`. Objects expire 30 days after writing;
reads do not refresh expiry.

The workflow selects the CMake launcher and stops the daemon in an `always()`
step, recording statistics in both the job log and summary. In sccache 0.18.0,
read-only misses count as `Cache write errors` even though no upload is attempted.
