# CERN S3 ccache snapshot pilot

`Builds / linux_ubuntu` restores a local ccache directory from S3, including
its direct-mode manifests. One action call restores the cache and registers a
post-job upload. Only successful `main` push jobs publish; other runs read anonymously.

The pilot reuses the `s3-sccache` environment and its `S3_ACCESS_KEY_ID` /
`S3_SECRET_ACCESS_KEY` secrets for main pushes. Other runs need no environment.
Storage: `https://s3.cern.ch`, region `cern`, bucket `cache`, prefix
`acts-sccache/ccache-snapshots/<repository>/linux_ubuntu/v1/`.

Main runs `ccache --cleanup` in the build environment, then uploads an immutable
tar of the 500 MB cache and updates `latest.json`. The transfer helper only needs
access to the cache directory, not a ccache installation. Cleanup and statistics
stay in the build environment. If cleanup fails, set `CCACHE_SNAPSHOT_SKIP_SAVE=true`
via `GITHUB_ENV` to skip the upload without failing the job.
Publication adds no concurrency group. Overlapping uploads may leave an older
valid snapshot as latest; ccache still validates entries when using them. Restores
verify the checksum and extract into staging. Failed transfers are nonfatal.
The tar is uncompressed since ccache already compresses entries. The existing
30-day lifecycle applies to both snapshots and the pointer.

The first main run seeds this namespace; earlier PR runs are cold. Compare
GitHub's setup, build, and post-job publish durations. Cache statistics appear in
both the log and summary.
