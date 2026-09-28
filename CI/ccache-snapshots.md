# CERN S3 ccache snapshot pilot

`Builds / linux_ubuntu` restores a local ccache directory from S3, including
its direct-mode manifests. Only `main` pushes publish; other runs read anonymously.

The pilot reuses the `s3-sccache` environment and its `S3_ACCESS_KEY_ID` /
`S3_SECRET_ACCESS_KEY` secrets. Other runs use `s3-sccache-read`.
Storage: `https://s3.cern.ch`, region `cern`, bucket `cache`, prefix
`acts-sccache/ccache-snapshots/<repository>/linux_ubuntu/v1/`.

Main uploads an immutable tar of the 500 MB cache, then updates `latest.json`.
Publication adds no concurrency group. Overlapping uploads may leave an older
valid snapshot as latest; ccache still validates entries when using them. Restores
verify the checksum and extract into staging. Failed transfers are nonfatal.
The tar is uncompressed since ccache already compresses entries. The existing
30-day lifecycle applies to both snapshots and the pointer.

The first main run seeds this namespace; earlier PR runs are cold. Compare
GitHub's restore, build, and publish step durations. Cache statistics appear in
both the log and summary.
