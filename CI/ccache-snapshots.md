# CERN S3 ccache snapshot pilot

`Builds / linux_ubuntu` restores a complete local ccache directory from S3.
Compilation uses ccache's local direct-mode lookups, including its manifests.
The cache remains writable locally on every branch; only `main` pushes upload.

The pilot reuses the `s3-sccache` environment and its `S3_ACCESS_KEY_ID` /
`S3_SECRET_ACCESS_KEY` secrets. Other runs use `s3-sccache-read` and anonymous
reads. Storage is `https://s3.cern.ch`, region `cern`, bucket `cache`, under
`acts-sccache/ccache-snapshots/<repository>/linux_ubuntu/v1/`.

Each main run restores the latest snapshot, builds, cleans the cache to the
configured 500 MB limit, uploads an immutable tar, then updates `latest.json`.
Main publishers are serialized, and older runs cannot replace newer snapshots.
The tar is uncompressed because ccache already compresses its entries. Downloads
are checksum-verified and extracted into staging before becoming the local cache.
Missing or failed restores fall back to a cold build; failed uploads retain the
previous pointer. The existing 30-day `acts-sccache/` lifecycle expires snapshots
and the pointer after inactivity. No bucket configuration changes are required.

The first main run seeds this separate namespace; earlier PR runs will be cold.
Logs and job summaries include transfer, extraction, packing, build timings, and
ccache hit statistics. Compare restore + build + publication time with the
sccache pilot, as well as the warm build time itself. PRs omit publication.
