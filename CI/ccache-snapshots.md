# CERN S3 compiler caches

Builds, Analysis, Detray CUDA, and PyPI wheels restore local ccache snapshots,
including direct-mode manifests. Restores are anonymous. Main pushes publish;
PyPI additionally permits scheduled/manual nightly builds of main. Explicit
release-ref wheel builds are read-only. Writers use the protected `s3-sccache`
environment with `S3_ACCESS_KEY_ID` and `S3_SECRET_ACCESS_KEY`; readers need no
environment.

Storage: `https://s3.cern.ch`, region `cern`, bucket `cache`. Prefixes under
`acts-sccache/ccache-snapshots/<repository>/` separate workflows, jobs, platforms,
and matrix variants. The original `linux_ubuntu/v1/` prefix is unchanged.
Each new variant starts cold until its first main publication.

Restore and publish remain explicit so jobs retain saves after failures and
before later tests. LCG nightly only publishes after a successful configure.
Cleanup and statistics run where ccache is available: the job container, sourced
LCG/Key4hep view, EIC Docker container, or wheel build environment. Wheels clean
before tests; failures before that point rely on ccache's automatic size cleanup.
Cache limits remain 500 MB, except Analysis's existing 1 GB limit.

Publication uploads an immutable tar, then updates `latest.json`. There is no
extra concurrency group; overlapping uploads may leave an older valid snapshot
as latest. Restores verify checksums and extract into staging. Transfer failures
are nonfatal. Entries are already compressed, so the tar is uncompressed. The
existing 30-day S3 lifecycle applies to snapshots and pointers. The old GitHub
cache retention workflow remains available for manual cleanup only.

Compare GitHub's restore, build, and publish step durations and ccache statistics.
Do not set `CCACHE_BASEDIR`: path rewriting breaks source-path FPE masks.
