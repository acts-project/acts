@page physmon Physics Monitoring in CI
@brief Monitoring physics performance in the ACTS CI

The ACTS CI runs a suite of physics performance monitoring jobs dubbed
*physmon*. The purpose is to monitor and detect changes in the physics
performance, both intentional and accidental.

The associated job will run a number of workflow combinations. Currently, this
includes the truth tracking and OpenDataDetector *full chain* workflows. The
latter is further split into configurations with full seeding, truth smeared or
truth estimated seeds. These jobs produce performance output files.

These performance output files are then evaluated using the programs dubbed
[](analysis_apps). After this step, the job will then run comparisons of the
diagnostics histograms created by these programs. If the comparisons indicate
that the histogram contents changed meaningfully, the job will fail and report
this on a pull request.

## Histogram comparison configuration

Histcmp uses the TOML files in `CI/physmon/config`. `standard.toml` preserves
histcmp's default checks for comparisons without a workflow-specific config;
`default.toml` and the track-fitting configs specify the existing physmon
thresholds and exceptions. Use `disabled = true` to disable an inherited check.
`info_only.toml` defines no checks for the visual GX2F/KF comparison.

The `[report.tags]` tables group histograms in the HTML reports by physics
quantity or diagnostic, such as efficiency, fakes, residuals, pulls, resolution,
and bias. Patterns are case-sensitive shell globs matched against the full
histogram key, and a histogram can have multiple tags. For example,
`pullwidth_d0_vs_eta` belongs to Pulls, Resolution, and Kinematics. Update the
relevant tag tables when adding new histogram families.

The YAML files in the same directory configure `generic_plotter.py` binning
and derived histograms; they are separate from the histcmp configuration.

## How do I investigate a physmon failure?

The physmon CI job attaches its results as an artifact to the CI run (also for successful runs)
From your pull request, you need to click on the *Checks* tab at the top:

![Checks button on GitHub](physmon/physmon_checks.png)

From there, click on the *Builds* workflow on the left:

![Builds runs the physmon job](physmon/physmon_run.png)

On the workflow overview, scroll down to find the attached artifacts, and locate the *physmon* artifact. You can click to download it:

![physmon artifact](physmon/physmon_artifact.png) {width=100%}

> [!note]
> GitHub Actions artifacts are deleted after a fixed amount of time. If the CI
> completed a while ago, it is possible this report is no longer available.

After the download, you need to unzip the archive, whose contents will look similar to this:

```
summary.md / summary.html     the results table, as on the pull request
data/<workflow>/...           performance output and histogram files
html/<workflow>/...           one report per comparison, with its plots
logs/                         one log per job
benchmarks/                   run time and peak memory of every job
reference-candidate.sha256    the manifest if all failing comparisons were accepted
```

Workflows with seeding variants have one directory per variant, for example
`data/trackfinding_1muon/GridTriplet`. The summary groups the comparisons by
algorithm and sample, with one row per comparison and one column per step
(producer, plot, histcmp), each linking to its log, so it shows which step
failed first. Its *Memory* section lists every job that used more memory than
it reserves. The `.html` files are single-file reports showing the detailed
results.
An example of an HTML report looks like this:

![Physmon report](physmon/physmon_report.png) {width=100%}

If you get a physmon job failure on your pull request, please investigate the
failing report(s), and try to understand if the change causing the
discrepancies is expected.

## How do I run physmon locally?

physmon is a [Snakemake](https://snakemake.readthedocs.io) workflow. It expects
the ACTS build in `build/` in the repository and activates it for every job
itself, so the shell needs no setup, and Snakemake can be installed with any
Python, for example with `uv tool install 'snakemake>=9,<10'`. From the
repository root:

```console
$ snakemake -s CI/physmon/Snakefile --cores 8 \
    --resources mem_mb=12000 --keep-incomplete all
```

This is the command the CI runs. `mem_mb` is the memory in MB that all jobs
running at the same time may reserve together; the CI uses 12 GB of its 16 GB.
The results are written to `physmon/`, starting with `physmon/summary.md`.

Instead of `all`, a target can be one group of comparisons (`simulation`, `kf`,
`gsf`, `gx2f`, `refit_kf`, `refit_gsf`, `fullchains`, `gx2f_vs_kf`) or a single
job such as `produce_trackfinding_1muon_Gbts`.

> [!note]
> A failed job does not stop the workflow: it leaves empty files in place of its
> outputs, so that the other jobs still run and the summary shows what failed.
> Snakemake then considers that job done. After fixing the cause, delete
> `physmon/` or add `--forcerun <job>` to run it again.

## How do I add a workflow or a seeding variant?

The catalog `CI/physmon/config/physmon.yaml` lists what runs; its header
describes every key.

- A new workflow needs a script `CI/physmon/workflows/physmon_<name>.py`, a
  block under `producers:` with the memory it reserves (`mem_mb`), and its
  comparisons under `comparisons:`.
- A new seeding variant of a track finding workflow needs its name in the
  script's `VARIANTS` list, the settings only it needs in the script's
  `seedingSettings()`, and a line under the workflow's `variants:` in the
  catalog.

The Snakemake rules are generated from the catalog. A new comparison has no
reference yet: it is reported as failed until its reference is published, see
[below](#how-do-i-update-the-reference-files).

## Where the reference files live

The reference histograms are **not** committed to the repository. They are
stored as content-addressed blobs in an OCI registry, and the repository records
only their hashes in `CI/physmon/reference.sha256`:

```
acf50e815ae89cb3053ffc04c6099ac877b78c7fe32cfd5ee8ae890abe56e28b  simulation/particles_fatras_hist.root
039fb280f421793ac92063ad96594c45d1900c91db91385728d8f5323aebc7fc  simulation/particles_geant4_hist.root
...
```

Committing the histograms themselves was costing the repository several MB of
permanent history per update, since ROOT files do not compress or delta against
each other. A reference update is now a one-line change per file, and reviewers
can see exactly which references a pull request touches instead of a list of
changed binaries.

The workflow's `fetch_references` rule populates `CI/physmon/reference/` from the
manifest, so local runs need no extra step: every reference is a normal output in
the DAG, and a comparison that needs one waits for the fetch. Downloaded blobs are
cached under `~/.cache/acts/physmon-references` (override with
`ACTS_PHYSMON_CACHE`), so switching between branches with different references
does not re-download everything.

To manage the directory by hand:

```console
$ uv run --no-project CI/physmon/reference.py pull     # populate CI/physmon/reference/
$ uv run --no-project CI/physmon/reference.py verify   # check it against the manifest
```

Set `ACTS_PHYSMON_NO_FETCH=1` (or `--config fetch_references=0`) to drop the
`fetch_references` rule from the workflow, if you want to point it at files of
your own: the references are then plain inputs that must already be on disk. The
manifest is in `sha256sum` format, so `sha256sum -c ../reference.sha256` works
from inside the reference directory too.

## How do I update the reference files?

Reference updates are published by a maintainer; you do not need registry
credentials to propose one.

1. Confirm the physmon differences on your pull request are expected, and say so
   on the pull request so a maintainer can act on it.
2. A maintainer runs the **Upload physmon references** workflow from the `main`
   branch, passing your pull request number or URL. The workflow finds the
   *Builds* run holding the physmon outputs of the pull request's head commit,
   uploads the new histograms, and reports the manifest to commit. A *Builds*
   run id works too, if a specific run is wanted.
3. Commit that manifest as `CI/physmon/reference.sha256` and **push it** to your
   branch.

The manifest is also attached to every physmon run as
`reference-candidate.sha256` inside the `physmon` artifact, so you can see in
advance exactly which entries an update would change.

The same resolution is available locally, for a maintainer with write access to
the registry:

```console
$ uv run --no-project CI/physmon/reference.py update 5736      # or the pull request URL
```

> [!note]
> The references are taken from the pull request's **head commit**. If the
> branch is pushed to after a physmon run, that run no longer describes the pull
> request, and the upload waits for the new run rather than publishing stale
> histograms.

> [!note]
> Step 3 must be a push, not a *Re-run jobs* on the existing run. Re-running
> replays the same commit, which still references the old hashes. Pushing the
> manifest commit is what makes CI resolve the new blobs.

> [!note]
> The blobs have to exist before the manifest naming them is committed, so the
> upload in step 2 has to happen before step 3. If the order slips, the physmon
> job fails while fetching rather than doing anything silently wrong.

Only the comparisons that actually failed are updated. ROOT files embed a
creation timestamp, so every physmon run produces a different hash for *every*
file even when the physics is identical; the changed set therefore comes from the
histcmp results rather than from comparing hashes. This is what keeps an update
to one workflow from rewriting all forty entries.
