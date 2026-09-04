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
acts_analysis_residuals_and_pulls.root
ckf_seeded.html
ckf_seeded_plots
ckf_truth_estimated.html
ckf_truth_estimated_plots
ckf_truth_smeared.html
ckf_truth_smeared_plots
ivf_seeded.html
ivf_seeded_plots
ivf_truth_estimated.html
ivf_truth_estimated_plots
ivf_truth_smeared.html
ivf_truth_smeared_plots
performance_ckf_seeded.root
performance_ckf_truth_estimated.root
performance_ckf_truth_smeared.root
performance_truth_tracking.root
performance_vertexing_seeded.root
performance_vertexing_seeded_hist.root
performance_vertexing_truth_estimated.root
performance_vertexing_truth_estimated_hist.root
performance_vertexing_truth_smeared.root
performance_vertexing_truth_smeared_hist.root
run.log
truth_tracking.html
truth_tracking_plots
```

The `.root` files are the performance output files and corresponding histogram
files. The physmon job log file is also available. Finally, the `_plots` folder
contain plots of all the histogram comparisons, and the `.html` files contain
single-file reports showing the detailed results. An example of an HTML report
looks like this:

![Physmon report](physmon/physmon_report.png) {width=100%}

If you get a physmon job failure on your pull request, please investigate the
failing report(s), and try to understand if the change causing the
discrepancies is expected.

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

`CI/physmon/phys_perf_mon.sh` populates `CI/physmon/reference/` from the manifest
before running any comparison, so local runs need no extra step. Downloaded blobs
are cached under `~/.cache/acts/physmon-references` (override with
`ACTS_PHYSMON_CACHE`), so switching between branches with different references
does not re-download everything.

To manage the directory by hand:

```console
$ uv run --no-project CI/physmon/reference.py pull     # populate CI/physmon/reference/
$ uv run --no-project CI/physmon/reference.py verify   # check it against the manifest
```

Set `ACTS_PHYSMON_NO_FETCH=1` to stop `phys_perf_mon.sh` from touching the
reference directory, if you want to point it at files of your own. The manifest
is in `sha256sum` format, so `sha256sum -c ../reference.sha256` works from inside
the reference directory too.

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
