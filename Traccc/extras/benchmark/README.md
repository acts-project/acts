# Benchmark tools

This directory contains some tools to help benchmark traccc:

* `benchmark.py` is designed to benchmark and compare the performance of traccc over time; it examines performance commit-by-commit and produces a CSV output.
* `plot.py` can be used to plot the output of `benchmark.py`.
* `parse_profile.py` takes an NSight Compute profile captured with `--section LaunchStats --section Occupancy --metrics gpu__time_duration.sum` and converts it into an easy-to-process CSV format.

All tools are designed to run with uv, e.g.:

```bash
$ uv run python3 benchmark.py ...
```
