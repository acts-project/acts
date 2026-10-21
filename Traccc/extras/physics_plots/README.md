# Physics plot generators

This directory contains code to produce traccc physics plots. To gather data,
run the `run_seeding_example.py` script with a compiled CUDA seeding example:

```bash
$ python run_seeding_example.py .../build/bin/traccc_seeding_example_cuda output_directory
```

This will produce CSV data in the directory `output_directory`. Use this data
as input to the plot generator:

```bash
$ python make_plots.py -i output_directory "Current commit" plot_directory
```

This will produce a series of plots in `plot_directory`. To make comparison
plots, provide multiple input directories, such as:

```bash
$ python make_plots.py -i output_directory1 "Current commit"  -i output_directory2 "Previous commit" plot_directory
```
