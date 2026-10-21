# CCL data generator

This little Python tool is designed to generate example data for connected
component analysis and labeling problems. While the data generated from this
tool is not meaningful from a physics perspective, it exposes the most
important parameters from a computational performance point of view.

The configurable parameters are as follows:

* Number of detector modules
* Size of each detector module (assumed to be square)
* Mean hits per module
* Variance of hits per module
* Mean cells per hit
* Variance of cells per hit

## Generating test files

An example run of this program looks like this:

```
$ python3 ccl_generator.py -N 1000 -S 5000 --Mm 1.65 --Ms 0.95 --Hm 1.8 --Hs 0.89 -C 100 -o my_run
```

This will generate one hundred files, in the format `my_run_0000000000.csv` and
so on. The files should in a CSV format readably by traccc.
