# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import matplotlib.pyplot
import pandas
import argparse
import logging
import pathlib
from datetime import datetime

log = logging.getLogger("traccc_plot")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input",
        type=pathlib.Path,
        help="input CSV database",
    )

    parser.add_argument(
        "output",
        type=pathlib.Path,
        help="output plot file",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        help="value in seconds under which to hide kernels",
        default=0.001,
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    log.info('Loading CSV database "%s"', args.input)
    df = pandas.read_csv(args.input)

    valid_kernels = set()

    ccs = []

    log.info("Filtering kernels with a threshold of %f seconds", args.threshold)

    for x in df.iloc:
        if x["rec_throughput"] >= args.threshold:
            valid_kernels.add(x["kernel"])
        if x["cc"] not in ccs:
            ccs.append(x["cc"])

    log.info(
        "Plotting %d kernels over %d Compute Capabilities", len(valid_kernels), len(ccs)
    )

    sorted_valid_kernels = sorted(list(valid_kernels))

    px = 1 / matplotlib.pyplot.rcParams["figure.dpi"]

    f, a = matplotlib.pyplot.subplots(1, 1, figsize=(806 * px, 600 * px))

    xrange = list(range(len(ccs)))

    for x in sorted_valid_kernels:
        x_data = []
        y_data = []
        y_error = []

        for y in df.iloc:
            if y["kernel"] == x:
                x_data.append(ccs.index(y["cc"]))
                y_data.append(y["rec_throughput"] * 1000)
                y_error.append(y["rec_throughput_dev"] * 1000)

        a.errorbar(
            x_data, y_data, yerr=y_error, label=x, marker="x", capsize=2, elinewidth=1
        )

    a.set_xticks(
        xrange,
        ccs,
    )
    a.set_xlabel("Compute Capability (PTX Target)")
    a.set_ylabel("Reciprocal throughput [ms]")
    a.legend(prop={"family": "monospace"})
    a.grid(color="lightgray", linewidth=0.5)
    f.tight_layout()

    log.info('Saving plot to "%s"', args.output)

    f.savefig(args.output)


if __name__ == "__main__":
    main()
