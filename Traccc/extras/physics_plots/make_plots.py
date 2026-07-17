import argparse
import pathlib
import logging
import traccc_physics_plots.plot


log = logging.getLogger("traccc_physics_plots.run_seeding_example")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    log.info("Starting traccc physics performance plotting")

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="data input directories, with optional names",
        nargs=2,
        action="append",
        dest="input",
        metavar=("directory", "name"),
    )

    parser.add_argument(
        "outdir",
        type=pathlib.Path,
        help="the directory in which to write output files",
    )

    args = parser.parse_args()

    plot_candidate_names = {
        i: (pathlib.Path(l[0]), l[1]) for i, l in enumerate(args.input)
    }

    traccc_physics_plots.plot.make_plots(plot_candidate_names, args.outdir)

    log.info("Data plotting complete")


if __name__ == "__main__":
    main()
