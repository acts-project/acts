import argparse
import pathlib
import tempfile
import logging
import os
import traccc_physics_plots.run


log = logging.getLogger("traccc_physics_plots.run_seeding_example")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    log.info("Starting data gathering for traccc physics performance plots")

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "exec",
        type=pathlib.Path,
        help="the traccc CUDA seeding example executable",
    )

    parser.add_argument(
        "outdir",
        type=pathlib.Path,
        help="the directory in which to write output files",
    )

    parser.add_argument(
        "--root-to-csv",
        help="the ROOT-to-CSV executable",
        type=pathlib.Path,
        default=None,
        dest="root_to_csv",
    )

    args = parser.parse_args()

    if args.root_to_csv is not None:
        log.info("Using provided ROOT-to-CSV converter at %s", str(args.root_to_csv))
        root_to_csv = args.root_to_csv
    else:
        root_to_csv = (
            pathlib.Path(os.path.abspath(__file__)).parent.parent
            / "root_to_csv"
            / "root_to_csv"
        )
        log.info("No ROOT-to-CSV converter provided; assuming %s", str(root_to_csv))

    assert root_to_csv.is_file()

    with tempfile.TemporaryDirectory() as tmpdirname:
        log.info("Running seeding example in %s", tmpdirname)
        traccc_physics_plots.run.run_convert_seeding_example(
            args.exec,
            traccc_physics_plots.run.SEEDING_EXAMPLE_ARGS,
            tmpdirname,
            args.outdir,
            root_to_csv,
        )

    log.info("Data gathering complete")


if __name__ == "__main__":
    main()
