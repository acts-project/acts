import pandas
import argparse
import logging
import pathlib


log = logging.getLogger("event_trimmer")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input",
        type=pathlib.Path,
        help="input event directory",
    )

    parser.add_argument(
        "output",
        type=pathlib.Path,
        help="output event directory",
    )

    parser.add_argument(
        "-i", "--event-id", help="event ID in input directory", default=0, type=int
    )

    parser.add_argument(
        "-p",
        "--particle-id",
        help="particle ID to filter",
        type=int,
        required=True,
        action="append",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    to_keep = args.particle_id

    log.info(
        "Keeping %d particles: %s", len(to_keep), ", ".join(str(x) for x in to_keep)
    )

    origin_event_prefix = "event%09d-" % args.event_id
    destination_event_prefix = "event%09d-" % 0

    # Logic for processing the particle initial states
    origin_particles_initial_file = args.input / (
        origin_event_prefix + "particles_initial.csv"
    )
    particles_initial_df = pandas.read_csv(origin_particles_initial_file)
    log.info(
        "Read data for %d initial input particles from %s",
        particles_initial_df.shape[0],
        origin_particles_initial_file,
    )
    filtered_particles_initial_df = particles_initial_df[
        particles_initial_df["particle_id"].isin(to_keep)
    ]
    destination_particles_initial_file = args.output / (
        destination_event_prefix + "particles_initial.csv"
    )
    filtered_particles_initial_df.to_csv(
        destination_particles_initial_file, index=False
    )
    log.info(
        "Wrote data for %d initial output particles to %s",
        filtered_particles_initial_df.shape[0],
        destination_particles_initial_file,
    )

    # Logic for processing the particle final states
    origin_particles_final_file = args.input / (
        origin_event_prefix + "particles_final.csv"
    )
    particles_final_df = pandas.read_csv(origin_particles_final_file)
    log.info(
        "Read data for %d final input particles from %s",
        particles_final_df.shape[0],
        origin_particles_final_file,
    )
    filtered_particles_final_df = particles_final_df[
        particles_final_df["particle_id"].isin(to_keep)
    ]
    destination_particles_final_file = args.output / (
        destination_event_prefix + "particles_final.csv"
    )
    filtered_particles_final_df.to_csv(destination_particles_final_file, index=False)
    log.info(
        "Wrote data for %d final output particles to %s",
        filtered_particles_final_df.shape[0],
        destination_particles_final_file,
    )

    # Logic for processing hits
    origin_hits_file = args.input / (origin_event_prefix + "hits.csv")
    hits_df = pandas.read_csv(origin_hits_file)
    log.info("Read data for %d input hits from %s", hits_df.shape[0], origin_hits_file)
    hits_filtered_df = hits_df[hits_df["particle_id"].isin(to_keep)]
    destination_hits_file = args.output / (destination_event_prefix + "hits.csv")
    hits_filtered_df.to_csv(destination_hits_file, index=False)
    log.info(
        "Wrote data for %d output hits to %s",
        hits_filtered_df.shape[0],
        destination_hits_file,
    )

    # Logic for processing measurements
    origin_measurements_file = args.input / (origin_event_prefix + "measurements.csv")
    measurements_df = pandas.read_csv(origin_measurements_file)
    log.info(
        "Read data for %d input measurements from %s",
        measurements_df.shape[0],
        origin_measurements_file,
    )
    measurements_filtered_df = measurements_df[hits_df["particle_id"].isin(to_keep)]
    measurement_ids = list(measurements_filtered_df.index)
    meas_id_map = {a: b for (b, a) in enumerate(measurement_ids)}
    measurements_df["measurement_id"] = measurements_df["measurement_id"].apply(
        lambda x: meas_id_map.get(x, -1)
    )
    measurements_filtered_df = measurements_df[hits_df["particle_id"].isin(to_keep)]
    destination_measurements_file = args.output / (
        destination_event_prefix + "measurements.csv"
    )
    measurements_filtered_df.to_csv(destination_measurements_file, index=False)
    log.info(
        "Wrote data for %d output measurements to %s",
        measurements_filtered_df.shape[0],
        destination_measurements_file,
    )

    # Logic for building the simhit map
    new_df = pandas.DataFrame(
        {
            "measurement_id": list(range(measurements_filtered_df.shape[0])),
            "hit_id": list(range(measurements_filtered_df.shape[0])),
        }
    )
    destination_simhit_map_file = args.output / (
        destination_event_prefix + "measurement-simhit-map.csv"
    )
    new_df.to_csv(destination_simhit_map_file, index=False)
    log.info(
        "Wrote data for %d output measurement-to-hit mappings to %s",
        new_df.shape[0],
        destination_simhit_map_file,
    )

    # Logic for processing cells
    origin_cells_file = args.input / (origin_event_prefix + "cells.csv")
    cells_df = pandas.read_csv(origin_cells_file)
    log.info(
        "Read data for %d input cells from %s", cells_df.shape[0], origin_cells_file
    )
    filter = cells_df["measurement_id"].isin(measurement_ids)
    cells_df["measurement_id"] = cells_df["measurement_id"].apply(
        lambda x: meas_id_map.get(x, -1)
    )
    cells_filtered_df = cells_df[filter]
    destination_cells_file = args.output / (destination_event_prefix + "cells.csv")
    cells_filtered_df.to_csv(destination_cells_file, index=False)
    log.info(
        "Wrote data for %d output cells to %s",
        cells_filtered_df.shape[0],
        destination_cells_file,
    )


if __name__ == "__main__":
    main()
