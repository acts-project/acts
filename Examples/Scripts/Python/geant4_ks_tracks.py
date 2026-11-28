"""Check K-short dipion masses with truth-matched CKF tracks from full_chain_odd.py.

The input CSV comes from geant4_ks_candidates.py. Generator IDs tag the decay;
only fitted perigee parameters enter the pair's closest approach and mass.
The helix extrapolation assumes full_chain_odd.py's uniform +2 T axial field.
"""

import argparse
from collections import Counter, defaultdict
import csv
import math
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

BARCODE = (
    "vertex_primary",
    "vertex_secondary",
    "particle",
    "generation",
    "sub_particle",
)
PION_MASS_GEV = 0.13957039


class Helix:
    def __init__(self, track):
        self.phi = track["ePHI_fit"]
        theta = track["eTHETA_fit"]
        d0, self.z0 = track["eLOC0_fit"], track["eLOC1_fit"]
        qop = track["eQOP_fit"]
        self.pt = abs(1 / qop) * math.sin(theta)
        self.cot_theta = 1 / math.tan(theta)
        # ACTS perigee d0 lies along z cross momentum. Positive charges turn
        # clockwise in the transverse plane. Coordinates and arc length are mm.
        normal = np.array([-math.sin(self.phi), math.cos(self.phi)])
        self.perigee_xy = d0 * normal
        self.curvature = 0.000299792458 * 2.0 * qop / math.sin(theta)
        self.center = self.perigee_xy - normal / self.curvature
        self.radius = abs(1 / self.curvature)

    def arc_to(self, xy):
        delta = self.curvature * (xy - self.perigee_xy)
        phi = math.atan2(math.sin(self.phi) - delta[0], math.cos(self.phi) + delta[1])
        return -math.remainder(phi - self.phi, 2 * math.pi) / self.curvature

    def position_tangent(self, arc):
        phi = self.phi - self.curvature * arc
        sin_phi, cos_phi = math.sin(phi), math.cos(phi)
        point = np.array(
            [
                self.perigee_xy[0] + (math.sin(self.phi) - sin_phi) / self.curvature,
                self.perigee_xy[1] + (cos_phi - math.cos(self.phi)) / self.curvature,
                self.z0 + arc * self.cot_theta,
            ]
        )
        return point, np.array([cos_phi, sin_phi, self.cot_theta])


def transverse_seeds(first, second):
    delta = second.center - first.center
    distance = np.linalg.norm(delta)
    if distance < 1e-9:
        return [first.perigee_xy, second.perigee_xy]
    unit = delta / distance
    along = (distance**2 + first.radius**2 - second.radius**2) / (2 * distance)
    height2 = first.radius**2 - along**2
    if height2 < 0:
        # Disjoint (possibly nested) circles. Start between their closest points.
        p1 = first.center + math.copysign(first.radius, along) * unit
        p2 = second.center + math.copysign(second.radius, along - distance) * unit
        return [(p1 + p2) / 2]
    base = first.center + along * unit
    offset = math.sqrt(height2) * np.array([-unit[1], unit[0]])
    return [base + offset, base - offset]


def closest_approach(first, second):
    """Minimize 3D separation on each track's turn nearest the perigee."""
    options = []
    for xy in transverse_seeds(first, second):
        arcs = np.array([first.arc_to(xy), second.arc_to(xy)])
        for _ in range(20):
            p1, t1 = first.position_tangent(arcs[0])
            p2, t2 = second.position_tangent(arcs[1])
            residual = p1 - p2
            step = np.linalg.lstsq(np.column_stack([t1, -t2]), -residual, rcond=None)[0]
            if np.linalg.norm(step) < 1e-6:
                break
            # Linearize in arc length. Backtrack if the linearization overshoots.
            for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
                trial = arcs + scale * step
                if (
                    max(
                        abs(first.curvature * trial[0]),
                        abs(second.curvature * trial[1]),
                    )
                    > math.pi
                ):
                    continue
                v1, _ = first.position_tangent(trial[0])
                v2, _ = second.position_tangent(trial[1])
                if np.linalg.norm(v1 - v2) < np.linalg.norm(residual):
                    arcs = trial
                    break
            else:
                break
        p1, t1 = first.position_tangent(arcs[0])
        p2, t2 = second.position_tangent(arcs[1])
        options.append(
            (np.linalg.norm(p1 - p2), (p1 + p2) / 2, first.pt * t1, second.pt * t2)
        )
    return min(options, key=lambda option: option[0])


def pair_mass_mev(plus, minus):
    energy = math.hypot(*plus, PION_MASS_GEV) + math.hypot(*minus, PION_MASS_GEV)
    total = plus + minus
    return 1000 * math.sqrt(energy**2 - total @ total)


def read_fits(filename):
    fit_fields = ("eLOC0_fit", "eLOC1_fit", "ePHI_fit", "eTHETA_fit", "eQOP_fit")
    labels = ["majorityParticleId_" + field for field in BARCODE]
    fields = [
        "trackClassification",
        "nMajorityHits",
        "nMeasurements",
        "hasFittedParams",
        "track_nr",
        *fit_fields,
        *labels,
    ]
    output = {}
    with filename.open("rb") as stream, uproot.open(stream) as root:
        for batch in root["tracksummary"].iterate(
            ["event_nr", *fields], step_size="40 MB"
        ):
            tracks = ak.zip({field: batch[field] for field in fields})
            selected = (
                (tracks.trackClassification == 1)
                & tracks.hasFittedParams
                & (2 * tracks.nMajorityHits > tracks.nMeasurements)
            )
            for event, rows in zip(
                batch.event_nr.to_list(), tracks[selected].to_list()
            ):
                matches = defaultdict(list)
                for row in rows:
                    if row["eQOP_fit"] != 0 and all(
                        math.isfinite(row[f]) for f in fit_fields
                    ):
                        matches[tuple(row[f] for f in labels)].append(row)
                output[event] = matches
    return output


def rows_for_pairs(candidates, directory, counts):
    fits = read_fits(directory / "tracksummary_ckf.root")
    fields = ["particle_type", "orig_part_idx", "process", *BARCODE]
    with (
        (directory / "particles_simulation.root").open("rb") as stream,
        uproot.open(stream) as root,
    ):
        for batch in root["particles"].iterate(
            ["event_id", *fields], step_size="20 MB"
        ):
            particles = ak.zip({field: batch[field] for field in fields})
            is_daughter = (
                (np.abs(particles.particle_type) == 211)
                & (particles.orig_part_idx > 0)
                & (particles.process == 1)
            )
            for event, daughters in zip(
                batch.event_id.to_list(), particles[is_daughter].to_list()
            ):
                lookup = {p["orig_part_idx"]: p for p in daughters}
                for row in candidates[event]:
                    plus, minus = [
                        lookup[int(row[key])] for key in ("plus_id", "minus_id")
                    ]
                    if (plus["particle_type"], minus["particle_type"]) != (211, -211):
                        raise ValueError(
                            f"Unexpected daughter PDG code in event {event}"
                        )
                    matched = [
                        fits[event].get(tuple(p[f] for f in BARCODE), [])
                        for p in (plus, minus)
                    ]
                    if any(len(m) == 1 for m in matched):
                        counts["at_least_one_fit"] += 1
                    if any(len(m) != 1 for m in matched):
                        continue
                    counts["two_unique_fits"] += 1
                    first, second = [m[0] for m in matched]
                    if first["track_nr"] == second["track_nr"]:
                        raise ValueError(f"Reused fitted track in event {event}")
                    if first["eQOP_fit"] <= 0 or second["eQOP_fit"] >= 0:
                        continue
                    counts["selected"] += 1
                    first_helix, second_helix = Helix(first), Helix(second)
                    dca, vertex, p1, p2 = closest_approach(first_helix, second_helix)
                    truth = np.array(
                        [float(row["g4_decay_" + axis + "_mm"]) for axis in "xyz"]
                    )
                    yield {
                        "event": event,
                        "kaon_id": row["kaon_id"],
                        "plus_track_id": first["track_nr"],
                        "minus_track_id": second["track_nr"],
                        "m_fit_at_dca_mev": pair_mass_mev(p1, p2),
                        "helix_dca_mm": dca,
                        "fitted_x_mm": vertex[0],
                        "fitted_y_mm": vertex[1],
                        "fitted_z_mm": vertex[2],
                        "fit_r_mm": np.linalg.norm(vertex[:2]),
                        "r_decay_mm": np.linalg.norm(truth[:2]),
                        "m_truth_mev": row["m_g4_mev"],
                        "fit_xyz_err_mm": np.linalg.norm(vertex - truth),
                    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates",
        required=True,
        type=Path,
        help="geant4_ks_candidates.py output CSV",
    )
    parser.add_argument(
        "--odd",
        required=True,
        type=Path,
        help="full_chain_odd.py ROOT output directory",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="CSV of truth-matched fitted pairs"
    )
    args = parser.parse_args()
    counts = Counter()
    pairs = defaultdict(list)
    with args.candidates.open() as stream:
        for row in csv.DictReader(stream):
            counts["generated"] += 1
            if int(row["pair"]):
                pairs[int(row["event"])].append(row)
                counts["truth_pairs"] += 1
                if min(int(row[k]) for k in ("plus_hits", "minus_hits")) > 0:
                    counts["both_have_hits"] += 1
    fields = (
        "event",
        "kaon_id",
        "plus_track_id",
        "minus_track_id",
        "m_fit_at_dca_mev",
        "helix_dca_mm",
        "fitted_x_mm",
        "fitted_y_mm",
        "fitted_z_mm",
        "fit_r_mm",
        "r_decay_mm",
        "m_truth_mev",
        "fit_xyz_err_mm",
    )
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_for_pairs(pairs, args.odd, counts))
    print(dict(counts))


if __name__ == "__main__":
    main()
