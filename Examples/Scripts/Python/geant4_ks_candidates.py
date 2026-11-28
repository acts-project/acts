"""Match prescribed HepMC3 K_S -> pi+ pi- decays to Geant4 particle truth.

Write one CSV row per generated charged-pion decay, including rejected roots
and parents whose prescribed decay did not occur. The input files must contain
the same event IDs. Momenta are evaluated at daughter creation. The proper
flight lengths use the elapsed lab time times c divided by the initial kaon
Lorentz factor. The helicity axes follow each pair's lab momentum.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import pyhepmc
import uproot

DECAY_PROCESS = 1  # ActsFatras::GenerationProcess
NUCLEAR_PROCESS = 4

BARCODE = (
    "vertex_primary",
    "vertex_secondary",
    "particle",
    "generation",
    "sub_particle",
)
KINEMATICS = (
    "px",
    "py",
    "pz",
    "m",
    "vx",
    "vy",
    "vz",
    "vt",
    "particle_type",
    "orig_part_idx",
    "process",
    "outcome",
    "number_of_hits",
)
PAIR_FIELDS = (
    "g4_L_mm",
    "g4_ctau_mm",
    "m_g4_mev",
    "pstar_g4_mev",
    "g4_alpha",
    "g4_qt_mev",
    "g4_costh",
    "g4_piplus_energy_gev",
    "plus_dp_mev",
    "minus_dp_mev",
    "parent_dp_mev",
    "vertex_delta_mm",
    "daughter_vertex_split_mm",
    "g4_decay_x_mm",
    "g4_decay_y_mm",
    "g4_decay_z_mm",
    "g4_decay_ct_mm",
    "plus_hits",
    "minus_hits",
)


def pair_observables(plus, minus):
    """Dipion mass, rest momentum, Armenteros variables and pi+ helicity angle."""
    pair = plus + minus
    p = np.linalg.norm(pair[:3])
    mass = np.sqrt((pair[3] - p) * (pair[3] + p))
    rest = (
        plus[:3]
        + ((plus[:3] @ pair[:3] / (pair[3] + mass) - plus[3]) / mass) * pair[:3]
    )
    pstar = np.linalg.norm(rest)
    if p < 1e-12:
        return mass * 1000, pstar * 1000, None, None, None
    axis = pair[:3] / p
    longitudinal = plus[:3] @ axis
    alpha = (2 * longitudinal - p) / p
    qt = np.linalg.norm(plus[:3] - longitudinal * axis)
    return mass * 1000, pstar * 1000, alpha, qt * 1000, rest @ axis / pstar


def generator_decays(path):
    events = {}
    with pyhepmc.open(str(path)) as stream:
        for event in stream:
            event.set_units(pyhepmc.Units.GEV, pyhepmc.Units.MM)
            if event.event_number in events:
                raise ValueError(f"Duplicate HepMC event {event.event_number}")
            decays = []
            for kaon in event.particles:
                if kaon.pid != 310 or kaon.status != 2 or kaon.end_vertex is None:
                    continue
                children = kaon.end_vertex.particles_out
                if sorted(p.pid for p in children) != [-211, 211]:
                    continue
                if len(kaon.end_vertex.particles_in) != 1 or any(
                    p.status != 1 for p in children
                ):
                    raise ValueError(
                        f"Invalid charged-pion decay at HepMC particle {kaon.id}"
                    )
                plus, minus = sorted(children, key=lambda p: -p.pid)
                if kaon.production_vertex is None:
                    raise ValueError(
                        f"Missing production vertex at HepMC particle {kaon.id}"
                    )
                decays.append(
                    (
                        (kaon.id, plus.id, minus.id),
                        [np.array(list(p.momentum)) for p in (kaon, plus, minus)],
                        np.array(list(kaon.production_vertex.position)),
                        np.array(list(kaon.end_vertex.position)),
                        kaon.generated_mass,
                    )
                )
            events[event.event_number] = decays
    return events


def vector(particles, fields, index):
    return np.array([particles[name][index] for name in fields], dtype=float)


def sim_momentum(particles, i):
    p = vector(particles, ("px", "py", "pz"), i)
    return np.append(p, np.hypot(np.linalg.norm(p), particles["m"][i]))


def candidate(event, decay, sim):
    (kid, plus_id, minus_id), (kaon, plus, minus), start, end, mass = decay
    beta_gamma = np.linalg.norm(kaon[:3]) / mass
    inv_gamma = mass / kaon[3]
    m, pstar, alpha, qt, costh = pair_observables(plus, minus)
    record = dict(
        event=event,
        kaon_id=kid,
        plus_id=plus_id,
        minus_id=minus_id,
        transported=0,
        pair=0,
        status="not_transported",
        outcome=None,
        direct_nuclear_children=0,
        gen_p_gev=np.linalg.norm(kaon[:3]),
        beta_gamma=beta_gamma,
        gen_L_mm=np.linalg.norm(end[:3] - start[:3]),
        gen_ctau_mm=(end[3] - start[3]) * inv_gamma,
        gen_mass_kaon_mev=mass * 1000,
        m_gen_mev=m,
        pstar_gen_mev=pstar,
        gen_alpha=alpha,
        gen_qt_mev=qt,
        gen_costh=costh,
        gen_piplus_energy_gev=plus[3],
        gen_decay_x_mm=end[0],
        gen_decay_y_mm=end[1],
        gen_decay_z_mm=end[2],
        gen_decay_ct_mm=end[3],
    )
    record.update(dict.fromkeys(PAIR_FIELDS))
    pdg, source, process = (
        sim[k] for k in ("particle_type", "orig_part_idx", "process")
    )
    roots = np.flatnonzero((source == kid) & (pdg == 310))
    plus_rows = np.flatnonzero((source == plus_id) & (pdg == 211))
    minus_rows = np.flatnonzero((source == minus_id) & (pdg == -211))
    if not roots.size:
        if plus_rows.size or minus_rows.size:
            raise ValueError(
                f"Missing transported parent in event {event} source {kid}"
            )
        return record
    if len(roots) != 1:
        raise ValueError(f"Multiple kaons matched event {event} source {kid}")
    root = roots[0]
    direct = np.logical_and.reduce(
        [sim["parent_" + col] == sim[col][root] for col in BARCODE]
    )
    outcome = int(sim["outcome"][root])
    nuclear = np.count_nonzero(direct & (process == NUCLEAR_PROCESS))
    record.update(transported=1, outcome=outcome, direct_nuclear_children=nuclear)
    if not (plus_rows.size or minus_rows.size):
        record["status"] = (
            "nuclear_products"
            if nuclear
            else {2: "volume_exit", 3: "time_cut"}.get(
                outcome, "no_prescribed_daughters"
            )
        )
        return record
    if plus_rows.size != 1 or minus_rows.size != 1:
        raise ValueError(
            f"Incomplete or duplicate pion pair in event {event} source {kid}"
        )
    i, j = plus_rows[0], minus_rows[0]
    for k in (i, j):
        if not direct[k] or process[k] != DECAY_PROCESS:
            raise ValueError(
                f"Incorrect pion parent or process in event {event} source {kid}"
            )
    a, b = sim_momentum(sim, i), sim_momentum(sim, j)
    coords = ("vx", "vy", "vz", "vt")
    va, vb, v0 = (vector(sim, coords, index) for index in (i, j, root))
    m, pstar, alpha, qt, costh = pair_observables(a, b)
    record.update(
        pair=1,
        status="decayed",
        g4_L_mm=np.linalg.norm(va[:3] - v0[:3]),
        g4_ctau_mm=(va[3] - v0[3]) * inv_gamma,
        m_g4_mev=m,
        pstar_g4_mev=pstar,
        g4_alpha=alpha,
        g4_qt_mev=qt,
        g4_costh=costh,
        g4_piplus_energy_gev=a[3],
        plus_dp_mev=np.linalg.norm(a[:3] - plus[:3]) * 1000,
        minus_dp_mev=np.linalg.norm(b[:3] - minus[:3]) * 1000,
        parent_dp_mev=np.linalg.norm(a[:3] + b[:3] - kaon[:3]) * 1000,
        vertex_delta_mm=np.linalg.norm(va[:3] - end[:3]),
        daughter_vertex_split_mm=np.linalg.norm(va[:3] - vb[:3]),
        g4_decay_x_mm=va[0],
        g4_decay_y_mm=va[1],
        g4_decay_z_mm=va[2],
        g4_decay_ct_mm=va[3],
        plus_hits=int(sim["number_of_hits"][i]),
        minus_hits=int(sim["number_of_hits"][j]),
    )
    return record


def match(hepmc, particles, output):
    events = generator_decays(hepmc)
    branches = (*KINEMATICS, *BARCODE, *("parent_" + col for col in BARCODE))
    rows = []
    with particles.open("rb") as stream, uproot.open(stream) as file:
        for batch in file["particles"].iterate(
            ["event_id", *branches], step_size="32 MB"
        ):
            for entry in batch:
                event = int(entry.event_id)
                if event not in events:
                    raise ValueError(f"Unmatched or duplicate ROOT event {event}")
                decays = events.pop(event)
                if decays:
                    sim = {col: np.asarray(entry[col]) for col in branches}
                    rows.extend(candidate(event, d, sim) for d in decays)
    if events:
        raise ValueError(f"{len(events)} HepMC events missing from ROOT")
    if not rows:
        raise ValueError("No generated K_S -> pi+ pi- decays found")
    rows.sort(key=lambda row: (row["event"], row["kaon_id"]))
    with output.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hepmc", type=Path, required=True)
    parser.add_argument("--particles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    match(args.hepmc, args.particles, args.output)
