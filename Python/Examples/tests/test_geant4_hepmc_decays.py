"""Transport a displaced pion and its pre-decayed muon through three sensors."""

import math
from pathlib import Path
import subprocess
import sys

import pytest
from helpers import geant4Enabled, hepmc3Enabled, rootEnabled


def _boost(momentum, parent, mass):
    """Boost a rest-frame (px, py, pz, E) into the parent's laboratory frame."""
    dot = sum(a * b for a, b in zip(momentum[:3], parent[:3]))
    scale = (momentum[3] + dot / (parent[3] + mass)) / mass
    return [momentum[i] + scale * parent[i] for i in range(3)] + [
        (parent[3] * momentum[3] + dot) / mass
    ]


def _event(event):
    # Momenta in MeV, space and ct in mm. Each parent lives for c*tau = 10 mm.
    pion_mass, mu_mass, e_mass = 139.57039, 105.6583755, 0.51099895
    pion = [0, 0, 500, math.hypot(500, pion_mass)]
    q = (pion_mass**2 - mu_mass**2) / (2 * pion_mass)
    muon = _boost([q, 0, 0, math.hypot(q, mu_mass)], pion, pion_mass)
    nu = _boost([-q, 0, 0, q], pion, pion_mass)
    # Event 2 decays inside the first silicon plane; event 3 leaves the world.
    pion_ctau = {2: 10 * pion_mass / pion[2], 3: 1e5}.get(event, 10)
    pion_end = [pion_ctau * p / pion_mass for p in pion]
    muon_end = [x + 10 * p / mu_mass for x, p in zip(pion_end, muon)]

    # Three-body decay with the rest-frame momenta separated by 120 degrees.
    q = (2 * mu_mass - math.sqrt(mu_mass**2 + 3 * e_mass**2)) / 3
    daughters = [
        _boost([q, 0, 0, math.hypot(q, e_mass)], muon, mu_mass),
        _boost([-q / 2, q * math.sqrt(3) / 2, 0, q], muon, mu_mass),
        _boost([-q / 2, -q * math.sqrt(3) / 2, 0, q], muon, mu_mass),
    ]
    # Alternate the units to exercise the file reader's event merge.
    p_scale, x_scale = (1, 0.1) if event % 2 else (0.001, 1)

    def vector(v, scale):
        return " ".join(f"{x * scale:.17g}" for x in v)

    def particle(idx, vertex, pdg, vec, mass, status):
        return f"P {idx} {vertex} {pdg} {vector(vec, p_scale)} {mass * p_scale:.17g} {status}"

    lines = [
        f"E {event} 3 7",
        "U MEV CM" if event % 2 else "U GEV MM",
        particle(1, 0, 22, [0, 0, 600, 600], 0, 4),
        "V -1 0 [1] @ 0 0 0 0",
        particle(2, -1, 211, pion, pion_mass, 2),
        f"V -2 0 [2] @ {vector(pion_end, x_scale)}",
        particle(3, -2, -13, muon, mu_mass, 2),
        particle(4, -2, 14, nu, 0, 1),
        f"V -3 0 [3] @ {vector(muon_end, x_scale)}",
        particle(5, -3, -11, daughters[0], e_mass, 1),
        particle(6, -3, 12, daughters[1], 0, 1),
        particle(7, -3, -14, daughters[2], 0, 1),
    ]
    if event == 4:
        # A neutral parent at rest spends 5 mm/c at its production point.
        kaon_mass = 497.611
        # Pythia rounds the charged pion mass to 139.57 MeV.
        pion_mass = 139.57
        q = math.sqrt((kaon_mass / 2) ** 2 - pion_mass**2)
        lines = [
            f"E {event} 2 4",
            "U GEV MM",
            particle(1, 0, 22, [0, 0, 600, 600], 0, 4),
            "V -1 0 [1] @ 0 0 0 0",
            particle(2, -1, 310, [0, 0, 0, kaon_mass], kaon_mass, 2),
            "V -2 0 [2] @ 0 0 0 5",
            particle(3, -2, 211, [0, 0, q, kaon_mass / 2], pion_mass, 1),
            particle(4, -2, -211, [0, 0, -q, kaon_mass / 2], pion_mass, 1),
        ]
    elif event == 5:
        # Stopped negative pion inside silicon. Nuclear capture competes with
        # a deliberately late muon decay; only the winning channel may appear.
        q = (pion_mass**2 - mu_mass**2) / (2 * pion_mass)
        lines = [
            f"E {event} 2 4",
            "U MEV CM",
            particle(1, 0, 22, [0, 0, 600, 600], 0, 4),
            "V -1 0 [1] @ 0 0 1 0",
            particle(2, -1, -211, [0, 0, 0, pion_mass], pion_mass, 2),
            "V -2 0 [2] @ 0 0 1 1000000",
            particle(3, -2, 13, [q, 0, 0, math.hypot(q, mu_mass)], mu_mass, 1),
            particle(4, -2, -14, [-q, 0, 0, q], 0, 1),
        ]
    return lines, [pion_end, muon_end], [pion, muon, nu, *daughters]


def _simulate(path, workers, events=6, field_tesla=0):
    import acts
    import acts.examples
    from acts.examples.hepmc3 import HepMC3InputConverter, HepMC3Reader
    from acts.examples.simulation import addGeant4

    level = acts.logging.WARNING
    detector = acts.examples.TelescopeDetector(
        positions=[10, 50, 90], stereos=[0] * 3, bounds=[100, 100]
    )
    s = acts.examples.Sequencer(
        events=events, numThreads=4 if workers else 1, trackFpes=False, logLevel=level
    )
    s.addReader(
        HepMC3Reader(level=level, inputPath=path / "decays.hepmc3", outputEvent="gen")
    )
    s.addAlgorithm(
        HepMC3InputConverter(
            level=level,
            inputEvent="gen",
            outputParticles="roots",
            outputVertices="vertices",
            outputDecayGraph="decays",
            decayPdgCodes=[211, -211, -13, 310],
        )
    )
    output = path / str(workers)
    addGeant4(
        s,
        detector,
        detector.trackingGeometry(),
        acts.ConstantBField(acts.Vector3(0, field_tesla * acts.UnitConstants.T, 0)),
        rnd=acts.examples.RandomNumbers(seed=42),
        inputParticles="roots",
        inputDecayGraph="decays",
        numGeant4Workers=workers,
        outputDirRoot=output,
    )
    s.run()

    # Restore the direct parent through the public ROOT reader and writer.
    from acts.examples.root import RootParticleReader, RootParticleWriter

    readback = acts.examples.Sequencer(events=events, numThreads=1, logLevel=level)
    readback.addReader(
        RootParticleReader(
            level=level,
            filePath=str(output / "particles_simulation.root"),
            outputParticles="reread",
        )
    )
    readback.addWriter(
        RootParticleWriter(
            level=level,
            filePath=str(output / "readback.root"),
            inputParticles="reread",
            writeParentBarcodes=True,
        )
    )
    readback.run()


def _transport(path, workers, events, field_tesla=0):
    records = ["HepMC::Version 3.03.01", "HepMC::Asciiv3-START_EVENT_LISTING"]
    for event in range(events):
        lines, _, _ = _event(event)
        records.extend(lines)
    records.append("HepMC::Asciiv3-END_EVENT_LISTING")
    (path / "decays.hepmc3").write_text("\n".join(records) + "\n")
    # Geant4 owns one run manager per process.
    run = subprocess.run(
        [
            sys.executable,
            __file__,
            str(path),
            str(workers),
            str(events),
            str(field_tesla),
        ],
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    return path / str(workers)


def _load(path, tree):
    import awkward as ak
    import uproot

    with path.open("rb") as stream, uproot.open(stream) as f:
        data = f[tree].arrays()
        # ROOT writers may finish events in a different order on each run.
        return data[ak.argsort(data.event_id, stable=True)]


pytestmark = pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled),
    reason="requires Geant4, HepMC3 and ROOT",
)


def test_prescribed_decay_chain(tmp_path):
    import awkward as ak
    import numpy as np
    from acts.examples.geant4 import nativeWorkersAvailable

    fields = [
        "vertex_primary",
        "vertex_secondary",
        "particle",
        "generation",
        "sub_particle",
    ]
    for workers in ([0, 2] if nativeWorkersAvailable else [0]):
        output = _transport(tmp_path, workers, events=6)
        particles = _load(output / "particles_simulation.root", "particles")
        hits = _load(output / "hits.root", "hits")
        restored = _load(output / "readback.root", "particles")
        for f in fields:
            assert ak.to_list(restored["parent_" + f]) == ak.to_list(
                particles["parent_" + f]
            )
        assert ak.to_list(particles.event_id) == list(range(6))

        for event, row in enumerate(particles):
            sources = list(row.orig_part_idx)
            if event == 3:
                assert row.particle_type[sources.index(2)] == 211
                assert all(i not in sources for i in (3, 4, 5, 6, 7))
                continue
            if event == 5:
                parent = sources.index(2)
                assert row.particle_type[parent] == -211
                assert 3 not in sources and 4 not in sources
                assert any(
                    row.process[i] == 4
                    and all(row["parent_" + f][i] == row[f][parent] for f in fields)
                    for i in range(len(sources))
                )
                continue
            if event == 4:
                k = sources.index(2)
                assert row.particle_type[k] == 310 and row.p[k] == 0
                assert all(
                    math.isnan(row[name][k])
                    for name in ("q_over_p", "eta", "theta", "phi")
                )
                expected_p = math.sqrt((0.497611 / 2) ** 2 - 0.13957**2)
                for source, pdg in ((3, 211), (4, -211)):
                    child = sources.index(source)
                    assert row.particle_type[child] == pdg
                    assert row.m[child] == pytest.approx(0.13957, abs=2e-8, rel=0)
                    assert row.p[child] == pytest.approx(expected_p, abs=3e-8, rel=0)
                    assert row.vt[child] == pytest.approx(5.0, abs=1e-5)
                    assert [row[f][child] for f in ("vx", "vy", "vz")] == pytest.approx(
                        [0] * 3
                    )
                    assert [row["parent_" + f][child] for f in fields] == [
                        row[f][k] for f in fields
                    ]
                continue
            _, vertices, momenta = _event(event)
            pdgs = [211, -13, 14, -11, 12, -14]
            indices = []
            for source, pdg in enumerate(pdgs, 2):
                matches = [
                    i
                    for i, value in enumerate(sources)
                    if value == source and row.particle_type[i] == pdg
                ]
                assert len(matches) == 1
                indices.append(matches[0])
            pion, muon, nu, positron, nue, antinu = indices
            for child, parent in (
                (muon, pion),
                (nu, pion),
                (positron, muon),
                (nue, muon),
                (antinu, muon),
            ):
                assert row.process[child] == 1  # G4Decay
                assert [row["parent_" + f][child] for f in fields] == [
                    row[f][parent] for f in fields
                ]
                vertex = vertices[0 if parent == pion else 1]
                assert [
                    row[f][child] for f in ("vx", "vy", "vz", "vt")
                ] == pytest.approx(vertex, abs=0.8)
            # The thin silicon planes allow a small amount of multiple scattering and energy loss.
            for i, expected in zip(indices, momenta):
                assert [row[f][i] * 1000 for f in ("px", "py", "pz")] == pytest.approx(
                    expected[:3], abs=6.0
                )
            event_hits = hits[hits.event_id == event]
            # A parent and daughter both leave a hit when the decay splits a sensor crossing.
            planes = ([10], [10], [50, 90]) if event == 2 else ([10], [50], [90])
            for i, expected_z in zip((pion, muon, positron), planes):
                matched = event_hits
                for f in fields:
                    matched = matched[matched["barcode_" + f] == row[f][i]]
                assert list(matched.tz) == pytest.approx(expected_z, abs=0.1)
                assert row.number_of_hits[i] == len(expected_z)
        if workers == 0:
            serial_particles, serial_hits = particles, hits
        else:
            for actual, expected in (
                (particles, serial_particles),
                (hits, serial_hits),
            ):
                for field in actual.fields:
                    # Compare NaNs (e.g. q/p for a neutral particle at rest) as equal.
                    np.testing.assert_array_equal(
                        ak.flatten(actual[field], axis=None),
                        ak.flatten(expected[field], axis=None),
                    )


def test_prescribed_decay_in_transverse_field(tmp_path):
    from acts.examples.geant4 import nativeWorkersAvailable

    # Both charged parents bend before decaying. Check proper time and the
    # boost of the daughters using the parent's actual production momentum.
    output = _transport(
        tmp_path, 2 if nativeWorkersAvailable else 0, events=2, field_tesla=5
    )
    particles = _load(output / "particles_simulation.root", "particles")
    assert list(particles.event_id) == [0, 1]

    ctau = 10.0  # mm
    k = 0.299792458 * 5.0  # q B in MeV/mm for charge +e and B along +y
    for row in particles:
        sources = list(row.orig_part_idx)

        def momentum(source):
            i = sources.index(source)
            return [1000 * row[f][i] for f in ("px", "py", "pz")]  # MeV

        def position(source):
            i = sources.index(source)
            return [row[f][i] for f in ("vx", "vy", "vz", "vt")]  # mm

        for parent, daughters, mass in (
            (2, (3, 4), 139.57039),
            (3, (5, 6, 7), 105.6583755),
        ):
            x, y, z, ct = position(parent)
            px, py, pz = momentum(parent)
            angle = k * ctau / mass
            sin, cos = math.sin(angle), math.cos(angle)
            # Integrate the helix in proper time and rotate the parent momentum.
            expected_vertex = [
                x + (px * sin + pz * (cos - 1)) / k,
                y + py * ctau / mass,
                z + (pz * sin + px * (1 - cos)) / k,
                ct + math.hypot(px, py, pz, mass) * ctau / mass,
            ]
            expected_momentum = [px * cos - pz * sin, py, pz * cos + px * sin]
            total_momentum = [
                sum(component) for component in zip(*(momentum(d) for d in daughters))
            ]
            # Allow for energy loss and multiple scattering in the thin sensors.
            for daughter in daughters:
                assert position(daughter) == pytest.approx(expected_vertex, abs=0.8)
            assert total_momentum == pytest.approx(expected_momentum, abs=6.0)


if __name__ == "__main__":
    _simulate(Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]))
