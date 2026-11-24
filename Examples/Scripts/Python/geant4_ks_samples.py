"""Generate Pythia8 K_S samples and transport a particle gun through Geant4."""

import argparse
from math import hypot
from pathlib import Path

KAON_DECAY = (
    "310:mayDecay = on",
    "310:onMode = off",
    "310:onIfMatch = 211 -211",
    "ParticleDecays:limitTau0 = off",
    "211:mayDecay = off",
)


def generate_gun(path, per_momentum, seed):
    import cppyy
    import pyhepmc

    # Use the installed Pythia headers and library through ROOT's cppyy.
    cppyy.load_library("libpythia8")
    cppyy.include("Pythia8/Pythia.h")
    pythia = cppyy.gbl.Pythia8.Pythia("", False)
    for setting in (
        "ProcessLevel:all = off",
        "Check:event = off",
        "Random:setSeed = on",
        f"Random:seed = {seed}",
        *KAON_DECAY,
        "Print:quiet = on",
    ):
        if not pythia.readString(setting):
            raise ValueError(f"Pythia rejected {setting}")
    if not pythia.init():
        raise RuntimeError("Pythia initialization failed")

    def copy_particle(source, status):
        particle = pyhepmc.GenParticle(
            (source.px(), source.py(), source.pz(), source.e()), source.id(), status
        )
        particle.generated_mass = source.m()
        return particle

    momenta = (0.5, 2.0, 5.0)  # GeV/c
    mass = pythia.particleData.m0(310)
    print(f"K_S mass {mass} GeV/c², cτ {pythia.particleData.tau0(310)} mm")
    with pyhepmc.io.WriterAscii(str(path)) as writer:
        writer.precision = 17
        for i in range(per_momentum * len(momenta)):
            p = momenta[i % len(momenta)]
            record = pythia.event
            record.reset()
            idx = record.append(310, 1, 0, 0, 0.0, 0.0, p, hypot(mass, p), mass)
            # Isolated particle gun with lifetime, as in Pythia8 example main234.
            record[idx].tau(record[idx].tau0() * pythia.rndm.exp())
            if not pythia.next():
                raise RuntimeError(f"Pythia decay failed at event {i}")
            kaon = record[idx]
            children = sorted(
                (record[j] for j in kaon.daughterList()),
                key=lambda child: child.id(),
                reverse=True,
            )
            if [child.id() for child in children] != [211, -211]:
                raise ValueError(f"Unexpected decay channel at event {i}")

            event = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
            event.event_number = i
            origin = pyhepmc.GenVertex()
            # Status 4 marks the incoming gun particle at the production vertex.
            origin.add_particle_in(copy_particle(kaon, 4))
            parent = copy_particle(kaon, 2)
            origin.add_particle_out(parent)
            event.add_vertex(origin)
            v = kaon.vDec()
            decay = pyhepmc.GenVertex((v.px(), v.py(), v.pz(), v.e()))
            decay.add_particle_in(parent)
            for child in children:
                decay.add_particle_out(copy_particle(child, 1))
            event.add_vertex(decay)
            writer.write_event(event)
            if writer.failed():
                raise OSError(f"Unable to write {path} at event {i}")


def generate_z(path, events, seed):
    import acts
    import acts.examples
    from acts.examples.hepmc3 import HepMC3Writer
    from acts.examples.simulation import addPythia8

    level = acts.logging.WARNING
    seq = acts.examples.Sequencer(
        events=events, numThreads=1, logLevel=level, trackFpes=False
    )
    addPythia8(
        seq,
        rnd=acts.examples.RandomNumbers(seed=seed),
        npileup=0,
        beam=[acts.PdgParticle.eElectron, acts.PdgParticle.ePositron],
        cmsEnergy=91.2,
        hardProcess=[
            "WeakSingleBoson:ffbar2gmZ = on",
            "23:onMode = off",
            "23:onIfAny = 1 2 3 4 5",
            *KAON_DECAY,
            "13:mayDecay = off",
        ],
    )
    seq.addWriter(
        HepMC3Writer(level=level, inputEvent="pythia8-event", outputPath=path)
    )
    seq.run()


def transport(path, output, workers):
    import acts
    import acts.examples
    from acts.examples.hepmc3 import HepMC3InputConverter, HepMC3Reader
    from acts.examples.simulation import addGeant4

    output.mkdir(parents=True, exist_ok=False)
    level = acts.logging.WARNING
    # All parents fly along +z in Galactic vacuum. Two 80 µm silicon planes
    # sit four metres off-axis, giving a large world without parent scattering.
    detector = acts.examples.TelescopeDetector(
        positions=[100, 200], offsets=[4000, 0], stereos=[0, 0], bounds=[25, 25]
    )
    seq = acts.examples.Sequencer(
        numThreads=max(1, workers), trackFpes=False, logLevel=level
    )
    seq.addReader(HepMC3Reader(level=level, inputPath=path, outputEvent="hepmc"))
    seq.addAlgorithm(
        HepMC3InputConverter(
            level=level,
            inputEvent="hepmc",
            outputParticles="roots",
            outputVertices="vertices",
            outputDecayGraph="decays",
            decayPdgCodes=[310],
        )
    )
    addGeant4(
        seq,
        detector,
        detector.trackingGeometry(),
        acts.ConstantBField(acts.Vector3(0, 0, 0)),
        acts.examples.RandomNumbers(seed=42),
        inputParticles="roots",
        inputDecayGraph="decays",
        numGeant4Workers=workers,
        outputDirRoot=output,
        logLevel=level,
    )
    seq.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    gen = sub.add_parser("generate-gun", help="Generate a K_S particle gun in HepMC3")
    gen.add_argument("output", type=Path)
    gen.add_argument("--per-momentum", type=int, default=4000)
    gen.add_argument("--seed", type=int, default=240510)
    z = sub.add_parser("generate-z", help="Generate hadronic Z events in HepMC3")
    z.add_argument("output", type=Path)
    z.add_argument("--events", type=int, required=True)
    z.add_argument("--seed", type=int, default=42)
    sim = sub.add_parser(
        "transport", help="Run ACTS Geant4 and write particle truth to ROOT"
    )
    sim.add_argument("input", type=Path)
    sim.add_argument("output", type=Path)
    sim.add_argument(
        "--workers", type=int, default=4, help="Zero selects serial Geant4"
    )
    args = parser.parse_args()
    if args.command == "transport":
        if args.workers < 0:
            parser.error("--workers must be nonnegative")
        transport(args.input, args.output, args.workers)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.exists():
            parser.error(f"Output file already exists: {args.output}")
        if args.command == "generate-gun":
            if args.per_momentum < 1:
                parser.error("--per-momentum must be positive")
            generate_gun(args.output, args.per_momentum, args.seed)
        else:
            if args.events < 1:
                parser.error("--events must be positive")
            generate_z(args.output, args.events, args.seed)
