#!/usr/bin/env python3

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import addGbtsTraining
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

u = acts.UnitConstants

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "odd_gbts_training",
)
parser.add_argument("--events", "-n", help="Number of events", type=int, default=400)
parser.add_argument("--particles", help="Muons per event", type=int, default=50)
parser.add_argument(
    "--jobs",
    "-j",
    help="Number of worker threads (-1 uses all cores)",
    type=int,
    default=-1,
)
parser.add_argument(
    "--prob-threshold",
    help="Drop transitions below this probability (-1 keeps all)",
    type=float,
    default=-1.0,
)
args = parser.parse_args()

actsDir = pathlib.Path(__file__).parent.parent.parent.parent
configDir = actsDir / "Examples/Configs"
outputDir = args.output
outputDir.mkdir(parents=True, exist_ok=True)

oddDir = getOpenDataDetectorDirectory()
oddMaterialDeco = acts.IMaterialDecorator.fromFile(
    oddDir / "data/odd-material-maps.root"
)

detector = getOpenDataDetector(materialDecorator=oddMaterialDeco)
trackingGeometry = detector.trackingGeometry()
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args.events, numThreads=args.jobs, outputDir=str(outputDir)
)

addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
    EtaConfig(-3.0, 3.0, uniform=True),
    ParticleConfig(args.particles, acts.PdgParticle.eMuon, randomizeCharge=True),
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
    ),
    rnd=rnd,
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=configDir / "odd-digi-smearing-config.json",
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-3.0, 3.0),
        measurements=(3, None),
        removeNeutral=True,
    ),
)

addGbtsTraining(
    s,
    selectedParticles="particles_selected",
    geometryFile=configDir / "odd-gbts-layer-geometry.txt",
    outputConnectionTable=outputDir / "odd-gbts-connection-table.txt",
    probThreshold=args.prob_threshold,
    doSymmetrization=True,
    useOldFormatting=True,
    logLevel=acts.logging.INFO,
)

s.run()
