#!/usr/bin/env python3
import argparse
import pathlib, acts, acts.examples
import acts.examples.dd4hep
from common import getOpenDataDetector, getOpenDataDetectorDirectory

parser = argparse.ArgumentParser(description="OpenDataDetector full chain example")
parser.add_argument(
    "-n", "--events", type=int, default=100, help="Number of events to run"
)
parser.add_argument(
    "-s", "--skip", type=int, default=0, help="Number of events to skip"
)
parser.add_argument(
    "-j",
    "--jobs",
    type=int,
    default=-1,
    help="Number of threads to use. Default: -1 i.e. number of cores",
)
parser.add_argument(
    "-o",
    "--output",
    type=pathlib.Path,
    default=pathlib.Path.cwd(),
    help="Output directories. Default: $PWD",
)
args = parser.parse_args()


u = acts.UnitConstants
outputDir = args.output

oddDir = getOpenDataDetectorDirectory()

oddMaterialMap = oddDir / "data/odd-material-maps.root"
oddDigiConfig = oddDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = oddDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(mdecorator=oddMaterialDeco)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

from particle_gun import addParticleGun, MomentumConfig, EtaConfig, ParticleConfig
from fatras import addFatras
from digitization import addDigitization
from seeding import addSeeding, SeedingAlgorithm, TruthSeedRanges
from ckf_tracks import addCKFTracks

s = acts.examples.Sequencer(events=args.events, numThreads=args.jobs, skip=args.skip)
s = addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, True),
    EtaConfig(-3.0, 3.0, True),
    ParticleConfig(1, acts.PdgParticle.eMuon, True),
    rnd=rnd,
)
s = addFatras(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-2.7, 2.7), nHits=(9, None)),
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
    initialVarInflation=[100, 100, 100, 100, 100, 100],
)
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(400.0 * u.MeV, None), nHits=(6, None)),
    outputDirRoot=outputDir,
)

s.run()
