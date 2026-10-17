#!/usr/bin/env python3

import argparse
import inspect
from pathlib import Path

import acts
import acts.examples

u = acts.UnitConstants


def runGbtsTrainingOdd(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    layerMapFile: Path,
    outputDir: Path,
    events: int = 2000,
    numParticles: int = 50,
    probThreshold: float = 1e-4,
    useLayerTilt: bool = False,
    numThreads: int = -1,
    decorators=[],
):
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

    outputDir.mkdir(parents=True, exist_ok=True)

    s = acts.examples.Sequencer(
        events=events, numThreads=numThreads, outputDir=str(outputDir)
    )

    for d in decorators:
        s.addContextDecorator(d)

    # not the seed of the physmon samples, so the table is not trained on them
    rnd = acts.examples.RandomNumbers(seed=1234)

    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        ParticleConfig(numParticles, acts.PdgParticle.eMuon, randomizeCharge=True),
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
        digiConfigFile=digiConfigFile,
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
        trackingGeometry,
        selectedParticles="particles_selected",
        geometryFile=layerMapFile,
        outputConnectionTable=outputDir / "odd-gbts-connection-table.json",
        probThreshold=probThreshold,
        doSymmetrization=True,
        useLayerTilt=useLayerTilt,
        logLevel=acts.logging.INFO,
    )

    return s


if "__main__" == __name__:
    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    defaults = inspect.signature(runGbtsTrainingOdd).parameters

    parser = argparse.ArgumentParser(
        description="Layer connection training tool for the GBTS on ODD"
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output directory",
        type=Path,
        default=Path.cwd() / "odd_gbts_training",
    )
    parser.add_argument(
        "--events",
        "-n",
        help="Number of events",
        type=int,
        default=defaults["events"].default,
    )
    parser.add_argument(
        "--particles",
        help="Muons per event",
        type=int,
        default=defaults["numParticles"].default,
    )
    parser.add_argument(
        "--jobs",
        "-j",
        help="Number of worker threads (-1 uses all cores)",
        type=int,
        default=defaults["numThreads"].default,
    )
    parser.add_argument(
        "--prob-threshold",
        help="Drop transitions below this probability (-1 keeps all)",
        type=float,
        default=defaults["probThreshold"].default,
    )
    parser.add_argument(
        "--use-layer-tilt",
        help="Cut the eta bins along the module tilt, as the seeding does with "
        "useLayerTilt",
        action="store_true",
    )
    args = parser.parse_args()

    configDir = (
        Path(__file__).resolve().parent.parent.parent.parent / "Examples/Configs"
    )

    oddDir = getOpenDataDetectorDirectory()
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(
        oddDir / "data/odd-material-maps.root"
    )
    detector = getOpenDataDetector(materialDecorator=oddMaterialDeco)

    runGbtsTrainingOdd(
        trackingGeometry=detector.trackingGeometry(),
        field=acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T)),
        digiConfigFile=configDir / "odd-digi-smearing-config.json",
        layerMapFile=configDir / "odd-gbts-layer-map.json",
        outputDir=args.output,
        events=args.events,
        numParticles=args.particles,
        probThreshold=args.prob_threshold,
        useLayerTilt=args.use_layer_tilt,
        numThreads=args.jobs,
        decorators=detector.contextDecorators(),
    ).run()
