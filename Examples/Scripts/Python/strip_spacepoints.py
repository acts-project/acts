#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples

u = acts.UnitConstants


def createStripSpacepoints(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    geoSelection: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputHitsPath: Optional[Path] = None,
    decorators=[],
    reverseFilteringMomThreshold=0 * u.GeV,
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    logger = acts.logging.getLogger("Truth tracking example")

    if inputParticlePath is None:
        # Note: We restrict the eta range to [-2,2] to get tracks with long-strip hits
        addParticleGun(
            s,
            ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            EtaConfig(-2.0, 2.0, uniform=True),
            MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=1,
            rnd=rnd,
        )
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        s.addReader(
            acts.examples.RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )
        s.addWhiteboardAlias("particles", "particles_generated")

    if inputHitsPath is None:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
            enableInteractions=True,
        )
    else:
        logger.info("Reading hits from %s", inputHitsPath.resolve())
        assert inputHitsPath.exists()
        s.addReader(
            acts.examples.RootSimHitReader(
                level=acts.logging.INFO,
                filePath=str(inputHitsPath.resolve()),
                outputSimHits="simhits",
            )
        )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    # Create strip selection of ODD:
    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=acts.logging.INFO,
            trackingGeometry=trackingGeometry,
            inputMeasurements="measurements",
            outputSpacePoints="spacepoints",
            stripGeometrySelection=acts.examples.readJsonGeometryList(
                str(geoSelection)
            ),
        )
    )

    s.addWriter(
        acts.examples.RootSpacepointWriter(
            level=acts.logging.INFO,
            inputSpacepoints="spacepoints",
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDir / "strip_spacepoints.root"),
        )
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    digiConfigFile = (
        srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
    )

    geoSelection = srcdir / "Examples/Configs/odd-strip-spacepoint-selection.json"
    print(geoSelection.resolve())
    assert geoSelection.exists()

    ## GenericDetector
    # detector = acts.examples.GenericDetector()
    # trackingGeometry = detector.trackingGeometry()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    createStripSpacepoints(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        geoSelection=geoSelection,
        outputDir=Path.cwd(),
    ).run()
