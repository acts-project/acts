#!/usr/bin/env python3

from pathlib import Path
import json

import acts
import acts.examples

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    addFatras,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
)

from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
    addTrackWriters,
)

u = acts.UnitConstants


if "__main__" == __name__:
    # GenericDetector
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    # Create minimal digitization config with only volume id 0 (broadcasts to all surfaces)
    digiConfig = {
        "acts-geometry-hierarchy-map": {
            "format-version": 0,
            "value-identifier": "digitization-configuration",
        },
        "entries": [
            {
                "volume": 0,
                "value": {
                    "smearing": [
                        {"index": 0, "mean": 0.0, "stddev": 0.01, "type": "Gauss"},
                        {"index": 1, "mean": 0.0, "stddev": 0.01, "type": "Gauss"},
                    ]
                },
            }
        ],
    }

    # Save config to file
    digiConfigFile = Path.cwd() / "digi-config-minimal.json"
    with open(digiConfigFile, "w") as f:
        json.dump(digiConfig, f, indent=4)
    logger = acts.getDefaultLogger("Truth tracking example", acts.logging.INFO)
    logger.info(f"Saved minimal digitization config to {digiConfigFile}")

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    s = acts.examples.Sequencer(events=10, numThreads=-1, logLevel=acts.logging.INFO)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
        PhiConfig(0.0, 360.0 * u.degree),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0, 0, 0, 0),
        ),
        multiplicity=1,
        rnd=rnd,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
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
            pt=(0.9 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
    )

    addTrackWriters(
        s,
        name="tracks",
        outputDirCsv="csv",
    )

    s.run()
