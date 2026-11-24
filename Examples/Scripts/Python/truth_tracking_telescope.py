#!/usr/bin/env python3

from pathlib import Path
import json
import tempfile

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:
    gen3 = False  # Set to True for Gen3 geometry, False for Gen1 geometry
    detector = acts.examples.TelescopeDetector(
        bounds=[200, 200],
        positions=[30, 60, 90, 120, 150, 180, 210, 240, 270],
        stereos=[0] * 9,
        rotDirection=2,
        gen3=gen3,
    )
    trackingGeometry = detector.trackingGeometry()

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # referenceSurface = acts.Surface.createPerigee(acts.Vector3(0, 0, 0))
    # Create reference surface of type Plane
    referenceSurface = detector.getReferenceSurface(acts.Vector3(0, 0, 0), 100, 100)

    if gen3:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "telescopeGen3-digi-smearing-config.json"

            nVolumes = 9
            config = {
                "acts-geometry-hierarchy-map": {
                    "format-version": 0,
                    "value-identifier": "digitization-configuration",
                },
                "entries": [
                    {
                        "volume": i,
                        "value": {
                            "smearing": [
                                {
                                    "index": 0,
                                    "mean": 0.0,
                                    "stddev": 0.025,
                                    "type": "Gauss",
                                },
                                {
                                    "index": 1,
                                    "mean": 0.0,
                                    "stddev": 0.1,
                                    "type": "Gauss",
                                },
                            ]
                        },
                    }
                    for i in range(1, nVolumes + 1)
                ],
            }

            # write JSON
            with open(config_path, "w", encoding="utf-8") as f:
                json.dump(config, f, indent=4)

            runTruthTrackingKalman(
                trackingGeometry,
                field,
                digiConfigFile=config_path,
                outputDir=Path.cwd(),
                referenceSurface=referenceSurface,
            ).run()
    else:
        runTruthTrackingKalman(
            trackingGeometry,
            field,
            digiConfigFile=srcdir
            / "Examples/Configs/telescope-digi-smearing-config.json",
            outputDir=Path.cwd(),
            referenceSurface=referenceSurface,
        ).run()
