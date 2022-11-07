#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[200, 200], positions=[30, 60, 90, 120, 150, 180, 210, 240, 270]
    )

    digiConfigFile = (
        Path(__file__).resolve().parent.parent.parent.parent
        / "Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json",
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
