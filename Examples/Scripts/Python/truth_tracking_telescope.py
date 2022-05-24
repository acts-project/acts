#!/usr/bin/env python3
import os
from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # detector, trackingGeometry, _ = getOpenDataDetector()
    # detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[200, 200], positions=[30, 60, 90, 120, 150, 180, 210, 240, 270]
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json",
        # / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        # "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        outputDir=Path.cwd(),
    ).run()
