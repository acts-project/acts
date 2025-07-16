#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:
    detector = acts.examples.TelescopeDetector(
        bounds=[200, 200],
        positions=[30, 60, 90, 120, 150, 180, 210, 240, 270],
        stereos=[0] * 9,
    )
    trackingGeometry = detector.trackingGeometry()

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=srcdir / "Examples/Configs/telescope-digi-smearing-config.json",
        outputDir=Path.cwd(),
    ).run()
