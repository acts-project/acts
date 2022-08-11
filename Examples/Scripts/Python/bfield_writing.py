#!/usr/bin/env python3
from pathlib import Path

import acts
import acts.examples

u = acts.UnitConstants


def runBFieldWriting(outputDir: Path, rewrites: int = 0):
    solenoid = acts.SolenoidBField(
        radius=1200 * u.mm, length=6000 * u.mm, bMagCenter=2 * u.T, nCoils=1194
    )
    field = acts.solenoidFieldMap(
        rlim=(0, 1200 * u.mm),
        zlim=(-5000 * u.mm, 5000 * u.mm),
        nbins=(10, 10),
        field=solenoid,
    )

    print("Solenoid ready")

    cfg = acts.examples.RootBFieldWriter.Config()
    cfg.bField = field
    cfg.gridType = acts.examples.RootBFieldWriter.GridType.rz
    cfg.fileName = str(outputDir / "solenoid.root")
    cfg.treeName = "solenoid"

    acts.examples.RootBFieldWriter.run(cfg, acts.logging.VERBOSE)

    for i in range(rewrites):
        print(f"Now read back {cfg.fileName}")

        field2 = acts.examples.MagneticFieldMapRz(cfg.fileName, tree="solenoid")

        cfg2 = acts.examples.RootBFieldWriter.Config()
        cfg2.bField = field2
        cfg2.gridType = acts.examples.RootBFieldWriter.GridType.rz
        cfg2.fileName = str(outputDir / f"solenoid{i+2}.root")
        cfg2.treeName = "solenoid"

        acts.examples.RootBFieldWriter.run(cfg2, acts.logging.VERBOSE)

        cfg = cfg2

    print("Done")


if "__main__" == __name__:
    runBFieldWriting(Path.cwd())
