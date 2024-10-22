#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[200, 200],
        positions=[30, 60, 90, 120, 150, 180, 210, 240, 270],
        stereos=[0] * 9,
    )
    
    pos=[10000,  10500, 11000, 19500, 20000,20500]
    n_layers = len(pos)
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[500, 1500], positions=pos, binValue=0,thickness=4,stereos=[0] * n_layers,
    )
    
    inputParticlePath = Path("")
    if not inputParticlePath.exists():
        inputParticlePath = None


    # inputParticlePath = Path("/data/atlassmallfiles/users/salin/Acts_x/GEN/HepMC_Root/Root_DarkPhoton_X_e_mu/Particles_DarkPhoton_m0.6457_mu_mu.root")
    # if not inputParticlePath.exists():
    #     inputParticlePath = None

    outputDir = Path(f"./Output_ttk/")
    outputDir.mkdir(parents=True, exist_ok=True)


    # field = acts.ConstantBField(acts.Vector3(0* u.T, 0, 1 * u.T))
    field = acts.RestrictedBField(acts.Vector3(0* u.T, 0, 1.0 * u.T))

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    #field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json",
        outputDir=Path.cwd(),
    ).run()
