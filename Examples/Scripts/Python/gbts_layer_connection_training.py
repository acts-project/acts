#!/usr/bin/env python3
import pathlib, acts, acts.root, acts.examples, acts.examples.itk
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    ParticleSelectorConfig,
    addGenParticleSelection,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import addGbtsTraining

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls()  # show acts.examples python binding calls

detector = acts.examples.itk.buildITkGeometry(geo_dir)
trackingGeometry = detector.trackingGeometry()
field = acts.root.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=1000000, numThreads=96, outputDir=str(outputDir))

if not ttbar_pu200:
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 1.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0 * u.mm, 28.0 * u.mm),
            absZ=(0.0 * u.mm, 1.0 * u.m),
            eta=(-4.0, 4.0),
            pt=(150 * u.MeV, None),
        ),
    )

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-4.0, 4.0),
        measurements=(9, None),
        removeNeutral=True,
    ),
)

addGbtsTraining(
    s,
    selectedParticles="particles_selected",
    geometryFile=pathlib.Path.cwd() / "gbts_layer_geometry.txt",
    outputConnectionTable=outputDir / "new_implementation_symm_1000k_singleMu.txt",
    probThreshold=-1.0,
    zMinTol=0.2340,
    zMaxTol=0.2340,
    rMinTol=2.5337,
    rMaxTol=2.5337,
    doSymmetrization=True,
    useOldFormatting=True,
    logLevel=acts.logging.INFO,
)

s.run()
