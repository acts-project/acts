import subprocess
import sys
import textwrap

import pytest

from helpers import geant4Enabled


def _region_creator(name, cut, volumes):
    from acts.examples.geant4 import RegionCreator

    cfg = RegionCreator.Config()
    cfg.name = name
    cfg.gammaCut = cut
    cfg.electronCut = cut
    cfg.positronCut = cut
    cfg.protonCut = 0.7
    cfg.volumes = volumes
    return RegionCreator(cfg)


@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_region_creator_construction_options():
    from acts.examples.geant4 import Geant4ConstructionOptions

    options = Geant4ConstructionOptions()
    options.regionCreators = [
        _region_creator("envelope", 1.0, ["Envelope Logic"]),
        _region_creator("layers", 0.01, ["Layer Logic"]),
    ]

    assert [r.config.name for r in options.regionCreators] == [
        "envelope",
        "layers",
    ]
    assert options.regionCreators[1].config.electronCut == pytest.approx(0.01)
    assert options.regionCreators[1].config.volumes == ["Layer Logic"]


# Geant4 keeps its regions in a process-wide store, so every simulation runs in
# its own process
_TELESCOPE_SIMULATION = textwrap.dedent("""
    import sys
    from pathlib import Path

    import acts
    import acts.examples
    from acts.examples.geant4 import Geant4ConstructionOptions, RegionCreator
    from acts.examples.simulation import (
        EtaConfig,
        ParticleConfig,
        PhiConfig,
        addGeant4,
        addParticleGun,
    )

    u = acts.UnitConstants
    outputDir = Path(sys.argv[1])
    layerCut = float(sys.argv[2])

    def regionCreator(name, cut, volumes):
        cfg = RegionCreator.Config()
        cfg.name = name
        cfg.gammaCut = cfg.electronCut = cfg.positronCut = cut
        cfg.protonCut = 0.7
        cfg.volumes = volumes
        return RegionCreator(cfg)

    detector = acts.examples.TelescopeDetector(
        acts.examples.TelescopeDetector.Config(
            bounds=[200, 200], positions=[30, 60, 90, 120, 150], stereos=[0] * 5
        )
    )
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    options = Geant4ConstructionOptions()
    options.regionCreators = [
        regionCreator("telescope_envelope", 0.7, ["Envelope Logic"]),
        regionCreator("telescope_layers", layerCut, ["Layer Logic"]),
    ]

    s = acts.examples.Sequencer(events=10, numThreads=1)
    addParticleGun(
        s,
        EtaConfig(4.0, 5.0),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(10, acts.PdgParticle.eMuon, True),
        multiplicity=5,
        rnd=rnd,
    )
    addGeant4(
        s,
        detector,
        trackingGeometry,
        field,
        rnd=rnd,
        outputDirCsv=outputDir,
        keepParticlesWithoutHits=True,
        killVolume=trackingGeometry.highestTrackingVolume,
        detectorConstructionOptions=options,
    )
    s.run()
    """)


def _count_simulated_particles(outputDir, layerCut):
    outputDir.mkdir()
    subprocess.check_call(
        [
            sys.executable,
            "-c",
            _TELESCOPE_SIMULATION,
            str(outputDir),
            str(layerCut),
        ]
    )
    files = sorted(outputDir.glob("*particles_simulated.csv"))
    assert len(files) == 10
    return sum(len(f.read_text().splitlines()) - 1 for f in files)


@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_region_creator_production_cuts(tmp_path):
    # 10 events with 50 muons each
    nPrimaries = 500

    # a small range cut in the silicon layers produces many delta rays, a range
    # cut much larger than the layer thickness produces none
    nSmallCut = _count_simulated_particles(tmp_path / "small_cut", 1e-3)
    nLargeCut = _count_simulated_particles(tmp_path / "large_cut", 1e5)

    assert nLargeCut >= nPrimaries
    assert nSmallCut > 2 * nLargeCut
