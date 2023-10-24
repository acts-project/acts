from pathlib import Path

import pytest

import acts.examples

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)

from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TruthEstimatedSeedingAlgorithmConfigArg,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorConfig,
)

from helpers import failure_threshold

u = acts.UnitConstants


@pytest.mark.parametrize(
    "seeding_algorithm",
    [
        SeedingAlgorithm.Default,
        SeedingAlgorithm.TruthSmeared,
        SeedingAlgorithm.TruthEstimated,
        SeedingAlgorithm.Orthogonal,
    ],
)
def test_ckf_tracking(seeding_algorithm, physmon: "Physmon"):
    s = acts.examples.Sequencer(
        events=10,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    for d in physmon.decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
        ),
        multiplicity=50,
        rnd=rnd,
    )

    addFatras(
        s,
        physmon.trackingGeometry,
        physmon.field,
        enableInteractions=True,
        rnd=rnd,
    )

    addDigitization(
        s,
        physmon.trackingGeometry,
        physmon.field,
        digiConfigFile=physmon.digiConfig,
        rnd=rnd,
    )

    addSeeding(
        s,
        physmon.trackingGeometry,
        physmon.field,
        TruthSeedRanges(pt=(500 * u.MeV, None), nHits=(9, None)),
        ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
        SeedFinderConfigArg(
            r=(33 * u.mm, 200 * u.mm),
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=5,
            radLengthPerSeed=0.1,
            minPt=500 * u.MeV,
            impactMax=3 * u.mm,
        ),
        SeedFinderOptionsArg(bFieldInZ=2 * u.T),
        TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(10.0 * u.mm, None)),
        seedingAlgorithm=seeding_algorithm,
        geoSelectionConfigFile=physmon.geoSel,
        rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
        outputDirRoot=physmon.tmp_path,
    )

    addCKFTracks(
        s,
        physmon.trackingGeometry,
        physmon.field,
        TrackSelectorConfig(
            pt=(500 * u.MeV, None),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=6,
        ),
        outputDirRoot=physmon.tmp_path,
    )

    stems = [
        "performance_ckf",
        "tracksummary_ckf",
    ]

    associatedParticles = "particles_input"
    if seeding_algorithm in [SeedingAlgorithm.Default, SeedingAlgorithm.Orthogonal]:
        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(maximumSharedHits=3),
            outputDirRoot=physmon.tmp_path,
        )

        associatedParticles = None

        stems += ["performance_ambi"]

    addVertexFitting(
        s,
        physmon.field,
        associatedParticles=associatedParticles,
        outputProtoVertices="ivf_protovertices",
        outputVertices="ivf_fittedVertices",
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=physmon.tmp_path / "ivf",
    )

    addVertexFitting(
        s,
        physmon.field,
        associatedParticles=associatedParticles,
        outputProtoVertices="amvf_protovertices",
        outputVertices="amvf_fittedVertices",
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=physmon.tmp_path / "amvf",
    )

    with failure_threshold(acts.logging.FATAL):
        s.run()
    del s

    # @TODO: Add plotting into ROOT file for vertexing and ckf tracksummary

    for vertexing in ["ivf", "amvf"]:
        target = f"performance_{vertexing}.root"
        physmon.add_output_file(
            f"{vertexing}/performance_vertexing.root",
            rename=target,
        )
        physmon.histogram_comparison(target, f"Performance {vertexing}")

    if seeding_algorithm in [
        SeedingAlgorithm.Default,
        SeedingAlgorithm.Orthogonal,
        SeedingAlgorithm.TruthEstimated,
    ]:
        stems += ["performance_seeding"]

    for stem in stems:
        target = f"{stem}.root"
        physmon.add_output_file(f"{stem}.root", rename=target)
        physmon.histogram_comparison(target, f"Performance {physmon.name} {stem}")

    #  ] + (
    #  ["performance_seeding", "performance_ambi"]
    #  if label in ["seeded", "orthogonal"]
    #  else ["performance_seeding"]
    #  if label == "truth_estimated"
    #  else []
    #  ):
    #  perf_file = physmon.tmp_path / f"{stem}.root"
    #  assert perf_file.exists(), "Performance file not found"
    #  shutil.copy(perf_file, physmon.outdir / f"{stem}_{label}.root")
