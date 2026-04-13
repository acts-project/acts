import pytest

import acts
from acts.examples import Sequencer

from helpers import failure_threshold

u = acts.UnitConstants


def assert_entries(root_file, tree_name, exp=None, non_zero=False):
    __tracebackhide__ = True
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    rf = ROOT.TFile.Open(str(root_file))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert tree_name in keys
    print("Entries:", rf.Get(tree_name).GetEntries())
    if non_zero:
        assert rf.Get(tree_name).GetEntries() > 0, f"{root_file}:{tree_name}"
    if exp is not None:
        assert rf.Get(tree_name).GetEntries() == exp, f"{root_file}:{tree_name}"


def assert_has_entries(root_file, tree_name):
    __tracebackhide__ = True
    assert_entries(root_file, tree_name, non_zero=True)


@pytest.mark.parametrize("revFiltMomThresh", [0 * u.GeV, 1 * u.TeV])
def test_truth_tracking_kalman(
    tmp_path, assert_root_hash, revFiltMomThresh, detector_config
):
    root_files = [
        ("trackstates_kf.root", "trackstates", 19),
        ("tracksummary_kf.root", "tracksummary", 10),
        ("performance_kf.root", None, -1),
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    with detector_config.detector:
        from truth_tracking_kalman import runTruthTrackingKalman

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        seq = Sequencer(events=10, numThreads=1)

        runTruthTrackingKalman(
            trackingGeometry=detector_config.trackingGeometry,
            field=field,
            digiConfigFile=detector_config.digiConfigFile,
            outputDir=tmp_path,
            reverseFilteringMomThreshold=revFiltMomThresh,
            s=seq,
        )

        seq.run()

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_has_entries(fp, tn)
            assert_root_hash(fn, fp)

    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    rf = ROOT.TFile.Open(str(tmp_path / "tracksummary_kf.root"))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert "tracksummary" in keys
    for entry in rf.Get("tracksummary"):
        assert entry.hasFittedParams


def test_python_track_access(generic_detector_config, tmp_path):
    with generic_detector_config.detector:
        from truth_tracking_kalman import runTruthTrackingKalman

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        seq = Sequencer(events=100, numThreads=-1)

        runTruthTrackingKalman(
            trackingGeometry=generic_detector_config.trackingGeometry,
            field=field,
            digiConfigFile=generic_detector_config.digiConfigFile,
            outputDir=tmp_path,
            numParticles=100,
            s=seq,
        )

        import hist
        import math
        import numpy as np

        class TrackAccess(acts.examples.IAlgorithm):
            def __init__(self):
                super().__init__("TrackAccess", acts.logging.INFO)

                self.tracks = acts.examples.ReadDataHandle(
                    self, acts.examples.ConstTrackContainer, "InputTracks"
                )
                self.tracks.initialize("selected-tracks")

                self.protoTracks = acts.examples.WriteDataHandle(
                    self, acts.examples.ProtoTrackContainer, "ProtoTracks"
                )
                self.protoTracks.initialize("proto-tracks")

                self.hists = {}
                self.hists["d0"] = hist.Hist(
                    hist.axis.Regular(10, -0.1, 0.1, name="d0"),
                    label="d0",
                )
                self.hists["z0"] = hist.Hist(
                    hist.axis.Regular(10, -1, 1, name="z0"),
                    label="z0",
                )
                self.hists["phi"] = hist.Hist(
                    hist.axis.Regular(10, math.pi, -math.pi, name="phi"),
                    label="phi",
                )
                self.hists["theta"] = hist.Hist(
                    hist.axis.Regular(10, 0, math.pi, name="theta"),
                    label="theta",
                )
                self.hists["eta"] = hist.Hist(
                    hist.axis.Regular(10, -5, 5, name="eta"),
                    label="eta",
                )
                self.hists["qop"] = hist.Hist(
                    hist.axis.Regular(10, -1, 1, name="qop"),
                    label="qop",
                )
                self.hists["nTracks"] = hist.Hist(
                    hist.axis.Regular(10, 9800, 10200, name="nTracks"),
                    label="nTracks",
                )

            def execute(self, context):
                self.logger.info("Track access")

                tracks = self.tracks(context.eventStore)
                assert isinstance(tracks, acts.examples.ConstTrackContainer)

                self.logger.info("Tracks: {}", len(tracks))

                self.hists["d0"].fill(tracks.parameters[:, 0])
                self.hists["z0"].fill(tracks.parameters[:, 1])
                self.hists["phi"].fill(tracks.parameters[:, 2])
                self.hists["theta"].fill(tracks.parameters[:, 3])
                self.hists["eta"].fill(-np.log(np.tan(tracks.parameters[:, 3] / 2)))
                self.hists["qop"].fill(tracks.parameters[:, 4])
                self.hists["nTracks"].fill(len(tracks))
                # for track in tracks:
                #     self.logger.info("Track: {}", track)
                #     self.logger.info("Track index: {}", track.index)
                #     self.logger.info("Track tip index: {}", track.tipIndex)
                #     self.logger.info("Track stem index: {}", track.stemIndex)
                #     self.logger.info(
                #         "Track reference surface: {}", track.referenceSurface
                #     )
                #     self.logger.info(
                #         "Track has reference surface: {}", track.hasReferenceSurface
                #     )
                #     self.logger.info("Track parameters: {}", track.parameters)
                #     self.logger.info("Track covariance: {}", track.covariance)
                #     self.logger.info(
                #         "Track particle hypothesis: {}", track.particleHypothesis
                #     )

                wb = context.eventStore
                assert not wb.exists("proto-tracks")
                myProtoTracks = acts.examples.ProtoTrackContainer()
                self.protoTracks(context, myProtoTracks)
                assert wb.exists("proto-tracks")

                return acts.examples.ProcessCode.SUCCESS

            def finalize(self):
                for h in self.hists.values():
                    print(h.label)
                    print(h)
                return acts.examples.ProcessCode.SUCCESS

        seq.addAlgorithm(TrackAccess())

        with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
            seq.run()


def test_python_track_state_access(generic_detector_config, tmp_path):
    with generic_detector_config.detector:
        from truth_tracking_kalman import runTruthTrackingKalman

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        seq = Sequencer(events=10, numThreads=1)

        runTruthTrackingKalman(
            trackingGeometry=generic_detector_config.trackingGeometry,
            field=field,
            digiConfigFile=generic_detector_config.digiConfigFile,
            outputDir=tmp_path,
            numParticles=10,
            linkForward=True,
            s=seq,
        )

        class TrackStateAccess(acts.examples.IAlgorithm):
            def __init__(self):
                super().__init__("TrackStateAccess", acts.logging.INFO)

                self.tracks = acts.examples.ReadDataHandle(
                    self, acts.examples.ConstTrackContainer, "InputTracks"
                )
                self.tracks.initialize("selected-tracks")

            def execute(self, context):
                tracks = self.tracks(context.eventStore)

                import numpy as np

                for track in tracks:
                    params = track.parameters
                    assert isinstance(params, acts.BoundVector)
                    cov = track.covariance
                    assert isinstance(cov, acts.BoundMatrix)
                    # parameters from per-proxy accessor must match the
                    # bulk numpy array at the same index
                    assert np.allclose(
                        [params[i] for i in range(6)],
                        tracks.parameters[track.index],
                    )

                    n_meas_from_summary = track.nMeasurements
                    n_meas_counted = 0
                    n_holes_counted = 0

                    for state in track.trackStatesReversed:
                        assert isinstance(state, acts.examples.ConstTrackStateProxy)

                        flags = state.typeFlags
                        # every state must be at least one of the known types
                        assert (
                            flags.isMeasurement
                            or flags.isOutlier
                            or flags.isHole
                            or flags.hasMaterial
                        )

                        if flags.isMeasurement:
                            n_meas_counted += 1

                        if flags.isHole:
                            n_holes_counted += 1

                        if state.hasPredicted:
                            pred = state.predicted
                            assert isinstance(pred, acts.BoundVector)

                        if state.hasFiltered:
                            filt = state.filtered
                            assert isinstance(filt, acts.BoundVector)

                        if state.hasSmoothed:
                            smth = state.smoothed
                            assert isinstance(smth, acts.BoundVector)

                    assert n_meas_counted == n_meas_from_summary
                    assert n_holes_counted == track.nHoles

                    assert track.isForwardLinked

                    rev_predicted = [
                        state.predicted
                        for state in track.trackStatesReversed
                        if state.hasPredicted
                    ]
                    fwd_predicted = [
                        state.predicted
                        for state in track.trackStates
                        if state.hasPredicted
                    ]
                    assert len(fwd_predicted) == len(rev_predicted)
                    for fwd, rev in zip(fwd_predicted, reversed(rev_predicted)):
                        assert all(fwd[i] == pytest.approx(rev[i]) for i in range(6))

                return acts.examples.ProcessCode.SUCCESS

        seq.addAlgorithm(TrackStateAccess())

        with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
            seq.run()


@pytest.mark.skip(reason="Needs updating after converter became unnecessary")
def test_python_space_point_access(generic_detector_config, tmp_path):
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
    )

    with generic_detector_config.detector:
        s = acts.examples.Sequencer(
            events=100, numThreads=1, logLevel=acts.logging.INFO
        )
        trackingGeometry = generic_detector_config.trackingGeometry
        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
        digiConfigFile = generic_detector_config.digiConfigFile

        rnd = acts.examples.RandomNumbers(seed=42)

        logger = acts.getDefaultLogger("Truth tracking example", acts.logging.INFO)

        addParticleGun(
            s,
            ParticleConfig(num=10, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
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
            seedingAlgorithm=SeedingAlgorithm.GridTriplet,
            geoSelectionConfigFile=generic_detector_config.geometrySelection,
            particleHypothesis=acts.ParticleHypothesis.muon,
            initialSigmas=[
                1 * u.mm,
                1 * u.mm,
                1 * u.degree,
                1 * u.degree,
                0 / u.GeV,
                1 * u.ns,
            ],
            initialSigmaQoverPt=0.1 / u.GeV,
            initialSigmaPtRel=0.1,
            initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        )

        spConverter = acts.examples.SpacePointConverter(
            inputSpacePoints="spacepoints",
            outputSpacePoints="spacepoints2",
            logger=logger,
        )
        s.addAlgorithm(spConverter)

        class SpacePointAccess(acts.examples.IAlgorithm):

            def __init__(self):
                super().__init__("SpacePointAccess", acts.logging.INFO)

                self.spacePoints = acts.examples.ReadDataHandle(
                    self, acts.SpacePointContainer2, "InputSpacePoints"
                )
                self.spacePoints.initialize("spacepoints2")

            def execute(self, context: acts.examples.AlgorithmContext):
                self.logger.info("Space point access")
                spacePoints: acts.SpacePointContainer2 = self.spacePoints(
                    context.eventStore
                )

                for sp in spacePoints:
                    self.logger.info("Space point: {}", sp.x)
                    self.logger.info("Space point: {}", sp.y)
                    self.logger.info("Space point: {}", sp.z)
                    self.logger.info("Space point: {}", sp.r)
                    self.logger.info("Space point: {}", sp.xy)
                    self.logger.info("Space point: {}", sp.xy[:, 0])
                    self.logger.info("Space point: {}", sp.xy[:, 1])
                    self.logger.info("Space point: {}", sp.xy[:, 2])
                    self.logger.info("Space point: {}", sp.xy[:, 3])

                return acts.examples.ProcessCode.SUCCESS

        s.addAlgorithm(SpacePointAccess())

        s.run()


def test_truth_tracking_gsf(tmp_path, assert_root_hash, detector_config):
    from truth_tracking_gsf import runTruthTrackingGsf

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = Sequencer(
        events=10,
        numThreads=1,
    )

    root_files = [
        ("trackstates_gsf.root", "trackstates"),
        ("tracksummary_gsf.root", "tracksummary"),
    ]

    for fn, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    with detector_config.detector:
        runTruthTrackingGsf(
            trackingGeometry=detector_config.trackingGeometry,
            decorators=detector_config.decorators,
            field=field,
            digiConfigFile=detector_config.digiConfigFile,
            outputDir=tmp_path,
            s=seq,
        )

        # See https://github.com/acts-project/acts/issues/1300
        with failure_threshold(acts.logging.FATAL):
            seq.run()

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_root_hash(fn, fp)


def test_refitting(tmp_path, detector_config, assert_root_hash):
    from truth_tracking_gsf_refitting import runRefittingGsf

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = Sequencer(
        events=10,
        numThreads=1,
    )

    with detector_config.detector:
        # Only check if it runs without errors right known
        # Changes in fitter behaviour should be caught by other tests
        runRefittingGsf(
            trackingGeometry=detector_config.trackingGeometry,
            field=field,
            digiConfigFile=detector_config.digiConfigFile,
            outputDir=tmp_path,
            s=seq,
        ).run()

    root_files = [
        ("trackstates_gsf_refit.root", "trackstates"),
        ("tracksummary_gsf_refit.root", "tracksummary"),
    ]

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_root_hash(fn, fp)
