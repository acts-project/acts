if False:
    from pathlib import Path

    import pytest

    import acts
    import acts.examples
    from acts.examples import Sequencer, GenericDetector, RootParticleWriter

    from helpers import (
        dd4hepEnabled,
        pythia8Enabled,
        AssertCollectionExistsAlg,
        assert_csv_output,
        assert_entries,
        assert_has_entries,
    )

    from acts.examples.odd import getOpenDataDetector
    from common import getOpenDataDetectorDirectory

    u = acts.UnitConstants

    @pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
    @pytest.mark.skipif(not pythia8Enabled, reason="Pythia8 not set up")
    @pytest.mark.slow
    @pytest.mark.odd
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_vertex_fitting(tmp_path):
        detector, trackingGeometry, decorators = getOpenDataDetector(
            getOpenDataDetectorDirectory()
        )

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        from vertex_fitting import runVertexFitting, VertexFinder

        s = Sequencer(events=100)

        runVertexFitting(
            field,
            vertexFinder=VertexFinder.Truth,
            outputDir=tmp_path,
            s=s,
        )

        alg = AssertCollectionExistsAlg(["fittedVertices"], name="check_alg")
        s.addAlgorithm(alg)

        s.run()
        assert alg.events_seen == s.config.events

    @pytest.mark.parametrize(
        "finder,inputTracks,entries",
        [
            ("Truth", False, 100),
            # ("Truth", True, 0), # this combination seems to be not working
            ("Iterative", False, 100),
            ("Iterative", True, 100),
            ("AMVF", False, 100),
            ("AMVF", True, 100),
        ],
    )
    @pytest.mark.filterwarnings("ignore::UserWarning")
    @pytest.mark.flaky(reruns=2)
    def test_vertex_fitting_reading(
        tmp_path, ptcl_gun, rng, finder, inputTracks, entries, assert_root_hash
    ):
        ptcl_file = tmp_path / "particles.root"

        detector, trackingGeometry, decorators = GenericDetector.create()
        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        from vertex_fitting import runVertexFitting, VertexFinder

        inputTrackSummary = None
        if inputTracks:
            from truth_tracking_kalman import runTruthTrackingKalman

            s2 = Sequencer(numThreads=1, events=100)
            runTruthTrackingKalman(
                trackingGeometry,
                field,
                digiConfigFile=Path(
                    Path(__file__).parent.parent.parent.parent
                    / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
                ),
                outputDir=tmp_path,
                s=s2,
            )
            s2.run()
            del s2
            inputTrackSummary = tmp_path / "tracksummary_fitter.root"
            assert inputTrackSummary.exists()
            assert ptcl_file.exists()
        else:
            s0 = Sequencer(events=100, numThreads=1)
            evGen = ptcl_gun(s0)
            s0.addWriter(
                RootParticleWriter(
                    level=acts.logging.INFO,
                    inputParticles=evGen.config.outputParticles,
                    filePath=str(ptcl_file),
                )
            )
            s0.run()
            del s0

            assert ptcl_file.exists()

        finder = VertexFinder[finder]

        s3 = Sequencer(numThreads=1)

        runVertexFitting(
            field,
            inputParticlePath=ptcl_file,
            inputTrackSummary=inputTrackSummary,
            outputDir=tmp_path,
            vertexFinder=finder,
            s=s3,
        )

        alg = AssertCollectionExistsAlg(["fittedVertices"], name="check_alg")
        s3.addAlgorithm(alg)

        s3.run()

        vertexing_file = tmp_path / "performance_vertexing.root"
        assert vertexing_file.exists()

        assert_entries(vertexing_file, "vertexing", entries)
        assert_root_hash(vertexing_file.name, vertexing_file)
