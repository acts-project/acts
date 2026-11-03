import math

import pytest

import acts
import acts.examples as ex

pytestmark = pytest.mark.pypi


def run_regional_quality(tmp_path, regions, root=False, outlier=False):
    class Source(ex.IAlgorithm):
        def __init__(self):
            super().__init__("Source", acts.logging.WARNING)
            self.handles = []
            for kind, name in (
                (ex.ConstTrackContainer, "tracks"),
                (ex.MeasurementContainer, "measurements"),
                (ex.SimParticleContainer, "particles"),
                (ex.MeasurementParticlesMap, "m2p"),
                (ex.ParticleMeasurementsMap, "p2m"),
            ):
                handle = ex.WriteDataHandle(self, kind, name)
                handle.initialize(name)
                self.handles.append(handle)

        def execute(self, ctx):
            particles = ex.SimParticleContainer()
            mapping = ex.MeasurementParticlesMap()
            measurements = ex.MeasurementContainer()
            for pid in range(1, 4):
                barcode = ex.SimBarcode()
                barcode.vertexPrimary = 1
                barcode.particle = pid
                particle = ex.SimParticle(barcode, acts.PdgParticle.eMuon)
                particle.direction = acts.Vector3(1, 0, 0)
                particle.absoluteMomentum = 20
                particles.insert(particle)
                indices = range(6) if pid == 1 else [pid + 4]
                for index in indices:
                    mapping.insert(index, barcode)
            for volume, layer in [
                (10, 2),
                (10, 2),
                (10, 4),
                (10, 4),
                (20, 2),
                (20, 2),
                (30, 2),
                (40, 2),
                (40, 2),
            ]:
                measurements.emplaceMeasurement(
                    acts.GeometryIdentifier(volume=volume, layer=layer, sensitive=1),
                    [0],
                    [0.0],
                    [1.0],
                )
            # Same globally matched particle, even in the contaminated region.
            tracks = ex.TrackContainer()
            for number, indices in enumerate(
                [[0, 1, 2, 4, 6], [0, 1, 2, 3], [0, 6, 7], [8]]
            ):
                track = tracks.makeTrack()
                track.nMeasurements = len(indices)
                for index in indices:
                    state = track.appendTrackState()
                    state.typeFlags.setIsMeasurement()
                    state.uncalibratedSourceLink = ex.IndexSourceLink(
                        measurements[index].geometryId,
                        index,
                    ).toSourceLink()
                if outlier and number == 0:
                    track.nOutliers = 1
                    state = track.appendTrackState()
                    state.typeFlags.isOutlier = True
                    state.uncalibratedSourceLink = ex.IndexSourceLink(
                        measurements[5].geometryId,
                        5,
                    ).toSourceLink()
                    track.nHoles = 1
                    state = track.appendTrackState()
                    state.typeFlags.isHole = True
            values = (
                tracks.makeConst(),
                measurements,
                particles,
                mapping,
                mapping.invert(),
            )
            for handle, value in zip(self.handles, values):
                handle(ctx, value)
            return ex.ProcessCode.SUCCESS

    s = ex.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    s.addAlgorithm(Source())
    s.addAlgorithm(
        ex.TrackTruthMatcher(
            level=acts.logging.WARNING,
            inputTracks="tracks",
            inputParticles="particles",
            inputMeasurementParticlesMap="m2p",
            outputTrackParticleMatching="tmatch",
            outputParticleTrackMatching="pmatch",
            matchingRatio=0.5,
            doubleMatching=True,
        )
    )
    kwargs = dict(
        level=acts.logging.WARNING,
        inputTracks="tracks",
        inputParticles="particles",
        inputTrackParticleMatching="tmatch",
        inputParticleTrackMatching="pmatch",
        inputParticleMeasurementsMap="p2m",
    )
    baseline = ex.PythonPatternRecognitionPerformanceWriter(**kwargs)
    if regions is not None:
        kwargs.update(
            inputMeasurements="measurements", trackQualityGeometrySelections=regions
        )
    writer = ex.PythonPatternRecognitionPerformanceWriter(**kwargs)
    s.addWriter(baseline)
    s.addWriter(writer)
    if root:
        from acts.examples.root import RootPatternRecognitionPerformanceWriter

        s.addWriter(
            RootPatternRecognitionPerformanceWriter(
                **kwargs,
                filePath=str(tmp_path / "regional.root"),
            )
        )
    s.run()
    return baseline, writer


def test_global_quality_fixture(tmp_path):
    baseline, writer = run_regional_quality(tmp_path, None)
    assert baseline.stats() == writer.stats()
    assert writer.stats()["nTotalQualityTracks"] == 2
    assert writer.stats()["nTotalFakeTracks"] == 1
    assert writer.stats()["nTotalUnmatchedTracks"] == 1
    assert writer.stats()["meanCompleteness"] == pytest.approx(2 / 3)
    assert writer.stats()["meanPurity"] == pytest.approx(0.9)


def test_disabled_regional_quality_preserves_outputs(tmp_path):
    baseline, writer = run_regional_quality(tmp_path, {})
    assert baseline.stats() == writer.stats()
    assert set(baseline.histograms()) == set(writer.histograms())


def test_regional_quality_preserves_global_matches_and_empty_denominators(tmp_path):
    bh = pytest.importorskip("boost_histogram")
    regions = {
        "pixel": [acts.GeometryIdentifier(volume=10)],
        "strip": [acts.GeometryIdentifier(volume=20)],
        "contamination": [acts.GeometryIdentifier(volume=30)],
        "empty": [acts.GeometryIdentifier(volume=99)],
        "layer": [acts.GeometryIdentifier(volume=10, layer=2)],
        "sensor": [acts.GeometryIdentifier(volume=10, layer=4, sensitive=1)],
        "overlap": [
            acts.GeometryIdentifier(volume=10),
            acts.GeometryIdentifier(volume=10, layer=2),
        ],
        "all": [acts.GeometryIdentifier()],
    }
    baseline, writer = run_regional_quality(tmp_path, regions)
    for key, value in baseline.stats().items():
        assert writer.stats()[key] == value
    assert writer.stats()["nTotalQualityTracks"] == 2
    expected = {
        "pixel": ((2, 0.875), (2, 1)),
        "strip": ((2, 0.25), (1, 1)),
        "contamination": ((0, None), (1, 0)),
        "empty": ((0, None), (0, None)),
        "layer": ((2, 1), (2, 1)),
        "sensor": ((2, 0.75), (2, 1)),
        "overlap": ((2, 0.875), (2, 1)),
        "all": ((2, 2 / 3), (2, 0.9)),
    }
    for region, values in expected.items():
        stats = writer.stats()["regionalQuality"][region]
        for metric, (count, mean) in zip(("completeness", "purity"), values):
            title = metric.capitalize()
            assert stats[f"n{title}"] == count
            assert stats[f"nUndefined{title}"] == 2 - count
            if mean is None:
                assert math.isnan(stats[f"mean{title}"])
            else:
                assert stats[f"mean{title}"] == pytest.approx(mean)
            for axis in ("pT", "eta", "phi"):
                h = bh.Histogram(writer.histograms()[f"{region}_{metric}_vs_{axis}"])
                assert h.counts().sum() == count
                if count:
                    assert h[bh.loc(20 if axis == "pT" else 0)].value == pytest.approx(
                        mean
                    )


@pytest.mark.root
def test_regional_quality_root_python_parity(tmp_path):
    bh = pytest.importorskip("boost_histogram")
    ROOT = pytest.importorskip("ROOT")
    _, writer = run_regional_quality(
        tmp_path, {"pixel": [acts.GeometryIdentifier(volume=10)]}, root=True
    )
    with ROOT.TFile.Open(str(tmp_path / "regional.root")) as file:
        for metric in ("purity", "completeness"):
            key = f"pixel_{metric}_vs_pT"
            root = file.Get(key)
            native = bh.Histogram(writer.histograms()[key])
            assert root.GetBinEntries(root.FindBin(20)) == 2
            assert root.GetBinContent(root.FindBin(20)) == pytest.approx(
                native[bh.loc(20)].value
            )
        assert file.Get("pixel_nCompleteness")[0] == 2
        assert file.Get("pixel_nUndefinedCompleteness")[0] == 0


def test_regional_quality_requires_measurement_geometry():
    with pytest.raises(ValueError, match="[Mm]easurements"):
        ex.PythonPatternRecognitionPerformanceWriter(
            level=acts.logging.WARNING,
            inputTracks="tracks",
            inputParticles="particles",
            inputTrackParticleMatching="tmatch",
            inputParticleTrackMatching="pmatch",
            inputParticleMeasurementsMap="p2m",
            trackQualityGeometrySelections={
                "pixel": [acts.GeometryIdentifier(volume=10)]
            },
        )


def test_regional_outliers_follow_global_quality_semantics(tmp_path):
    baseline, writer = run_regional_quality(
        tmp_path,
        {"all": [], "strip": [acts.GeometryIdentifier(volume=20)]},
        outlier=True,
    )
    global_stats = baseline.stats()
    regions = writer.stats()["regionalQuality"]
    for metric in ("Completeness", "Purity"):
        assert regions["all"][f"mean{metric}"] == pytest.approx(
            global_stats[f"mean{metric}"]
        )
    assert regions["strip"]["meanCompleteness"] == pytest.approx(0.25)
    assert regions["strip"]["meanPurity"] == pytest.approx(0.5)
    assert regions["strip"]["nUndefinedPurity"] == 1
