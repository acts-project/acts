import pytest
from helpers import failure_threshold

import acts
import acts.examples as ex

pytestmark = pytest.mark.pypi


@pytest.mark.parametrize("geometry_mode", ["provided", "omitted", "unknown_id"])
def test_proto_track_state_surfaces(geometry_mode):
    detector = ex.GenericDetector()
    geometry = detector.trackingGeometry()
    surfaces = []
    geometry.visitSurfaces(surfaces.append)
    geo_id = surfaces[0].geometryId
    if geometry_mode == "unknown_id":
        geo_id = acts.GeometryIdentifier(volume=255, sensitive=1)

    class Source(ex.IAlgorithm):
        def __init__(self):
            super().__init__("Source", acts.logging.WARNING)
            self.measurements = ex.MeasurementContainer()
            self.measurements.emplaceMeasurement(geo_id, [0], [0.0], [1.0])
            self.output_measurements = ex.WriteDataHandle(
                self, ex.MeasurementSubset, "Measurements"
            )
            self.output_measurements.initialize("measurements")
            self.output_protos = ex.WriteDataHandle(
                self, ex.ProtoTrackContainer, "Protos"
            )
            self.output_protos.initialize("protos")

        def execute(self, ctx):
            self.output_measurements(ctx, ex.MeasurementSubset(self.measurements, [0]))
            self.output_protos(ctx, ex.ProtoTrackContainer([ex.ProtoTrack([0])]))
            return ex.ProcessCode.SUCCESS

    class Check(ex.IAlgorithm):
        def __init__(self):
            super().__init__("Check", acts.logging.WARNING)
            self.input_tracks = ex.ReadDataHandle(
                self, ex.ConstTrackContainer, "Tracks"
            )
            self.input_tracks.initialize("tracks")
            self.checked = False

        def execute(self, ctx):
            tracks = self.input_tracks(ctx.eventStore)
            assert len(tracks) == 1
            track = tracks.getTrack(0)
            assert track.nMeasurements == 1
            state = next(iter(track.trackStatesReversed))
            assert state.hasReferenceSurface == (geometry_mode == "provided")
            if geometry_mode == "provided":
                assert state.referenceSurface.geometryId.value == geo_id.value
            self.checked = True
            return ex.ProcessCode.SUCCESS

    sequencer = ex.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    sequencer.addAlgorithm(Source())
    kwargs = {} if geometry_mode == "omitted" else {"trackingGeometry": geometry}
    sequencer.addAlgorithm(
        ex.ProtoTracksToTracks(
            level=acts.logging.WARNING,
            inputProtoTracks="protos",
            inputMeasurements="measurements",
            outputTracks="tracks",
            **kwargs,
        )
    )
    check = Check()
    sequencer.addAlgorithm(check)
    if geometry_mode == "unknown_id":
        with (
            failure_threshold(acts.logging.MAX),
            pytest.raises(RuntimeError, match="No surface found"),
        ):
            sequencer.run()
    else:
        sequencer.run()
        assert check.checked
