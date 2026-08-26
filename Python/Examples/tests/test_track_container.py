"""Coverage for the track/track-state container Python bindings that the
existing simulation-driven tests (test_truth_tracking.py,
test_histogram_fit_backends.py) don't exercise: parity between the mutable
and const proxies, and the const <-> mutable TrackContainer round trip. Pure
Python, no detector or simulation needed.
"""

import numpy as np
import pytest

import acts
import acts.examples as ae


def _fill_track(track):
    surface = acts.Surface.createPerigee(acts.Vector3(0, 0, 0))
    geo_id = acts.GeometryIdentifier()

    track.referenceSurface = surface
    track.parameters = acts.BoundVector(0.1, 0.2, 0.3, 1.4, 0.01, 0.0)
    track.covariance = acts.BoundMatrix.Identity()
    track.particleHypothesis = acts.ParticleHypothesis.muon
    track.nMeasurements = 3
    track.nHoles = 1
    track.nOutliers = 2
    track.nSharedHits = 1
    track.chi2 = 4.5
    track.nDoF = 5

    state = track.appendTrackState(ae.TrackStatePropMask.All)
    state.typeFlags.setIsMeasurement()
    state.referenceSurface = surface
    state.uncalibratedSourceLink = ae.IndexSourceLink(geo_id, 0).toSourceLink()
    state.predicted = acts.BoundVector(0.1, 0.2, 0.3, 1.4, 0.01, 0.0)
    state.predictedCovariance = acts.BoundMatrix.Identity()
    state.filtered = acts.BoundVector(0.15, 0.25, 0.3, 1.4, 0.01, 0.0)
    state.filteredCovariance = acts.BoundMatrix.Identity()
    state.smoothed = acts.BoundVector(0.12, 0.22, 0.3, 1.4, 0.01, 0.0)
    state.smoothedCovariance = acts.BoundMatrix.Identity()
    state.jacobian = acts.BoundMatrix.Identity()
    state.chi2 = 1.5
    state.pathLength = 12.0
    state.allocateCalibrated(2)
    state.effectiveCalibrated = [1.0, 2.0]
    state.effectiveCalibratedCovariance = [[1.0, 0.0], [0.0, 1.0]]
    state.setProjectorSubspaceIndices([0, 1])

    return geo_id, state


def _check_track(track, state, *, linked):
    assert track.index == 0
    assert track.tipIndex == state.index
    assert track.stemIndex == (state.index if linked else ae.kTrackIndexInvalid)
    assert track.hasReferenceSurface
    assert track.parameters[0] == pytest.approx(0.1)
    assert track.covariance[0, 0] == pytest.approx(1.0)
    assert track.particleHypothesis.absolutePdg == acts.PdgParticle.eMuon
    assert track.nMeasurements == 3
    assert track.nHoles == 1
    assert track.nOutliers == 2
    assert track.nSharedHits == 1
    assert track.chi2 == pytest.approx(4.5)
    assert track.nDoF == 5
    assert track.nTrackStates == 1
    assert track.isForwardLinked is linked


def _check_state(state, geo_id, *, is_const):
    assert state.hasReferenceSurface
    assert state.referenceSurface is not None
    assert state.hasPredicted
    assert state.hasFiltered
    assert state.hasSmoothed
    assert state.hasJacobian
    assert state.hasProjector
    assert state.hasUncalibratedSourceLink
    assert state.hasCalibrated
    assert state.calibratedSize == 2
    assert state.effectiveCalibrated == pytest.approx([1.0, 2.0])
    assert np.allclose(state.effectiveCalibratedCovariance, [[1.0, 0.0], [0.0, 1.0]])
    assert list(state.projectorSubspaceIndices) == [0, 1]
    assert state.chi2 == pytest.approx(1.5)
    assert state.pathLength == pytest.approx(12.0)
    assert state.predicted[0] == pytest.approx(0.1)
    assert state.filtered[0] == pytest.approx(0.15)
    assert state.smoothed[0] == pytest.approx(0.12)
    assert all(
        state.jacobian[i, j] == pytest.approx(1.0 if i == j else 0.0)
        for i in range(6)
        for j in range(6)
    )
    assert state.parameters[0] == pytest.approx(state.smoothed[0])
    assert state.typeFlags.isMeasurement

    # the originally reported gap: uncalibratedSourceLink and referenceSurface
    # must be readable on the const proxy, not just the mutable one.
    isl = ae.IndexSourceLink.fromSourceLink(state.uncalibratedSourceLink)
    assert isl.geometryId() == geo_id

    proxy_type = ae.ConstTrackStateProxy if is_const else ae.TrackStateProxy
    assert isinstance(state, proxy_type)


def test_mutable_track_proxy():
    tc = ae.TrackContainer()
    track = tc.makeTrack()
    geo_id, state = _fill_track(track)

    _check_track(track, state, linked=False)
    _check_state(state, geo_id, is_const=False)

    # trackStates()/trackStatesReversed() are only meaningful once linked
    track.linkForward()
    fwd = list(track.trackStates)
    rev = list(track.trackStatesReversed)
    assert len(fwd) == len(rev) == 1
    assert fwd[0].index == rev[0].index == state.index

    assert track.hasColumn("doesNotExist") is False


def test_const_mutable_parity():
    """The const and mutable track/track-state proxies must read back
    identical values for every property they share -- the parity the
    templated binder is supposed to guarantee by construction."""
    tc = ae.TrackContainer()
    track = tc.makeTrack()
    geo_id, _ = _fill_track(track)
    track.linkForward()

    const_tc = tc.makeConst()
    const_track = const_tc[0]
    const_state = next(iter(const_track.trackStatesReversed))

    _check_track(const_track, const_state, linked=True)
    _check_state(const_state, geo_id, is_const=True)

    assert isinstance(const_track, ae.ConstTrackProxy)
    assert const_track.hasColumn("doesNotExist") is False


def test_const_to_mutable_round_trip():
    tc = ae.TrackContainer()
    track = tc.makeTrack()
    _fill_track(track)
    track.linkForward()

    const_tc = tc.makeConst()

    # acts.examples.TrackContainer(const_tc) and const_tc.makeMutable() are
    # both independent, fully mutable copies.
    for mutable_copy in (ae.TrackContainer(const_tc), const_tc.makeMutable()):
        assert len(mutable_copy) == 1
        copy_track = mutable_copy[0]
        assert copy_track.chi2 == pytest.approx(4.5)
        assert copy_track.nMeasurements == 3

        copy_track.chi2 = 99.0
        assert copy_track.chi2 == pytest.approx(99.0)
        # the source const container is untouched by mutating the copy
        assert const_tc[0].chi2 == pytest.approx(4.5)


def test_track_container_iteration_and_getitem():
    tc = ae.TrackContainer()
    for i in range(3):
        t = tc.makeTrack()
        t.chi2 = float(i)

    assert len(tc) == 3
    assert [t.chi2 for t in tc] == [0.0, 1.0, 2.0]
    assert tc[1].chi2 == pytest.approx(1.0)
    assert tc.getTrack(2).chi2 == pytest.approx(2.0)

    const_tc = tc.makeConst()
    assert len(const_tc) == 3
    assert [t.chi2 for t in const_tc] == [0.0, 1.0, 2.0]
    assert const_tc[1].chi2 == pytest.approx(1.0)
    assert const_tc.getTrack(2).chi2 == pytest.approx(2.0)


def test_track_container_soa_numpy_views():
    tc = ae.TrackContainer()
    for i in range(3):
        t = tc.makeTrack()
        t.parameters = acts.BoundVector(float(i), 0.0, 0.0, 1.4, 0.01, 0.0)
        t.chi2 = float(i) + 0.5
        t.nOutliers = i
        t.nSharedHits = 2 * i

    const_tc = tc.makeConst()

    for i, track in enumerate(const_tc):
        assert const_tc.parameters[i, 0] == pytest.approx(track.parameters[0])
        assert const_tc.chi2[i] == pytest.approx(track.chi2)
        assert const_tc.nOutliers[i] == track.nOutliers
        assert const_tc.nSharedHits[i] == track.nSharedHits


def test_ensure_dynamic_columns_and_copy_from():
    tc = ae.TrackContainer()
    track = tc.makeTrack()
    _fill_track(track)
    track.linkForward()
    const_tc = tc.makeConst()

    dst = ae.TrackContainer()
    dst.ensureDynamicColumns(const_tc)
    dst_track = dst.makeTrack()
    dst_track.copyFrom(const_tc[0])

    assert dst_track.chi2 == pytest.approx(4.5)
    assert dst_track.nMeasurements == 3
    assert len(list(dst_track.trackStates)) == 1


def test_any_proxy():
    """AnyMutableTrackProxy/AnyConstTrackProxy/AnyMutableTrackStateProxy/
    AnyConstTrackStateProxy are type-erased views onto the same underlying
    storage as the concrete proxy they were constructed from -- reading them
    must agree, and mutating through the Any* handle must be visible back on
    the original proxy.

    Note: an Any*Proxy constructed from a mutable TrackProxy/TrackStateProxy
    must not outlive a TrackContainer.makeConst() call on its container --
    that moves the backing storage out, and the Any*Proxy would dangle just
    like the original proxy would.
    """
    tc = ae.TrackContainer()
    track = tc.makeTrack()
    geo_id, state = _fill_track(track)

    any_mut_track = ae.AnyMutableTrackProxy(track)
    assert any_mut_track.chi2 == pytest.approx(4.5)
    assert any_mut_track.tipIndex == track.tipIndex
    assert any_mut_track.parameters[0] == pytest.approx(track.parameters[0])
    any_mut_track.chi2 = 99.0
    assert track.chi2 == pytest.approx(99.0)  # same storage
    track.chi2 = 4.5

    any_mut_state = ae.AnyMutableTrackStateProxy(state)
    assert any_mut_state.chi2 == pytest.approx(1.5)
    assert any_mut_state.predicted[0] == pytest.approx(state.predicted[0])
    isl = ae.IndexSourceLink.fromSourceLink(any_mut_state.uncalibratedSourceLink)
    assert isl.geometryId() == geo_id
    assert any_mut_state.effectiveCalibrated == pytest.approx([1.0, 2.0])
    any_mut_state.chi2 = 88.0
    assert state.chi2 == pytest.approx(88.0)  # same storage
    state.chi2 = 1.5

    const_tc = tc.makeConst()
    const_track = const_tc[0]
    const_state = next(iter(const_track.trackStatesReversed))

    any_const_track = ae.AnyConstTrackProxy(const_track)
    assert any_const_track.chi2 == pytest.approx(4.5)

    any_const_state = ae.AnyConstTrackStateProxy(const_state)
    assert any_const_state.chi2 == pytest.approx(1.5)
