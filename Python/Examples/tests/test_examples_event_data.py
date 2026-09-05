"""Pure-Python coverage for the Examples-module EDM bindings (measurements,
tracks/track states, track containers) that needs no detector or simulation.
Mirrors Python/Core/tests/test_event_data.py, but for acts.examples types.
"""

import numpy as np
import pytest

import acts
import acts.examples as ae


def test_measurement_creation():
    meas_properties = [
        {
            "geometryId": acts.GeometryIdentifier(798),
            "indices": [0],
            "parameters": [1.0],
            "covariance": [0.1],
        },
        {
            "geometryId": acts.GeometryIdentifier(123),
            "indices": [0, 1],
            "parameters": [1.0, 2.0],
            "covariance": [0.1, 0.1],
        },
        {
            "geometryId": acts.GeometryIdentifier(456),
            "indices": [0, 1, 4],
            "parameters": [3.0, 4.0, 5.0],
            "covariance": [0.2, 0.2, 0.2],
        },
    ]

    container = acts.examples.MeasurementContainer()
    container.reserve(3)
    for meas_prop in meas_properties:
        meas = container.emplaceMeasurement(**meas_prop)

    for i in range(len(meas_properties)):
        meas = container[i]
        meas_prop = meas_properties[i]

        dim = len(meas_prop["indices"])
        assert meas.geometryId.value == meas_prop["geometryId"].value
        assert [meas_prop["indices"][i] == meas.subspaceIndices[i] for i in range(dim)]
        indices = meas_prop["indices"]
        assert [
            meas_prop["parameters"][i] == meas.fullParameters[indices[i]]
            for i in range(dim)
        ]
        assert [
            meas_prop["covariance"][i] == meas.fullCovariance[indices[i], indices[i]]
            for i in range(dim)
        ]

    assert len(container) == 3

    # Build a subset from indices 0 and 2 and verify it mirrors the container data
    subset = acts.examples.MeasurementSubset(container, [0, 2])
    assert len(subset) == 2

    # Iteration covers exactly the selected measurements in order
    subset_list = list(subset)
    assert len(subset_list) == 2
    assert subset_list[0].index == 0
    assert subset_list[1].index == 2

    # __getitem__ by subset position
    assert subset[0].index == 0
    assert subset[1].index == 2

    # getMeasurement by original-container index
    assert (
        subset.getMeasurement(0).geometryId.value
        == meas_properties[0]["geometryId"].value
    )
    assert (
        subset.getMeasurement(2).geometryId.value
        == meas_properties[2]["geometryId"].value
    )

    # Measurement data is consistent with the container entries
    for pos, orig_idx in enumerate([0, 2]):
        meas = subset[pos]
        meas_prop = meas_properties[orig_idx]
        dim = len(meas_prop["indices"])
        assert meas.geometryId.value == meas_prop["geometryId"].value
        indices = meas_prop["indices"]
        assert [meas.subspaceIndices[i] == meas_prop["indices"][i] for i in range(dim)]
        assert [
            meas.fullParameters[indices[i]] == meas_prop["parameters"][i]
            for i in range(dim)
        ]


def test_measurement_map_creation():
    from acts.examples import (
        MeasurementParticlesMap,
        MeasurementSimHitsMap,
        ParticleMeasurementsMap,
        SimBarcode,
        SimHitMeasurementsMap,
    )

    # MeasurementSimHitsMap: meas 0 → simhits {10, 11}, meas 1 → simhit {20}
    m = MeasurementSimHitsMap()
    m.insert(0, 10)
    m.insert(0, 11)  # same key — multi-map
    m.insert(1, 20)
    assert len(m) == 3

    assert 0 in m
    assert 2 not in m

    vals = m.valuesFor(0)
    assert sorted(vals) == [10, 11]
    assert m.valuesFor(1) == [20]
    assert m.valuesFor(99) == []

    pairs = list(m)
    assert len(pairs) == 3
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in pairs)

    inv = m.invert()
    assert isinstance(inv, SimHitMeasurementsMap)
    assert len(inv) == 3
    assert inv.valuesFor(10) == [0]
    assert inv.valuesFor(11) == [0]
    assert inv.valuesFor(20) == [1]

    # MeasurementParticlesMap: meas 0 came from two particles, meas 1 from one
    bc0 = SimBarcode()
    bc0.particle = 1
    bc1 = SimBarcode()
    bc1.particle = 2
    mp = MeasurementParticlesMap()
    mp.insert(0, bc0)
    mp.insert(0, bc1)  # same measurement, two particles
    mp.insert(1, bc0)
    assert len(mp) == 3

    assert mp.valuesFor(0) == [bc0, bc1]

    inv_p = mp.invert()
    assert isinstance(inv_p, ParticleMeasurementsMap)
    assert len(inv_p) == 3


# --- Track / track state / track container -------------------------------
#
# Coverage for the track/track-state container bindings that the
# simulation-driven tests (test_truth_tracking.py, test_histogram_fit_backends.py)
# don't exercise: parity between the mutable and const proxies, and the
# const <-> mutable TrackContainer round trip.


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
