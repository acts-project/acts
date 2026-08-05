"""Cross-checks of the synthetic detector layouts against real geometries.

These live here rather than next to the other Fatras tests because they need a
detector, and the ODD one needs DD4hep.
"""

import pytest

import acts
from acts.fatras import synthetic as syn

# ODD pixel volumes: 17 is the barrel, 16 and 18 are the negative and positive
# endcap
ODD_PIXEL_VOLUMES = {16, 17, 18}
# generic detector pixel volumes: 8 is the barrel, 7 and 9 the two endcaps
GENERIC_PIXEL_VOLUMES = {7, 8, 9}


def _convert(tracking_geometry, volumes=None, beam_pipe_radius=None):
    options = syn.TrackingGeometryLayoutOptions()
    if volumes is not None:
        options.setSurfaceSelector(lambda s: s.geometryId.volume in volumes)
    if beam_pipe_radius is not None:
        options.passiveBeamPipeRadius = beam_pipe_radius
    return syn.makeLayoutFromTrackingGeometry(
        tracking_geometry,
        acts.GeometryContext.dangerouslyDefaultConstruct(),
        options,
    )


def _split(layout):
    """Split a layout into its sensitive cylinders and discs.

    Returns (cylinders, discs), each entry (refCoord, minBound, maxBound).
    """
    cylinders, discs = [], []
    for surface in layout.surfaces:
        if surface.passive:
            continue
        layers = [layout.layers[i] for i in surface.layers]
        entry = (
            surface.refCoord,
            min(l.minBound for l in layers),
            max(l.maxBound for l in layers),
        )
        if surface.shape == syn.SurfaceShape.Cylinder:
            cylinders.append(entry)
        else:
            discs.append(entry)
    return cylinders, discs


@pytest.mark.odd
def test_odd_pixel_layout_matches_real_geometry(odd_detector):
    """The hardcoded ODD pixel preset has to agree with the real geometry.

    The preset exists so that the layout is available without DD4hep. This is
    what keeps it honest. Note that the preset describes where the sensitive
    silicon is, which is 1.8 mm inside the nominal layer radius the ODD XML
    declares, because a module carries its sensor at an offset.
    """
    layout = _convert(
        odd_detector.trackingGeometry(), ODD_PIXEL_VOLUMES, beam_pipe_radius=24.0
    )
    cylinders, discs = _split(layout)
    description = syn.openDataDetectorPixelDescription()

    # four barrel cylinders, and seven endcap layers per side that each resolve
    # into two rings staggered along z
    assert len(cylinders) == 4
    assert len(discs) == 2 * 14
    assert len(description.discs) == 14

    for converted, preset in zip(
        sorted(r for r, _, _ in cylinders), sorted(description.barrelRadii)
    ):
        assert converted == pytest.approx(preset, abs=0.1)

    half_length = max(max(abs(lo), abs(hi)) for _, lo, hi in cylinders)
    assert half_length == pytest.approx(max(description.barrelHalfLengthsZ), abs=0.1)

    positive = sorted(z for z, _, _ in discs if z > 0)
    for converted, preset in zip(positive, sorted(d.absZ for d in description.discs)):
        assert converted == pytest.approx(preset, abs=0.1)
    # and the endcaps are mirror images of each other
    negative = sorted(-z for z, _, _ in discs if z < 0)
    assert negative == pytest.approx(positive, abs=1e-3)

    # the rings themselves, not just the envelope around them
    preset_rings = sorted(
        (d.absZ, r.rMin, r.rMax) for d in description.discs for r in d.rings
    )
    converted_rings = sorted((z, lo, hi) for z, lo, hi in discs if z > 0)
    for converted, preset in zip(converted_rings, preset_rings):
        assert converted == pytest.approx(preset, abs=0.1)

    # the beam pipe is there, passive, and in front of the innermost layer
    passive = [s for s in layout.surfaces if s.passive]
    assert len(passive) == 1
    assert passive[0].refCoord < min(r for r, _, _ in cylinders)


def test_generic_detector_pixel_layout_matches_real_geometry(trk_geo):
    """Same check as for the ODD, on a detector that needs no external source."""
    layout = _convert(trk_geo, GENERIC_PIXEL_VOLUMES, beam_pipe_radius=19.0)
    cylinders, discs = _split(layout)
    description = syn.genericDetectorPixelDescription()

    # unlike the ODD's, a generic endcap layer is one plane of one radial band,
    # so there is no ring structure for the reduction to find
    assert len(cylinders) == 4
    assert len(discs) == 2 * 7
    assert len(description.discs) == 7
    assert all(len(d.rings) == 1 for d in description.discs)

    for converted, preset in zip(
        sorted(r for r, _, _ in cylinders), sorted(description.barrelRadii)
    ):
        assert converted == pytest.approx(preset, abs=0.1)

    half_length = max(max(abs(lo), abs(hi)) for _, lo, hi in cylinders)
    assert half_length == pytest.approx(max(description.barrelHalfLengthsZ), abs=0.1)

    positive = sorted(z for z, _, _ in discs if z > 0)
    for converted, preset in zip(positive, sorted(d.absZ for d in description.discs)):
        assert converted == pytest.approx(preset, abs=0.1)

    assert min(lo for _, lo, _ in discs) == pytest.approx(
        min(r.rMin for d in description.discs for r in d.rings), abs=0.1
    )
    assert max(hi for _, _, hi in discs) == pytest.approx(
        max(r.rMax for d in description.discs for r in d.rings), abs=0.1
    )


def test_converter_classifies_shapes(trk_geo):
    """A cylinder is thin in r and long in z, a disc the other way round."""
    layout = _convert(trk_geo, GENERIC_PIXEL_VOLUMES)
    cylinders, discs = _split(layout)

    # cylinders come out ordered from the interaction point outwards
    radii = [r for r, _, _ in cylinders]
    assert radii == sorted(radii)
    for _, lo, hi in cylinders:
        assert hi - lo > 100.0
    for refCoord, lo, hi in discs:
        assert hi - lo > 30.0
        assert abs(refCoord) > 100.0


def test_converter_handles_the_strip_system(trk_geo):
    """The converter is not pixel-specific.

    The generic detector's strip systems are volumes 12-14 and 16-18. Unlike its
    pixel endcap, a strip endcap disc really does have rings, and it interleaves
    them along z, so resolving them turns the twelve layers per side into
    eighteen discs.
    """
    layout = _convert(trk_geo, {12, 13, 14, 16, 17, 18})
    cylinders, discs = _split(layout)

    assert len(cylinders) == 4 + 2
    assert len(discs) == 2 * 18
    assert max(r for r, _, _ in cylinders) > 1000.0
    assert max(abs(z) for z, _, _ in discs) > 2900.0

    # No two of the discs of one disc cover the same radii: that would be one ring
    # counted twice, and a track crossing it would get two space points instead of
    # one. The comparison is band by band, not over the envelope of a disc, since
    # a disc carrying an inner and an outer band has another disc's band between
    # them by construction.
    max_overlap = syn.TrackingGeometryLayoutOptions().maxRingOverlap
    bands = [
        (s.refCoord, layout.layers[i].minBound, layout.layers[i].maxBound)
        for s in layout.surfaces
        if not s.passive and s.shape == syn.SurfaceShape.Disc
        for i in s.layers
    ]
    for z, lo, hi in bands:
        for other_z, other_lo, other_hi in bands:
            if z == other_z or abs(z - other_z) > 10.0:
                continue
            overlap = min(hi, other_hi) - max(lo, other_lo)
            assert overlap < max_overlap * min(hi - lo, other_hi - other_lo)


def test_generic_detector_whole_layout(trk_geo):
    """No selector means every sensitive surface is taken."""
    everything = _convert(trk_geo)
    pixel_only = _convert(trk_geo, GENERIC_PIXEL_VOLUMES)
    assert len(everything.surfaces) > len(pixel_only.surfaces)


def test_material_surfaces_are_off_by_default(trk_geo):
    """Material surfaces multiply the crossings, so they are opt-in."""
    without = _convert(trk_geo, GENERIC_PIXEL_VOLUMES)
    assert all(not s.passive for s in without.surfaces)

    options = syn.TrackingGeometryLayoutOptions()
    options.setSurfaceSelector(lambda s: s.geometryId.volume in GENERIC_PIXEL_VOLUMES)
    options.includeMaterialSurfaces = True
    with_material = syn.makeLayoutFromTrackingGeometry(
        trk_geo, acts.GeometryContext.dangerouslyDefaultConstruct(), options
    )
    assert len(with_material.surfaces) >= len(without.surfaces)


def test_layout_drives_the_generator(trk_geo):
    """A converted layout is a layout: the generator runs on it."""
    layout = _convert(trk_geo, GENERIC_PIXEL_VOLUMES, beam_pipe_radius=20.0)

    config = syn.EventConfig()
    config.generation.pileup = 2
    event = syn.generateEvent(layout, config)

    assert len(event.spacePoints) > 0
    assert len(event.particles) > 0
    summary = syn.summarize(event, 1.0)
    assert summary.primaries == config.generation.numPrimaries()
    # the beam pipe is the only material in front of the first layer, so without
    # it there would be no secondaries at all on the innermost one
    assert summary.secondaries > 0

    layer_ids = event.layerIds
    particle_ids = event.particleIds
    assert len(layer_ids) == len(event.spacePoints)
    assert len(particle_ids) == len(event.spacePoints)
    assert max(layer_ids) < len(layout.layers)
    assert max(particle_ids) < len(event.particles)
