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


def _options(volumes=None, beam_pipe_radius=None):
    options = syn.TrackingGeometryLayoutOptions()
    if volumes is not None:
        options.setSurfaceSelector(lambda s: s.geometryId.volume in volumes)
    if beam_pipe_radius is not None:
        options.passiveBeamPipeRadius = beam_pipe_radius
    return options


def _describe(tracking_geometry, volumes=None, beam_pipe_radius=None):
    return syn.makeDescriptionFromTrackingGeometry(
        tracking_geometry,
        acts.GeometryContext.dangerouslyDefaultConstruct(),
        _options(volumes, beam_pipe_radius),
    )


def _convert(tracking_geometry, volumes=None, beam_pipe_radius=None):
    return syn.makeLayout(_describe(tracking_geometry, volumes, beam_pipe_radius))


def _cylinders(description):
    """Every cylinder of a description, whichever subsystem and barrel it is in."""
    return [
        cylinder
        for subsystem in description.subsystems
        for barrel in subsystem.barrels
        for cylinder in barrel.cylinders
    ]


def _one_subsystem(subsystem):
    """A one-subsystem description, so that the helpers above apply to it."""
    description = syn.DetectorDescription()
    description.subsystems = [subsystem]
    return description


def _discs(description):
    """Every disc of a description, and how many sides its endcap is built on."""
    return [
        (disc, 2 if endcap.placement == syn.EndcapPlacement.Mirrored else 1)
        for subsystem in description.subsystems
        for endcap in subsystem.endcaps
        for disc in endcap.discs
    ]


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
    """The shipped ODD pixel description has to agree with the real geometry.

    The files exist so that the detector is available without DD4hep. This is
    what keeps them honest. Note that they describe where the sensitive silicon
    is, which is 1.8 mm inside the nominal layer radius the ODD XML declares,
    because a module carries its sensor at an offset.

    The comparison is between the two *layouts* rather than the two
    descriptions: the reduction measures each side of the detector separately
    and the shipped file mirrors one, which is the same detector said two ways.
    """
    reduced = _convert(
        odd_detector.trackingGeometry(), ODD_PIXEL_VOLUMES, beam_pipe_radius=24.0
    )
    shipped = syn.makeLayout(
        syn.readDetectorDescription(syn.dataPath("odd-description.json"))
    )

    reduced_cylinders, reduced_discs = _split(reduced)
    shipped_cylinders, shipped_discs = _split(shipped)

    # four barrel cylinders, and seven endcap layers per side that each resolve
    # into two discs staggered along z
    assert len(reduced_cylinders) == 4
    assert len(reduced_discs) == 2 * 14
    assert len(shipped_cylinders) == len(reduced_cylinders)
    assert len(shipped_discs) == len(reduced_discs)

    for converted, preset in zip(sorted(reduced_cylinders), sorted(shipped_cylinders)):
        assert converted == pytest.approx(preset, abs=0.1)
    # the rings themselves, not just the envelope around them
    for converted, preset in zip(sorted(reduced_discs), sorted(shipped_discs)):
        assert converted == pytest.approx(preset, abs=0.1)

    # and the endcaps really are mirror images of each other, which is what lets
    # the shipped file state one side
    positive = sorted(z for z, _, _ in reduced_discs if z > 0)
    negative = sorted(-z for z, _, _ in reduced_discs if z < 0)
    assert negative == pytest.approx(positive, abs=1e-3)

    # the beam pipe is there, passive, and in front of the innermost layer
    reduced_passive = [s for s in reduced.surfaces if s.passive]
    assert len(reduced_passive) == 1
    assert reduced_passive[0].refCoord < min(r for r, _, _ in reduced_cylinders)


def test_generic_detector_pixel_layout_matches_real_geometry(trk_geo):
    """Same check as for the ODD, on a detector that needs no external source."""
    layout = _convert(trk_geo, GENERIC_PIXEL_VOLUMES, beam_pipe_radius=19.0)
    cylinders, discs = _split(layout)
    description = syn.genericDetectorPixelDescription()
    preset_cylinders = _cylinders(description)
    preset_discs = [disc for disc, _ in _discs(description)]

    # unlike the ODD's, a generic endcap layer is one plane of one radial band,
    # so there is no ring structure for the reduction to find
    assert len(cylinders) == 4
    assert len(discs) == 2 * 7
    assert len(preset_discs) == 7
    assert all(len(d.rings) == 1 for d in preset_discs)

    for converted, preset in zip(
        sorted(r for r, _, _ in cylinders), sorted(c.radius for c in preset_cylinders)
    ):
        assert converted == pytest.approx(preset, abs=0.1)

    half_length = max(max(abs(lo), abs(hi)) for _, lo, hi in cylinders)
    assert half_length == pytest.approx(
        max(c.halfLengthZ for c in preset_cylinders), abs=0.1
    )

    positive = sorted(z for z, _, _ in discs if z > 0)
    for converted, preset in zip(positive, sorted(d.absZ for d in preset_discs)):
        assert converted == pytest.approx(preset, abs=0.1)

    assert min(lo for _, lo, _ in discs) == pytest.approx(
        min(r.rMin for d in preset_discs for r in d.rings), abs=0.1
    )
    assert max(hi for _, _, hi in discs) == pytest.approx(
        max(r.rMax for d in preset_discs for r in d.rings), abs=0.1
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


def test_reduction_names_subsystems_by_volume(trk_geo):
    """A geometry says which volume a surface is in and not much else, so that is
    what a subsystem is called unless the caller knows better."""
    description = _describe(trk_geo, GENERIC_PIXEL_VOLUMES)

    names = [subsystem.name for subsystem in description.subsystems]
    assert sorted(names) == [f"volume-{v}" for v in sorted(GENERIC_PIXEL_VOLUMES)]
    # the barrel volume holds the cylinders and the endcap volumes the discs
    barrel = next(s for s in description.subsystems if s.name == "volume-8")
    assert len(_cylinders(_one_subsystem(barrel))) == 4
    assert not barrel.endcaps


def test_caller_can_group_volumes_into_subsystems(trk_geo):
    """A caller that knows the detector calls its systems what they are, which is
    what selection and material then key onto."""
    options = _options(GENERIC_PIXEL_VOLUMES)
    options.setSubsystemName(lambda s: "generic-pixel")
    description = syn.makeDescriptionFromTrackingGeometry(
        trk_geo, acts.GeometryContext.dangerouslyDefaultConstruct(), options
    )

    assert [s.name for s in description.subsystems] == ["generic-pixel"]
    assert len(_cylinders(description)) == 4
    assert len(_discs(description)) == 2 * 7

    # and the layout says which subsystem each of its layers belongs to
    layout = syn.makeLayout(description)
    assert list(layout.subsystems) == ["generic-pixel"]
    assert all(layer.subsystem == 0 for layer in layout.layers)


def test_subsystems_can_be_selected_out_of_a_reduction(trk_geo):
    """The whole detector is reduced once and built a system at a time."""
    options = _options()
    options.setSubsystemName(
        lambda s: "pixel" if s.geometryId.volume in GENERIC_PIXEL_VOLUMES else "rest"
    )
    options.passiveBeamPipeRadius = 19.0
    whole = syn.makeDescriptionFromTrackingGeometry(
        trk_geo, acts.GeometryContext.dangerouslyDefaultConstruct(), options
    )
    assert sorted(s.name for s in whole.subsystems) == ["pixel", "rest"]

    pixels = syn.selectSubsystems(whole, ["pixel"])
    assert [s.name for s in pixels.subsystems] == ["pixel"]
    # the beam pipe belongs to no subsystem and stays, as does the containment of
    # the whole tracker
    assert len(pixels.passives) == len(whole.passives) == 1
    assert pixels.escapeRadius == whole.escapeRadius

    assert len(syn.makeLayout(pixels).surfaces) < len(syn.makeLayout(whole).surfaces)
    with pytest.raises(ValueError):
        syn.selectSubsystems(whole, ["pixels"])


def test_a_reduction_round_trips_through_a_file(trk_geo, tmp_path):
    """What the reduction found can be written out and read back, which is how the
    shipped descriptions were produced."""
    options = _options(GENERIC_PIXEL_VOLUMES, beam_pipe_radius=19.0)
    options.setSubsystemName(lambda s: "generic-pixel")
    description = syn.makeDescriptionFromTrackingGeometry(
        trk_geo, acts.GeometryContext.dangerouslyDefaultConstruct(), options
    )

    # where the layers are and what they are made of are two files
    decoration = syn.extractMaterial(description)
    assert len(decoration) > 0
    syn.stripMaterial(description)
    syn.writeDetectorDescription(str(tmp_path / "description.json"), description)
    syn.writeMaterialDecoration(str(tmp_path / "material.json"), decoration)

    bare = syn.readDetectorDescription(str(tmp_path / "description.json"))
    assert [s.name for s in bare.subsystems] == ["generic-pixel"]
    assert all(not c.material.bands for c in _cylinders(bare))

    syn.decorate(bare, syn.readMaterialDecoration(str(tmp_path / "material.json")))
    assert all(c.material.bands for c in _cylinders(bare))

    reread = syn.makeLayout(bare)
    assert len(reread.surfaces) == len(syn.makeLayout(description).surfaces)
    for surface, expected in zip(reread.surfaces, syn.makeLayout(description).surfaces):
        assert surface.refCoord == expected.refCoord


def test_shipped_configurations_are_files(tmp_path):
    """A fitted configuration is read from file, every field of it."""
    config = syn.readEventConfig(syn.dataPath("itk-ttbar-pu200.json"))
    assert config.generation.pileup == 200

    syn.writeEventConfig(str(tmp_path / "config.json"), config)
    again = syn.readEventConfig(str(tmp_path / "config.json"))
    assert again.generation.pileup == config.generation.pileup
    assert again.simulation.secondaries.electronRate == pytest.approx(
        config.simulation.secondaries.electronRate
    )
