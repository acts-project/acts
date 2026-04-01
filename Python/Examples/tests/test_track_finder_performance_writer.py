# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Tests for the Python-only TrackFinderPerformanceWriter."""

from pathlib import Path

import pytest

import acts
import acts.examples

u = acts.UnitConstants


# ---------------------------------------------------------------------------
# Smoke tests: class and config accessibility
# ---------------------------------------------------------------------------


def test_class_and_config_exist():
    assert hasattr(acts.examples, "PythonTrackFinderPerformanceWriter")
    assert hasattr(acts.examples.PythonTrackFinderPerformanceWriter, "Config")


def test_config_fields():
    cfg = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg.inputTracks = "tracks"
    cfg.inputParticles = "particles"
    cfg.inputTrackParticleMatching = "track_particle_matching"
    cfg.inputParticleTrackMatching = "particle_track_matching"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    assert cfg.inputTracks == "tracks"
    assert cfg.inputParticles == "particles"


def test_instantiation_raises_on_empty_config():
    """Writer should raise on empty required string fields."""
    cfg = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg.inputTracks = "tracks"
    # Missing all other required fields — should raise
    with pytest.raises(Exception):
        acts.examples.PythonTrackFinderPerformanceWriter(cfg, acts.logging.WARNING)


def test_instantiation_with_valid_config():
    cfg = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg.inputTracks = "tracks"
    cfg.inputParticles = "particles"
    cfg.inputTrackParticleMatching = "track_particle_matching"
    cfg.inputParticleTrackMatching = "particle_track_matching"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    writer = acts.examples.PythonTrackFinderPerformanceWriter(cfg, acts.logging.WARNING)
    assert writer is not None
    assert writer.config.inputTracks == "tracks"


# ---------------------------------------------------------------------------
# Integration test: run a minimal tracking chain and inspect histograms()
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_histograms_after_truth_kalman_run(tmp_path):
    """Run truth-seeded Kalman tracking, add the writer, check histograms()."""
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        MomentumConfig,
        PhiConfig,
        addFatras,
        addDigitization,
        addDigiParticleSelection,
        ParticleSelectorConfig,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TrackSmearingSigmas,
        addKalmanTracks,
    )

    srcdir = Path(__file__).parent.parent.parent.parent

    # Use Generic detector (no DD4hep required)
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=5, numThreads=1, logLevel=acts.logging.WARNING)

    addParticleGun(
        s,
        ParticleConfig(num=4, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        PhiConfig(0.0, 360.0 * u.degree),
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
        digiConfigFile=srcdir / "Examples/Configs/generic-digi-smearing-config.json",
        rnd=rnd,
    )

    addDigiParticleSelection(s, ParticleSelectorConfig())

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        rnd=rnd,
        trackSmearingSigmas=TrackSmearingSigmas(),
        initialSigmas=None,
        initialSigmaQoverPt=None,
        initialSigmaPtRel=None,
        initialVarInflation=None,
        particleHypothesis=None,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
    )

    # Configure the writer under test
    cfg = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg.inputTracks = "tracks"
    cfg.inputParticles = "particles"
    cfg.inputTrackParticleMatching = "track_particle_matching"
    cfg.inputParticleTrackMatching = "particle_track_matching"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"

    writer = acts.examples.PythonTrackFinderPerformanceWriter(cfg, acts.logging.WARNING)
    s.addWriter(writer)

    s.run()

    hists = writer.histograms()

    # Basic sanity: dict must be non-empty
    assert len(hists) > 0, "histograms() returned an empty dict"

    # Check for a known set of expected histogram names
    expected_keys = [
        "trackeff_vs_eta",
        "fakeRatio_vs_eta",
        "nMeasurements_vs_eta",
    ]
    for key in expected_keys:
        assert (
            key in hists
        ), f"Expected histogram key '{key}' not found in {list(hists.keys())}"

    # Check types: efficiency histograms are Efficiency1, profile are ProfileHistogram1
    assert isinstance(hists["trackeff_vs_eta"], acts.Efficiency1)
    assert isinstance(hists["nMeasurements_vs_eta"], acts.ProfileHistogram1)

    # Check that numpy arrays have sensible shapes (non-zero number of bins)
    eff = hists["trackeff_vs_eta"]
    accepted = eff.accepted().values()
    total = eff.total().values()
    assert accepted.ndim == 1
    assert accepted.shape == total.shape
    assert accepted.shape[0] > 0

    prof = hists["nMeasurements_vs_eta"]
    means = prof.histogram().means()
    assert means.ndim == 1
    assert means.shape[0] > 0
