# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import filecmp
import os

import pytest

import detray.core
import detray.io
import detray.tests


def test_material_scan_config():
    cfg = detray.tests.MaterialScanConfig()
    cfg.track_generator.eta_steps = 9
    assert cfg.track_generator.eta_steps == 9


def test_material_validation_config():
    cfg = detray.tests.MaterialValidationConfig()
    cfg.material_file = "navigation_output"
    cfg.n_tracks = 42
    cfg.propagation.stepping.path_limit = 1234.0

    assert cfg.material_file == "navigation_output"
    assert cfg.n_tracks == 42
    assert cfg.propagation.stepping.path_limit == pytest.approx(1234.0, rel=1e-6)


def _material_validation_configs(datadir):
    """Material scan and validation configs that write into @p datadir"""
    scan_cfg = detray.tests.MaterialScanConfig()
    scan_cfg.track_generator.phi_steps = 2
    scan_cfg.track_generator.eta_steps = 2
    scan_cfg.track_generator.uniform_eta = True
    scan_cfg.material_file = os.path.join(datadir, "material_scan")

    validation_cfg = detray.tests.MaterialValidationConfig()
    validation_cfg.relative_error = 0.01
    validation_cfg.material_file = os.path.join(datadir, "navigation_material_trace")

    return scan_cfg, validation_cfg


def test_material_validation_end_to_end(toy_detector, tmp_path):
    det, names = toy_detector
    datadir = str(tmp_path)

    scan_cfg, validation_cfg = _material_validation_configs(datadir)

    result = detray.tests.run_material_validation(det, names, scan_cfg, validation_cfg)
    assert result.passed
    assert result.n_tracks == 4
    assert result.n_failed == 0

    scan_csv = os.path.join(datadir, "toy_detector_material_scan.csv")
    trace_csv = os.path.join(datadir, "toy_detector_navigation_material_trace_cpu.csv")
    assert os.path.getsize(scan_csv) > 0
    assert os.path.getsize(trace_csv) > 0


def test_material_validation_truth_from_file(toy_detector, tmp_path):
    det, names = toy_detector
    datadir = str(tmp_path)

    # Run the material scan once and keep the navigation trace
    scan_cfg, validation_cfg = _material_validation_configs(datadir)
    first = detray.tests.run_material_validation(det, names, scan_cfg, validation_cfg)

    trace_csv = os.path.join(datadir, "toy_detector_navigation_material_trace_cpu.csv")
    first_trace_csv = os.path.join(datadir, "first_trace.csv")
    os.rename(trace_csv, first_trace_csv)

    # Read the truth material back from file (in the same process)
    scan_cfg.truth_from_file = True
    second = detray.tests.run_material_validation(
        det, names, scan_cfg, validation_cfg
    )

    assert second.passed
    assert second.n_tracks == first.n_tracks
    assert filecmp.cmp(first_trace_csv, trace_csv, shallow=False)

    # The truth data must match the configured tracks
    scan_cfg.track_generator.eta_steps = 3
    with pytest.raises(ValueError):
        detray.tests.run_material_validation(det, names, scan_cfg, validation_cfg)
