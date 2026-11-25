# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import os
import subprocess
import sys

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


def test_material_validation_end_to_end(toy_detector, tmp_path):
    det, names = toy_detector
    datadir = str(tmp_path)

    # Fill material scan config.
    scan_cfg = detray.tests.MaterialScanConfig()
    scan_cfg.track_generator.phi_steps = 2
    scan_cfg.track_generator.eta_steps = 2
    scan_cfg.track_generator.uniform_eta = True
    scan_cfg.material_file = os.path.join(datadir, "material_scan")

    # Fill material validation config.
    validation_cfg = detray.tests.MaterialValidationConfig()
    validation_cfg.relative_error = 0.01
    validation_cfg.material_file = os.path.join(datadir, "navigation_material_trace")

    assert 0 == detray.tests.run_material_validation(
        det, names, scan_cfg, validation_cfg
    )

    scan_csv = os.path.join(datadir, "toy_detector_material_scan.csv")
    trace_csv = os.path.join(datadir, "toy_detector_navigation_material_trace_cpu.csv")
    assert os.path.getsize(scan_csv) > 0
    assert os.path.getsize(trace_csv) > 0
