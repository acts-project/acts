# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import os

import pytest

import detray.core
import detray.io
import detray.tests


def test_material_validation_end_to_end(detector, tmp_path):
    det, names = detector
    datadir = str(tmp_path)

    # Fill material scan config.
    scan_cfg = detray.tests.MaterialScanConfig()
    scan_cfg.trackGenerator.phiSteps = 2
    scan_cfg.trackGenerator.etaSteps = 2
    scan_cfg.trackGenerator.uniformEta = True
    scan_cfg.materialFile = os.path.join(datadir, "material_scan")

    # Fill material validation config.
    validation_cfg = detray.tests.MaterialValidationConfig()
    validation_cfg.relativeError = 0.01
    validation_cfg.materialFile = os.path.join(datadir, "navigation_material_trace")

    assert 0 == detray.tests.runMaterialValidation(det, names, scan_cfg, validation_cfg)

    scan_csv = os.path.join(datadir, "toy_detector_material_scan.csv")
    trace_csv = os.path.join(datadir, "toy_detector_navigation_material_trace_cpu.csv")
    assert os.path.getsize(scan_csv) > 0
    assert os.path.getsize(trace_csv) > 0


def test_detector_outlives_memory_resource(toy_detector_files):
    reader_config = (
        detray.io.DetectorReaderConfig()
        .addFile(toy_detector_files["geometry"])
        .addFile(toy_detector_files["material"])
    )
    memory_resource = detray.core.HostMemoryResource()
    det, names = detray.io.readDetector(memory_resource, reader_config)

    del memory_resource
    gc.collect()

    assert len(det.volumes) == 22
    assert det.name(names) == "toy_detector"


def test_detector_volume_container_bounds(detector):
    det, _ = detector
    n = len(det.volumes)
    assert n == 22

    det.volumes[0]
    det.volumes[-1]  # Wraps to last element.

    with pytest.raises(IndexError):
        det.volumes[n]
    with pytest.raises(IndexError):
        det.volumes[-(n + 1)]
