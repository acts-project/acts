# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import json
import weakref

import pytest

import detray
import detray.io


def test_common_type_and_metadata_override(toy_detector_files, tmp_path):
    cfg = detray.io.DetectorReaderConfig().addFile(toy_detector_files["geometry"])
    mr = detray.HostMemoryResource()
    toy, names = detray.io.readDetector(mr, cfg)
    default, _ = detray.io.readDetector(mr, cfg, metadataId="default")
    assert type(toy) is type(default) is detray.Detector
    assert toy.metadataId == "toy"
    assert default.metadataId == "default"
    assert len(toy.surfaces) == len(default.surfaces)
    assert toy.surfaces[0].identifier == default.surfaces[0].identifier
    assert len(list(toy.volumes)) == len(toy.volumes)

    writer = detray.io.DetectorWriterConfig()
    writer.path = str(tmp_path)
    detray.io.writeDetector(toy, names, writer)
    output = tmp_path / "toy_detector_geometry.json"
    data = json.loads(output.read_text())
    assert data["header"]["common"]["metadata"] == "toy"
    data["header"]["common"]["metadata"] = "unavailable_family"
    output.write_text(json.dumps(data))
    bad_cfg = detray.io.DetectorReaderConfig().addFile(str(output))
    with pytest.raises(ValueError, match="metadata ID"):
        detray.io.readDetector(mr, bad_cfg)


def test_inspection_collection_retains_allocator(toy_detector_files):
    cfg = detray.io.DetectorReaderConfig().addFile(toy_detector_files["geometry"])
    mr = detray.HostMemoryResource()
    mr_ref = weakref.ref(mr)
    detector, _ = detray.io.readDetector(mr, cfg)
    surfaces = detector.surfaces
    record = surfaces[0]
    del detector, mr
    gc.collect()
    assert mr_ref() is not None
    assert surfaces[0].identifier == record.identifier
    del surfaces
    gc.collect()
    assert mr_ref() is None
    assert record.index >= 0  # Inspection records are independent snapshots.
