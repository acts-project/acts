import os
import tempfile

import detray.core
import detray.examples


def test_generate_and_read_detector():
    # Use the toy generator to produce detector files to read back.
    out_dir = tempfile.mkdtemp() + os.sep
    detray.examples.generateToyDetector(outputDir=out_dir)

    geometry_file = os.path.join(out_dir, "toy_detector_geometry.json")
    assert os.path.exists(geometry_file)

    reader_config = detray.core.DetectorReaderConfig().addFile(geometry_file)

    detector, names = detray.core.readDetector(reader_config)
    assert isinstance(names, detray.core.NameMap)
    assert len(detector.volumes) == 22
    assert len(detector.surfaces) == 3230
