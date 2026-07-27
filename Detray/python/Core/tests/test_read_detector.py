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

    detector, names = detray.core.readDetector(geometry_file)
    assert isinstance(names, detray.core.NameMap)
    assert detector.n_volumes() == 22
    assert detector.n_surfaces() == 3230
